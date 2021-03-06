use std::fs;
use std::io::{Read, Write};
use std::process::Command;

use assert_cmd::prelude::*;
use predicates::str::contains;
use sourmash::signature::Signature;
use tempfile::TempDir;

#[test]
fn search() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("decoct")?;

    cmd.arg("search")
        .arg("tests/data/SRR2060939_1.sig")
        .arg("tests/data/v5.sbt.json")
        .assert()
        .success()
        .stdout(contains("SRR2060939_1.fastq.gz"))
        .stdout(contains("SRR2060939_2.fastq.gz"))
        .stdout(contains("SRR2255622_1.fastq.gz"));

    Ok(())
}

#[test]
#[ignore]
fn search_only_leaves() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("decoct")?;

    cmd.arg("search")
        .arg("tests/data/SRR2060939_1.sig")
        .arg("tests/data/leaves.sbt.json")
        .assert()
        .success()
        .stdout(contains("SRR2060939_1.fastq.gz"))
        .stdout(contains("SRR2060939_2.fastq.gz"))
        .stdout(contains("SRR2255622_1.fastq.gz"));

    Ok(())
}

#[test]
#[ignore]
#[cfg(unix)]
fn compute_index_and_search() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy("tests/data/short.fa", tmp_dir.path().join("short.fa"))?;
    fs::copy("tests/data/short2.fa", tmp_dir.path().join("short2.fa"))?;

    assert!(tmp_dir.path().join("short.fa").exists());
    assert!(tmp_dir.path().join("short2.fa").exists());

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .args(&["short.fa", "short2.fa"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("short.fa.sig").exists());
    assert!(tmp_dir.path().join("short2.fa.sig").exists());

    let mut cmd = Command::new("sourmash");
    //let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("index")
        .args(&["-k", "31"])
        //.args(&["-o", "zzz.sbt.json"])
        .arg("zzz.sbt.json")
        .args(&["short.fa.sig", "short2.fa.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("zzz.sbt.json").exists());

    let cmds = vec![Command::new("sourmash"), Command::cargo_bin("decoct")?];

    for mut cmd in cmds {
        cmd.arg("search")
            .args(&["-k", "31"])
            .arg("short.fa.sig")
            .arg("zzz.sbt.json")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("short.fa"))
            .stdout(contains("short2.fa"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn index_and_search() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy("tests/data/SRR2060939_1.sig", tmp_dir.path().join("1.sig"))?;
    fs::copy("tests/data/SRR2060939_2.sig", tmp_dir.path().join("2.sig"))?;

    assert!(tmp_dir.path().join("1.sig").exists());
    assert!(tmp_dir.path().join("2.sig").exists());

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("index")
        .args(&["-k", "31"])
        .args(&["-o", "zzz.sbt.json"])
        .args(&["1.sig", "2.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("zzz.sbt.json").exists());

    let cmds = vec![Command::new("sourmash"), Command::cargo_bin("decoct")?];

    for mut cmd in cmds {
        cmd.arg("search")
            .args(&["-k", "31"])
            .arg("1.sig")
            .arg("zzz.sbt.json")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("2 matches:"))
            .stdout(contains("SRR2060939_1.fastq.gz"))
            .stdout(contains("SRR2060939_2.fastq.gz"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/data/ecoli.genes.fna",
        tmp_dir.path().join("ecoli.fna"),
    )?;

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna")
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli.fna.sig").exists());

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_sourmash.fna.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_sourmash.fna.sig").exists());

    for k in &["21", "31", "51"] {
        let mut cmd = Command::new("sourmash");
        cmd.arg("compare")
            .args(&["-k", k])
            .arg("ecoli.fna.sig")
            .arg("ecoli_sourmash.fna.sig")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("min similarity in matrix: 1.000"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute_compressed_input() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    let file = tmp_dir.path().join("ecoli.fna.gz");
    {
        let mut writer = niffler::get_writer(
            Box::new(std::fs::File::create(file)?),
            niffler::compression::Format::Gzip,
            niffler::compression::Level::Nine,
        )?;
        let mut input = std::fs::File::open("tests/data/ecoli.genes.fna")?;
        let mut buffer = Vec::new();
        input.read_to_end(&mut buffer)?;

        writer.write_all(&buffer)?;
    }

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna.gz")
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli.fna.gz.sig").exists());

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .arg("ecoli.fna.gz")
        .args(&["-o", "ecoli_sourmash.fna.gz.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_sourmash.fna.gz.sig").exists());

    for k in &["21", "31", "51"] {
        let mut cmd = Command::new("sourmash");
        cmd.arg("compare")
            .args(&["-k", k])
            .arg("ecoli.fna.gz.sig")
            .arg("ecoli_sourmash.fna.gz.sig")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("min similarity in matrix: 1.000"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute_output() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/data/ecoli.genes.fna",
        tmp_dir.path().join("ecoli.fna"),
    )?;

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_decoct.fna.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_decoct.fna.sig").exists());

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_sourmash.fna.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_sourmash.fna.sig").exists());

    for k in &["21", "31", "51"] {
        let mut cmd = Command::new("sourmash");
        cmd.arg("compare")
            .args(&["-k", k])
            .arg("ecoli_decoct.fna.sig")
            .arg("ecoli_sourmash.fna.sig")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("min similarity in matrix: 1.000"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute_merge() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/data/ecoli.genes.fna",
        tmp_dir.path().join("ecoli.fna"),
    )?;

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_decoct.fna.sig"])
        .args(&["--merge", "decoct"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_decoct.fna.sig").exists());
    // TODO: parse sig and assert name == decoct

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_sourmash.fna.sig"])
        .args(&["--merge", "sourmash"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_sourmash.fna.sig").exists());
    // TODO: parse sig and assert name == sourmash

    for k in &["21", "31", "51"] {
        let mut cmd = Command::new("sourmash");
        cmd.arg("compare")
            .args(&["-k", k])
            .arg("ecoli_decoct.fna.sig")
            .arg("ecoli_sourmash.fna.sig")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("0-decoct"))
            .stdout(contains("1-sourmash"))
            .stdout(contains("min similarity in matrix: 1.000"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute_merge_no_output() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/data/ecoli.genes.fna",
        tmp_dir.path().join("ecoli.fna"),
    )?;

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["--merge", "decoct"])
        .current_dir(&tmp_dir)
        .assert()
        .failure()
        .code(255);

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute_singleton() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/data/ecoli.genes.fna",
        tmp_dir.path().join("ecoli.fna"),
    )?;

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_decoct.fna.sig"])
        .arg("--singleton")
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_decoct.fna.sig").exists());

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .arg("ecoli.fna")
        .args(&["-o", "ecoli_sourmash.fna.sig"])
        .arg("--singleton")
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli_sourmash.fna.sig").exists());

    for k in &["21", "31", "51"] {
        let mut cmd = Command::new("sourmash");
        cmd.arg("compare")
            .args(&["-k", k])
            .arg("ecoli_decoct.fna.sig")
            .arg("ecoli_sourmash.fna.sig")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("0-gi|556503834:33...\t[1. 0. 1. 0.]"))
            .stdout(contains("1-gi|556503834:28...\t[0. 1. 0. 1.]"))
            .stdout(contains("2-gi|556503834:33...\t[1. 0. 1. 0.]"))
            .stdout(contains("3-gi|556503834:28...\t[0. 1. 0. 1.]"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn compute_name_from_first() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/data/ecoli.genes.fna",
        tmp_dir.path().join("ecoli.fna"),
    )?;

    let mut cmd = Command::cargo_bin("decoct")?;
    cmd.arg("compute")
        .arg("ecoli.fna")
        .arg("--name-from-first")
        .args(&["-k", "31"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("ecoli.fna.sig").exists());

    let sigs = Signature::from_path(tmp_dir.path().join("ecoli.fna.sig"))?;
    assert_eq!(sigs.len(), 1);

    let sig = sigs.first().unwrap();
    dbg!(sig);
    assert_eq!(
        sig.name(),
        "gi|556503834:337-2799 Escherichia coli str. K-12 substr. MG1655, complete genome"
    );

    Ok(())
}
