// Note for compute: Use something like glium::DrawParameters to deal with
// parameters complexity...
// https://docs.rs/glium/0.26.0-alpha5/glium/draw_parameters/struct.DrawParameters.html
//
// let params = glium::DrawParameters {
//     line_width: Some(0.02),
//     point_size: Some(0.02),
//     .. Default::default()
// };
//
// target.draw(..., &params).unwrap();

use std::path::Path;

use failure::Error;
use needletail::parse_sequence_path;
use niffler::{get_output, CompressionFormat};
use sourmash::index::storage::ToWriter;
use sourmash::signature::Signature;
use sourmash::sketch::minhash::{max_hash_for_scaled, HashFunctions, KmerMinHash};
use sourmash::sketch::Sketch;

pub fn compute<P: AsRef<Path>>(filenames: Vec<P>, params: &ComputeParameters) -> Result<(), Error> {
    let template = build_template(&params);

    filenames.iter().for_each(|filename| {
        let sigfile = format!("{}.sig", filename.as_ref().to_str().unwrap());

        if params.singleton {
            // TODO: implement one sig per record
            unimplemented!()
        } else if params.input_is_10x {
            // TODO: implement 10x parsing
            unimplemented!()
        } else {
            // make minhashes for the whole file
            let mut sig = Signature::builder()
                .hash_function("0.murmur64")
                .name(Some(filename.as_ref().to_str().unwrap().into()))
                .signatures(template.clone())
                .build();

            parse_sequence_path(
                filename,
                |_| {},
                |record| {
                    // if there is anything other than ACGT in sequence,
                    // it is replaced with A.
                    // This matches khmer and screed behavior
                    //
                    // NOTE: sourmash is different! It uses the force flag to drop
                    // k-mers that are not ACGT
                    /*
                     for record in records {
                         if args.input_is_protein {
                             sig.add_protein(record.seq)
                         } else {
                             sig.add_sequence(record.seq, check_sequence)
                         }
                     }
                    */
                    let seq: Vec<u8> = record
                        .seq
                        .iter()
                        .map(|&x| match x as char {
                            'A' | 'C' | 'G' | 'T' => x,
                            'a' | 'c' | 'g' | 't' => x.to_ascii_uppercase(),
                            _ => b'A',
                        })
                        .collect();

                    sig.add_sequence(&seq, false)
                        .expect("Error adding sequence");
                },
            )
            .unwrap();

            let mut output = get_output(&sigfile, CompressionFormat::No).unwrap();
            sig.to_writer(&mut output).unwrap();
        }
    });
    Ok(())
}

pub struct ComputeParameters {
    pub ksizes: Vec<u32>,
    pub check_sequence: bool,
    pub dna: bool,
    pub dayhoff: bool,
    pub hp: bool,
    pub singleton: bool,
    pub count_valid_reads: usize,
    pub barcodes_file: Option<String>, // TODO: check
    pub line_count: usize,
    pub rename_10x_barcodes: Option<bool>,    // TODO: check
    pub write_barcode_meta_csv: Option<bool>, // TODO: check
    pub save_fastas: Option<bool>,            // TODO: check
    pub scaled: u64,
    pub force: bool,
    pub output: Option<String>, // TODO: check
    pub num_hashes: u32,
    pub protein: bool,
    pub name_from_first: bool,
    pub seed: usize,
    pub input_is_protein: bool,
    pub merge: Option<bool>, // TODO: check
    pub track_abundance: bool,
    pub randomize: bool,
    pub license: String,
    pub input_is_10x: bool,
    pub processes: usize,
}

impl Default for ComputeParameters {
    fn default() -> Self {
        ComputeParameters {
            ksizes: vec![21, 31, 51],
            check_sequence: false,
            dna: true,
            dayhoff: false,
            hp: false,
            singleton: false,
            count_valid_reads: 0,
            barcodes_file: None,
            line_count: 1500,
            rename_10x_barcodes: None,
            write_barcode_meta_csv: None,
            save_fastas: None,
            scaled: 10000,
            force: false,
            output: None,
            num_hashes: 500,
            protein: false,
            name_from_first: false,
            seed: 42,
            input_is_protein: false,
            merge: None,
            track_abundance: false,
            randomize: false,
            license: "CC0".into(),
            input_is_10x: false,
            processes: 2,
        }
    }
}

fn build_template(params: &ComputeParameters) -> Vec<Sketch> {
    let max_hash = max_hash_for_scaled(params.scaled).unwrap_or(0);

    params
        .ksizes
        .iter()
        .flat_map(|k| {
            let mut ksigs = vec![];

            if params.protein {
                ksigs.push(Sketch::MinHash(
                    KmerMinHash::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_protein)
                        .max_hash(max_hash)
                        .abunds(if params.track_abundance {
                            Some(vec![])
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            if params.dayhoff {
                ksigs.push(Sketch::MinHash(
                    KmerMinHash::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_dayhoff)
                        .max_hash(max_hash)
                        .abunds(if params.track_abundance {
                            Some(vec![])
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            if params.hp {
                ksigs.push(Sketch::MinHash(
                    KmerMinHash::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_hp)
                        .max_hash(max_hash)
                        .abunds(if params.track_abundance {
                            Some(vec![])
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            if params.dna {
                ksigs.push(Sketch::MinHash(
                    KmerMinHash::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_DNA)
                        .max_hash(max_hash)
                        .abunds(if params.track_abundance {
                            Some(vec![])
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            ksigs
        })
        .collect()
}