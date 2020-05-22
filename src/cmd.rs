use std::path::Path;

use failure::Error;
use log::info;
use needletail::parse_sequence_path;
use sourmash::cmd::ComputeParameters;
use sourmash::index::storage::ToWriter;
use sourmash::signature::Signature;

pub fn compute<P: AsRef<Path>>(
    filenames: Vec<P>,
    params: &ComputeParameters,
) -> Result<Vec<Signature>, Error> {
    if params.merge.is_some() {
        // make one signature for all files
        let mut n = 0;
        let mut total_seq = 0;

        let mut sig = Signature::from_params(&params);
        sig.set_filename(filenames.iter().last().unwrap().as_ref().to_str().unwrap());
        sig.set_name(&params.merge.clone().unwrap());

        filenames.iter().for_each(|filename| {
            // consume & calculate signatures
            info!(
                "... reading sequences from {}",
                filename.as_ref().to_str().unwrap()
            );

            parse_sequence_path(
                filename,
                |_| {},
                |record| {
                    if params.input_is_protein {
                        sig.add_protein(&record.seq).expect("Error adding sequence");
                    } else {
                        sig.add_sequence(&record.seq, !params.check_sequence)
                            .expect("Error adding sequence");
                    }
                    n += 1;
                },
            )
            .unwrap();

            total_seq += n + 1;
        });

        info!(
            "calculated {} signatures for {} sequences taken from {} files",
            sig.size(),
            total_seq,
            filenames.len()
        );

        let mut output = niffler::to_path(
            params.output.as_ref().unwrap(),
            niffler::compression::Format::No,
            niffler::compression::Level::One,
        )?;
        sig.to_writer(&mut output).unwrap();

        return Ok(vec![sig]);
    }

    // Not merging
    let mut siglist: Vec<Signature> = vec![];

    filenames.iter().for_each(|filename| {
        let sigfile = params
            .output
            .clone()
            .unwrap_or(format!("{}.sig", filename.as_ref().to_str().unwrap()));

        if params.singleton {
            parse_sequence_path(
                filename,
                |_| {},
                |record| {
                    let mut sig = Signature::from_params(&params);
                    sig.set_name(&String::from_utf8(record.id.into_owned()).unwrap());
                    sig.set_filename(&filename.as_ref().to_str().unwrap());

                    if params.input_is_protein {
                        sig.add_protein(&record.seq).expect("Error adding sequence");
                    } else {
                        sig.add_sequence(&record.seq, !params.check_sequence)
                            .expect("Error adding sequence");
                    }
                    siglist.push(sig);
                },
            )
            .unwrap();
        } else if params.input_is_10x {
            // TODO: implement 10x parsing
            unimplemented!()
        } else {
            // make minhashes for the whole file

            let fname: String = filename.as_ref().to_str().unwrap().into();
            info!("... reading sequences from {}", &fname);

            let mut sig = Signature::from_params(&params);
            sig.set_name(&fname);
            sig.set_filename(&fname);

            let mut name = None;

            parse_sequence_path(
                filename,
                |_| {},
                |record| {
                    if params.input_is_protein {
                        sig.add_protein(&record.seq).expect("Error adding sequence");
                    } else {
                        sig.add_sequence(&record.seq, !params.check_sequence)
                            .expect("Error adding sequence");
                    }

                    if params.name_from_first && name.is_none() {
                        name = Some(String::from_utf8(record.id.into_owned()).unwrap())
                    };
                },
            )
            .unwrap();

            if let Some(n) = name {
                sig.set_name(&n)
            };

            if params.output.is_none() {
                let mut output = niffler::to_path(
                    &sigfile,
                    niffler::compression::Format::No,
                    niffler::compression::Level::One,
                )
                .expect("Error creating output file");
                sig.to_writer(&mut output).unwrap();
                siglist = vec![sig];
            } else {
                siglist.push(sig);
            }
        }

        if let Some(ref output_name) = params.output {
            let mut output = niffler::to_path(
                &output_name,
                niffler::compression::Format::No,
                niffler::compression::Level::One,
            )
            .expect("Error creating output file");
            serde_json::to_writer(&mut output, &siglist).unwrap();
        }
    });

    Ok(siglist)
}

pub fn compare(
    signatures: &[Signature],
    params: &CompareParameters,
) -> Result<Vec<Vec<f64>>, Error> {
    unimplemented!();
}

pub struct CompareParameters {
    pub ksize: u32,
    pub modhash: bool,

    pub dna: bool,
    pub dayhoff: bool,
    pub hp: bool,
    pub protein: bool,

    pub output: Option<String>,
    pub csv: Option<String>,

    pub ignore_abundance: bool,
    pub traverse_directory: bool,

    pub processes: usize,
}

impl Default for CompareParameters {
    fn default() -> Self {
        CompareParameters {
            ksize: 31,
            modhash: false,
            dna: true,
            dayhoff: false,
            hp: false,
            protein: false,
            output: None,
            csv: None,
            ignore_abundance: false,
            traverse_directory: false,
            processes: 1,
        }
    }
}
