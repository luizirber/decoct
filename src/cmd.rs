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
use log::info;
use needletail::parse_sequence_path;
use niffler::{get_output, CompressionFormat};
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
        sig.filename = Some(
            filenames
                .iter()
                .last()
                .unwrap()
                .as_ref()
                .to_str()
                .unwrap()
                .into(),
        );

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
                    n += 1;
                },
            )
            .unwrap();

            total_seq += n + 1;
        });

        info!(
            "calculated {} signatures for {} sequences taken from {} files",
            sig.signatures.len(),
            total_seq,
            filenames.len()
        );

        let mut output =
            get_output(params.output.as_ref().unwrap(), CompressionFormat::No).unwrap();
        sig.to_writer(&mut output).unwrap();

        Ok(vec![sig])
    } else {
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
                        let fname = Some(filename.as_ref().to_str().unwrap().into());
                        let mut sig = Signature::from_params(&params);
                        sig.name = Some(String::from_utf8(record.id.into_owned()).unwrap());
                        sig.filename = fname.clone();
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

                        sig.add_sequence(&seq, params.check_sequence)
                            .expect("Error adding sequence");
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
                sig.name = Some(fname.clone());
                sig.filename = Some(fname.clone());

                let mut name = None;

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

                        sig.add_sequence(&seq, params.check_sequence)
                            .expect("Error adding sequence");

                        if params.name_from_first && name.is_none() {
                            name = Some(String::from_utf8(record.id.into_owned()).unwrap())
                        };
                    },
                )
                .unwrap();

                dbg!(&name);
                if let Some(n) = name {
                    sig.name = Some(n)
                };

                if params.output.is_none() {
                    let mut output = get_output(&sigfile, CompressionFormat::No).unwrap();
                    sig.to_writer(&mut output).unwrap();
                    siglist = vec![sig];
                } else {
                    siglist.push(sig);
                }
            }

            if let Some(ref output_name) = params.output {
                let mut output = get_output(&output_name, CompressionFormat::No).unwrap();
                serde_json::to_writer(&mut output, &siglist).unwrap();
            }
        });

        Ok(siglist)
    }
}
