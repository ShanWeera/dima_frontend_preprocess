mod utils;

use wasm_bindgen::prelude::*;
use seq_io::fasta::{Reader, Record, OwnedRecord};
use std::borrow::Borrow;
use serde::{Serialize, Deserialize};

#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
pub struct Sequence {
    index: usize,
    header: String,
    sequence: String,
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize)]
pub struct ProblemHeader {
    index: usize,
    header: String,
}

#[wasm_bindgen]
pub struct MSA {
    sequences: Option<Vec<Sequence>>,
}

#[wasm_bindgen]
impl MSA {
    #[wasm_bindgen(constructor)]
    pub fn new() -> MSA {
        MSA {
            sequences: None,
        }
    }

    pub fn set_seqs(&mut self, sequences: String) -> Result<(), JsValue> {
        let mut all_sequences: Vec<Sequence> = vec![];
        let mut reader = Reader::new(sequences.as_bytes());
        let mut counter = 0;

        while let Some(result) = reader.next() {
            counter += 1;

            if let Ok(record) = result {
                if let Ok((header_id, header_desc)) = record.id_desc() {
                    let corrected_header: String;

                    if let Some(desc) = header_desc {
                        corrected_header = [header_id, desc].join("");
                    } else {
                        corrected_header = header_id.to_string();
                    }

                    all_sequences.push(
                        Sequence {
                            index: counter,
                            sequence: String::from_utf8_lossy(record
                                .full_seq()
                                .borrow())
                                .parse()
                                .unwrap(),
                            header: corrected_header,
                        });
                }
            } else {
                return Err(
                    JsValue::from_str("Not a valid FASTA Multiple Sequence Alignment.")
                );
            };
        }

        let sequence_length = all_sequences[0].sequence.len();

        for sequence in all_sequences.iter() {
            if sequence.sequence.len() != sequence_length {
                return Err(
                    JsValue::from_str(&*format!("Not a valid MSA: {}", sequence.header))
                );
            };
        }

        self.sequences = Some(all_sequences);
        Ok(())
    }

    pub fn check_headers(&self, item_count: usize) -> JsValue {
         let problem_headers = self.
            sequences
            .as_ref()
            .unwrap()
            .iter()
            .filter(|sequence| {
                let header_parts = sequence.header.split("|");
                header_parts.count() != item_count
            })
            .map(|sequence| ProblemHeader {
                index: sequence.index,
                header: sequence.header.to_string(),
            })
            .collect::<Vec<ProblemHeader>>();

        JsValue::from_serde(&problem_headers).unwrap()
    }

    pub fn get_seq_count(&self) -> Result<usize, JsValue> {
        if let Some(sequences) = self.sequences.as_ref() {
            Ok(sequences.len())
        } else {
            Err(
                JsValue::from_str("No sequences set. Please upload a valid MSA.")
            )
        }
    }

    pub fn get_seqs(&self) -> String {
        self
            .sequences
            .as_ref()
            .unwrap()
            .iter()
            .map(|sequence| {
                let header = [">", sequence.header.as_str()].join("");
                [header, sequence.sequence.to_string()].join("\n")
            })
            .collect::<Vec<String>>()
            .join("\n")
    }
}
