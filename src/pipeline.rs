use crate::hashing::SequenceHasher;
use crate::sequence::{get_record_accession, SequenceProcessor};
use needletail::{parse_fastx_file, parse_fastx_stdin};
use std::io::{self, Write};
use std::process;
use std::str;

use clio::Input;

pub fn pipeline(input: &Input, hasher: &SequenceHasher, sequence_processor: &SequenceProcessor) {
    let reader = match input.is_std() {
        true => parse_fastx_stdin(),
        false => {
            if input.is_empty().unwrap() {
                eprintln!("Error: the input file is empty");
                process::exit(1);
            }
            parse_fastx_file(input.path().to_path_buf())
        }
    };

    let mut reader = match reader {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Error: {}", e);
            process::exit(1);
        }
    };

    // Iterate over the sequence records
    while let Some(record) = reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(e) => {
                eprintln!("Error: {}", e);
                process::exit(1);
            }
        };

        // Get the accession of the record
        let record_header = record.id().to_vec();
        let accession = match get_record_accession(&record_header) {
            Some(acc) => acc,
            None => {
                eprintln!("Error: a record with an invalid header was found");
                process::exit(1);
            }
        };

        // Check if the record is shorter than the k-mer size
        if hasher.multi_kmer_hashing && record.num_bases() < hasher.k as usize {
            eprintln!(
                "Error: record {} is shorter than the k-mer size",
                str::from_utf8(accession).unwrap_or("'NA'")
            );
            process::exit(1);
        }

        // Normalize the sequence: capitalize all bases and remove newlines.
        // Handle circular sequences as follows:
        // - If `circular_kmers` is true, add k-mers formed by wrapping around the sequence.
        // - If `circular_rotation` is true, rotate to the lexicographically minimal form.
        let processed_seq = sequence_processor.process_sequence(record);

        // Compute the hash of the sequence
        let hash_seq = match hasher.compute_hash(&processed_seq) {
            Ok(hash) => hash,
            Err(e) => {
                eprintln!(
                    "Error: failed to compute the hash for record {}: {}",
                    str::from_utf8(accession).unwrap_or("'NA'"),
                    e
                );
                process::exit(1);
            }
        };

        // Print the record accession and the hash of the sequence
        let output = format!(
            "{}\t{}\n",
            str::from_utf8(accession).unwrap_or("'NA'"),
            hex::encode(hash_seq.to_be_bytes())
        );

        // Write to stdout and handle potential errors
        if let Err(e) = io::stdout().write_all(output.as_bytes()) {
            // Check if it's a broken pipe error
            if e.kind() == io::ErrorKind::BrokenPipe {
                // Exit gracefully
                std::process::exit(0);
            }
            // For other errors, you might want to handle them differently
            eprintln!("Error writing to stdout: {}", e);
            std::process::exit(1);
        }
    }
}
