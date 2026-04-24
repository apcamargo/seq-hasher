use crate::hashing::SequenceHasher;
use crate::sequence::{get_record_accession, SequenceProcessor};
use needletail::{parse_fastx_reader, parser::FastxReader};
use std::io::{self, Write};
use std::process;
use std::str;

use clio::Input;

const HEX_DIGITS: &[u8; 16] = b"0123456789abcdef";

pub fn create_fasta_reader(input: Input) -> Result<Box<dyn FastxReader>, String> {
    if input.can_seek() && input.is_empty() == Some(true) {
        return Err("the input file is empty".to_string());
    }
    parse_fastx_reader(input).map_err(|e| e.to_string())
}

pub fn pipeline(
    mut reader: Box<dyn FastxReader>,
    writer: &mut impl Write,
    hasher: &SequenceHasher,
    sequence_processor: &SequenceProcessor,
    print_sequence: bool,
) {
    // Iterate over the sequence records
    while let Some(record) = reader.next() {
        let record = match record {
            Ok(record) => record,
            Err(e) => {
                eprintln!("Error: {e}");
                exit_after_flush(writer, 1);
            }
        };

        // Get the accession of the record
        let accession = match get_record_accession(record.id()) {
            Some(acc) => acc,
            None => {
                eprintln!("Error: a record with an invalid header was found");
                exit_after_flush(writer, 1);
            }
        };

        // Check if the record is shorter than the k-mer size
        if hasher.multi_kmer_hashing && record.num_bases() < usize::from(hasher.k.get()) {
            eprintln!(
                "Error: record {} is shorter than the k-mer size",
                str::from_utf8(accession).unwrap_or("'NA'")
            );
            exit_after_flush(writer, 1);
        }

        // Normalize the sequence: capitalize all bases and remove newlines.
        // Handle circular sequences as follows:
        // - If `circular_kmers` is true, add k-mers formed by wrapping around the sequence.
        // - If `circular_rotation` is true, rotate to the lexicographically minimal form.
        let processed_seq = sequence_processor.process_sequence(&record);

        // Compute the hash of the sequence
        let hash_seq = match hasher.compute_hash(processed_seq.as_ref()) {
            Ok(hash) => hash,
            Err(e) => {
                eprintln!(
                    "Error: failed to compute the hash for record {}: {}",
                    str::from_utf8(accession).unwrap_or("'NA'"),
                    e
                );
                exit_after_flush(writer, 1);
            }
        };

        if let Err(error) = write_output_record(
            writer,
            accession,
            hash_seq,
            processed_seq.as_ref(),
            print_sequence,
        ) {
            handle_output_error(error);
        }
    }
}

pub fn flush_writer_or_exit(writer: &mut impl Write) {
    if let Err(error) = writer.flush() {
        handle_output_error(error);
    }
}

pub fn exit_after_flush(writer: &mut impl Write, code: i32) -> ! {
    flush_writer_or_exit(writer);
    process::exit(code);
}

fn write_output_record(
    writer: &mut impl Write,
    accession: &[u8],
    hash_seq: u128,
    processed_seq: &[u8],
    print_sequence: bool,
) -> io::Result<()> {
    writer.write_all(str::from_utf8(accession).unwrap_or("'NA'").as_bytes())?;
    writer.write_all(b"\t")?;
    write_hash_hex(writer, hash_seq)?;
    if print_sequence {
        writer.write_all(b"\t")?;
        writer.write_all(processed_seq)?;
    }
    writer.write_all(b"\n")
}

fn write_hash_hex(writer: &mut impl Write, hash_seq: u128) -> io::Result<()> {
    let mut hex = [0_u8; 32];
    for (idx, byte) in hash_seq.to_be_bytes().iter().enumerate() {
        hex[idx * 2] = HEX_DIGITS[usize::from(byte >> 4)];
        hex[idx * 2 + 1] = HEX_DIGITS[usize::from(byte & 0x0f)];
    }
    writer.write_all(&hex)
}

fn handle_output_error(error: io::Error) -> ! {
    if error.kind() == io::ErrorKind::BrokenPipe {
        process::exit(0);
    }
    eprintln!("Error writing to stdout: {error}");
    process::exit(1);
}
