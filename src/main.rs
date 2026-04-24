mod hashing;
mod pipeline;
mod sequence;

use crate::hashing::SequenceHasher;
use crate::pipeline::{create_fasta_reader, exit_after_flush, flush_writer_or_exit, pipeline};
use crate::sequence::SequenceProcessor;
use clap::{
    builder::styling::{AnsiColor, Style, Styles},
    CommandFactory, Parser,
};
use clio::Input;
use std::io::{self, IsTerminal, LineWriter};
use std::num::NonZeroU8;
use std::process;

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Cyan.on_default().bold())
    .usage(AnsiColor::Yellow.on_default().bold())
    .literal(AnsiColor::Yellow.on_default().bold())
    .placeholder(Style::new().dimmed());

/// Compute hash digests for sequences in a FASTA file
#[derive(Parser)]
#[command(version, about, max_term_width = 79, styles = STYLES)]
struct Cli {
    /// Input file(s). Use '-' for stdin
    #[arg(default_value = "-")]
    input: Vec<Input>,

    /// Print sequences in a third column
    #[arg(short = 's', long, default_value = "false", help_heading = "Output")]
    print_sequence: bool,

    /// Instead of hashing the entire sequence at once, hash each k-mer
    /// individually and then combine the resulting hashes
    #[arg(short = 'm', long, default_value = "false", help_heading = "Hashing")]
    multi_kmer_hashing: bool,

    /// Replace ntHash with xxHash for hashing k-mers. Works only with
    /// --multi-kmer-hashing
    #[arg(
        short = 'x',
        long,
        requires = "multi_kmer_hashing",
        default_value = "false",
        help_heading = "Hashing"
    )]
    xxhash: bool,

    /// Size of the k-mers to hash when using --multi-kmer-hashing
    #[arg(
        short = 'k',
        long = "kmer-size",
        requires = "multi_kmer_hashing",
        default_value = "31",
        help_heading = "Hashing"
    )]
    k: NonZeroU8,

    /// Make hashing robust to circular permutations via deterministic rotation
    /// to the lexicographically minimal sequence
    #[arg(
        short = 'r',
        long,
        default_value = "false",
        conflicts_with = "circular_kmers",
        help_heading = "Circular sequences"
    )]
    circular_rotation: bool,

    /// Make hashing robust to circular permutations via addition of the
    /// k-mers that wrap around the end of the sequence. Works only with
    /// --multi-kmer-hashing
    #[arg(
        short = 'w',
        long,
        default_value = "false",
        requires = "multi_kmer_hashing",
        conflicts_with_all = ["circular_rotation"],
        help_heading = "Circular sequences"
    )]
    circular_kmers: bool,
}

fn main() {
    let cli = Cli::parse();
    let input_count = cli.input.len();

    let sequence_processor =
        SequenceProcessor::new(cli.circular_rotation, cli.circular_kmers, cli.k);
    let hasher = SequenceHasher::new(cli.multi_kmer_hashing, cli.xxhash, cli.k);

    // If it's an interactive session with no data piped to stdin and files provided,
    // show help and exit
    if input_count == 1 && cli.input[0].is_std() && io::stdin().is_terminal() {
        Cli::command().print_help().unwrap();
        process::exit(0);
    }

    let stdout = io::stdout();
    let mut writer = LineWriter::new(stdout.lock());

    for input in cli.input {
        let input_display = input.to_string();

        let reader = match create_fasta_reader(input) {
            Ok(reader) => reader,
            Err(error_msg) => {
                eprintln!(
                    "Error: failed to create reader for {}: {}",
                    input_display, error_msg
                );
                exit_after_flush(&mut writer, 1);
            }
        };
        pipeline(
            reader,
            &mut writer,
            &hasher,
            &sequence_processor,
            cli.print_sequence,
        );
    }

    flush_writer_or_exit(&mut writer);
}
