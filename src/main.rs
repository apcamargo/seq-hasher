mod hashing;
mod pipeline;
mod sequence;

use crate::hashing::SequenceHasher;
use crate::pipeline::pipeline;
use crate::sequence::SequenceProcessor;
use clap::Parser;
use clio::Input;

/// Compute hash digests for sequences in a FASTA file
#[derive(Parser)]
#[command(version, about, max_term_width = 79)]
struct Cli {
    /// Input file(s). Use '-' for stdin
    #[clap(value_parser, default_value = "-")]
    input: Vec<Input>,

    /// Instead of hashing the entire sequence at once, hash each k-mer
    /// individually and then combine the resulting hashes
    #[clap(
        short = 'm',
        long,
        value_parser,
        default_value = "false",
        help_heading = "Hashing"
    )]
    multi_kmer_hashing: bool,

    /// Replace ntHash with xxHash for hashing k-mers. Works only with
    /// --multi-kmer-hashing
    #[clap(
        short = 'x',
        long,
        value_parser,
        requires = "multi_kmer_hashing",
        default_value = "false",
        help_heading = "Hashing"
    )]
    xxhash: bool,

    /// Size of the k-mers to hash when using --multi-kmer-hashing
    #[clap(
        short = 'k',
        long = "kmer-size",
        value_parser,
        requires = "multi_kmer_hashing",
        default_value = "31",
        help_heading = "Hashing"
    )]
    k: u8,

    /// Make hashing robust to circular permutations via deterministic rotation
    /// to the lexicographically minimal sequence
    #[clap(
        short = 'r',
        long,
        value_parser,
        default_value = "false",
        conflicts_with = "circular_kmers",
        help_heading = "Circular sequences"
    )]
    circular_rotation: bool,

    /// Make hashing robust to circular permutations via addition of the
    /// k-mers that wrap around the end of the sequence. Works only with
    /// --multi-kmer-hashing
    #[clap(
        short = 'w',
        long,
        value_parser,
        default_value = "false",
        requires = "multi_kmer_hashing",
        conflicts_with_all = ["circular_rotation"],
        help_heading = "Circular sequences"
    )]
    circular_kmers: bool,
}

fn main() {
    let cli = Cli::parse();

    let multi_kmer_hashing = cli.multi_kmer_hashing;
    let use_xxhash = cli.xxhash;
    let k: u8 = cli.k;
    let circular_rotation = cli.circular_rotation;
    let circular_kmers = cli.circular_kmers;

    let sequence_processor = SequenceProcessor::new(circular_rotation, circular_kmers, k);
    let hasher = SequenceHasher::new(multi_kmer_hashing, use_xxhash, k);

    for input in cli.input {
        pipeline(&input, &hasher, &sequence_processor);
    }
}
