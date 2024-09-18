# seq-hasher

Compute hash digests for DNA sequences in a FASTA file, with support for circular permutations.

## Usage

```
seq-hasher [OPTIONS] [INPUT]...
```

### Arguments

| Argument | Description |
|:---------|:------------|
| `[INPUT]...` | Input file(s). Use `-` for stdin [default: `-`] |

### Options

| Option | Description |
|:-------|:------------|
| `-h`, `--help` | Print help |
| `-V`, `--version` | Print version |

### Hashing

| Option | Description |
|:-------|:------------|
| `-m`, `--multi-kmer-hashing` | Instead of hashing the entire sequence at once, hash each k-mer individually and then combine the resulting hashes |
| `-x`, `--xxhash` | Replace ntHash with xxHash for hashing k-mers. Works only with `--multi-kmer-hashing` |
| `-k`, `--kmer-size <K>` | Size of the k-mers to hash when using `--multi-kmer-hashing` [default: 31] |

### Circular sequences

| Option | Description |
|:-------|:------------|
| `r`, `--circular-rotation` | Make hashing robust to circular permutations via deterministic rotation to the lexicographically minimal sequence |
| `-w`, `--circular-kmers` | Make hashing robust to circular permutations via addition of the k-mers that wrap around the end of the sequence. Works only with `--multi-kmer-hashing` |

## Explanation

`seq-hasher` is a tool that processes sequences from FASTA files and generates a hash digest for each sequence. The default hashing method employs the `XXH3_128bits` algorithm from [xxHash](https://github.com/Cyan4973/xxHash) to hash the entire sequence at once.

By default, the entire sequence is hashed at once. However, you can switch to an alternative (and slower) approach using the `--multi-k-mer-hashing/-m` option. In this mode, the hash digest is created by hashing all k-mers of the sequence, sorting them, and then hashing the resulting hashes. By default, the k-mer size is set to 31 and the default hash function used for k-mers is [ntHash](https://github.com/bcgsc/ntHash). You can also switch the k-mer hashing algorithm to xxHash's `XXH3_64bit` using the `--xxhash/-x` option.

To ensure that circularly permuted sequences (e.g., "CATTTAA" and "TTAACAT") produce the same hash digest, `seq-hasher` provides the `--circular-rotation/-r` option, which deterministically rotates the sequences to their [lexicographically minimal sequence](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation) before hashing. Alternatively, if you use `--multi-k-mer-hashing`, you can make hashing robust to circular permutation through the `--circular-kmers/-w` option, which will include the extra k-mers that are formed by wrapping around the sequence in the hash computation.

## Examples

`seq-hasher` generates identical hash digests for sequences that differ only in character case or are reverse complements of each other.

```
$ echo -e ">seq_a\nCGAAACGTTCTT\n>seq_b\ncgaaacgttctt\n>seq_c\nAAGAACGTTTCG" | \
  seq-hasher
```

```
seq_a	8142a4bc45a6db851be0157ffafaaa68
seq_b	8142a4bc45a6db851be0157ffafaaa68
seq_c	8142a4bc45a6db851be0157ffafaaa68
```

By default, circularly permuted sequences produce different hash digests.

```
$ echo -e ">seq_a\nCGAAACGTTCTT\n>seq_d\nGTTCTTCGAAAC" | \
  seq-hasher
```

```
seq_a	8142a4bc45a6db851be0157ffafaaa68
seq_d	b3c4b8bcd66e6b51e237b2c40905dfb6
```

To make the hashing robust to circular permutations, use the `--rotate-circular/-r` option.

```
$ echo -e ">seq_a\nCGAAACGTTCTT\n>seq_d\nGTTCTTCGAAAC" | \
  seq-hasher -r
```

```
seq_a	790418b4049777375089c2537f77074a
seq_d	790418b4049777375089c2537f77074a
```

Compressed FASTA files in the `.gz`, `.bz2`, and `.xz` formats are supported.

```
$ curl -sfSLJ https://www.ebi.ac.uk/ena/browser/api/fasta/CP108542,AP024942 | \
  gzip > sequences.fna.gz
$ seq-hasher sequences.fna.gz
```

```
ENA|CP108542|CP108542.1	60e6f9f3d1f15fd44b63964f18f2fd97
ENA|AP024942|AP024942.1	18905b3c4ff019d2e8cc70839453edd1
```
