use itertools::sorted;
use needletail::Sequence;
use nthash::NtHashIterator;
use std::error::Error;
use std::iter::Iterator;
use std::{fmt, hash::Hasher};
use twox_hash::xxh3::{hash64, Hash128, HasherExt};

enum KmerHasher<'a> {
    NtHash(NtHashIterator<'a>),
    XxHash(Box<dyn Iterator<Item = u64> + 'a>),
}

impl<'a> Iterator for KmerHasher<'a> {
    type Item = u64;
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            KmerHasher::NtHash(iter) => iter.next(),
            KmerHasher::XxHash(iter) => iter.next(),
        }
    }
}

#[derive(Debug)]
struct KmerHasherError {
    source: Box<dyn Error>,
}

impl fmt::Display for KmerHasherError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "HashIterator error: {}", self.source)
    }
}

impl From<String> for KmerHasherError {
    fn from(err: String) -> Self {
        KmerHasherError {
            source: Box::new(std::io::Error::new(std::io::ErrorKind::Other, err)),
        }
    }
}

struct KmerHashIterator<'a> {
    seq: &'a [u8],
    k: u8,
    use_xxhash: bool,
}

impl<'a> KmerHashIterator<'a> {
    fn new(seq: &'a [u8], k: u8, use_xxhash: bool) -> Self {
        KmerHashIterator { seq, k, use_xxhash }
    }

    fn get_kmer_hashes(&self) -> Result<impl Iterator<Item = u64> + 'a, KmerHasherError> {
        if self.use_xxhash {
            let rc = self.seq.reverse_complement();
            let iter = self
                .produce_xxhash_iterator(&rc)
                .map_err(|e| format!("Error: {}", e))?;
            Ok(sorted(iter))
        } else {
            let iter = self
                .produce_nthash_iterator()
                .map_err(|e| format!("Error: {}", e))?;
            Ok(sorted(iter))
        }
    }

    fn produce_nthash_iterator(&self) -> Result<KmerHasher<'a>, String> {
        match NtHashIterator::new(self.seq, self.k as usize) {
            Ok(iter) => Ok(KmerHasher::NtHash(iter)),
            Err(e) => Err(format!("Error: {}", e))?,
        }
    }

    fn produce_xxhash_iterator(&self, rc: &'a [u8]) -> Result<KmerHasher<'a>, String> {
        let iter = self
            .seq
            .canonical_kmers(self.k, rc)
            .map(|(_, kmer, _)| kmer)
            .map(hash64);
        Ok(KmerHasher::XxHash(Box::new(iter)))
    }
}

pub struct SequenceHasher {
    pub multi_kmer_hashing: bool,
    pub use_xxhash: bool,
    pub k: u8,
}

impl SequenceHasher {
    pub fn new(multi_kmer_hashing: bool, use_xxhash: bool, k: u8) -> Self {
        SequenceHasher {
            multi_kmer_hashing,
            use_xxhash,
            k,
        }
    }

    pub fn compute_hash(&self, seq: &[u8]) -> Result<u128, String> {
        match self.multi_kmer_hashing {
            true => self.compute_sequence_hash_multi_kmer(seq, self.k),
            false => self.compute_sequence_hash_single_kmer(seq),
        }
    }

    fn compute_sequence_hash_multi_kmer(&self, seq: &[u8], k: u8) -> Result<u128, String> {
        // Get the hashes of the k-mers in the sequence using either xxhash or nthash
        let kmer_hash_iterator = KmerHashIterator::new(seq, k, self.use_xxhash);
        let kmer_hashes = match kmer_hash_iterator.get_kmer_hashes() {
            Ok(hashes) => hashes,
            Err(e) => Err(format!("Error: {}", e))?,
        };
        // Combine the hashes of the k-mers to produce a single hash for the sequence
        let hash_seq = kmer_hashes
            .fold(Hash128::default(), |mut acc, h| {
                acc.write_u64(h);
                acc
            })
            .finish_ext();
        Ok(hash_seq)
    }

    fn compute_sequence_hash_single_kmer(&self, seq: &[u8]) -> Result<u128, String> {
        // Get the hash of the sequence using xxhash128
        let mut hasher = Hash128::default();
        hasher.write(seq);
        let hash_seq = hasher.finish_ext();
        Ok(hash_seq)
    }
}
