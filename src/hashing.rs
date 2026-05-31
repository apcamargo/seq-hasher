use needletail::Sequence;
use nthash::NtHashIterator;
use std::num::NonZeroU8;
use twox_hash::{XxHash3_128, XxHash3_64};

pub struct SequenceHasher {
    pub multi_kmer_hashing: bool,
    pub use_xxhash: bool,
    pub k: NonZeroU8,
}

impl SequenceHasher {
    pub fn new(multi_kmer_hashing: bool, use_xxhash: bool, k: NonZeroU8) -> Self {
        SequenceHasher {
            multi_kmer_hashing,
            use_xxhash,
            k,
        }
    }

    pub fn compute_hash(&self, seq: &[u8]) -> Result<u128, String> {
        if self.multi_kmer_hashing {
            self.compute_sequence_hash_multi_kmer(seq)
        } else {
            Ok(Self::compute_sequence_hash_single_kmer(seq))
        }
    }

    fn compute_sequence_hash_multi_kmer(&self, seq: &[u8]) -> Result<u128, String> {
        let kmer_hashes = if self.use_xxhash {
            Self::collect_xxhash_kmer_hashes(seq, self.k)
        } else {
            Self::collect_nthash_kmer_hashes(seq, self.k)?
        };
        Ok(Self::combine_kmer_hashes(kmer_hashes))
    }

    fn collect_nthash_kmer_hashes(seq: &[u8], k: NonZeroU8) -> Result<Vec<u64>, String> {
        let mut hashes: Vec<u64> = NtHashIterator::new(seq, usize::from(k.get()))
            .map_err(|e| format!("Error: {e}"))?
            .collect();
        hashes.sort_unstable();
        Ok(hashes)
    }

    fn collect_xxhash_kmer_hashes(seq: &[u8], k: NonZeroU8) -> Vec<u64> {
        let rc = seq.reverse_complement();
        let mut hashes: Vec<u64> = seq
            .canonical_kmers(k.get(), &rc)
            .map(|(_, kmer, _)| XxHash3_64::oneshot(kmer))
            .collect();
        hashes.sort_unstable();
        hashes
    }

    fn combine_kmer_hashes(kmer_hashes: Vec<u64>) -> u128 {
        kmer_hashes
            .into_iter()
            .fold(XxHash3_128::default(), |mut acc, hash| {
                acc.write(&hash.to_ne_bytes());
                acc
            })
            .finish_128()
    }

    fn compute_sequence_hash_single_kmer(seq: &[u8]) -> u128 {
        XxHash3_128::oneshot(seq)
    }
}
