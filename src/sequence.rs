use needletail::{parser::SequenceRecord, Sequence};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::num::NonZeroU8;

pub fn get_record_accession(record_header: &[u8]) -> Option<&[u8]> {
    let accession = record_header
        .split(|&b| matches!(b, b' ' | b'\t' | b'\n' | b'\x0C' | b'\r'))
        .next();
    match accession {
        Some(acc) if !acc.is_empty() => Some(acc),
        _ => None,
    }
}

pub struct SequenceProcessor {
    pub circular_rotation: bool,
    pub circular_kmers: bool,
    pub k: NonZeroU8,
}

impl SequenceProcessor {
    pub fn new(circular_rotation: bool, circular_kmers: bool, k: NonZeroU8) -> Self {
        SequenceProcessor {
            circular_rotation,
            circular_kmers,
            k,
        }
    }

    pub fn process_sequence<'a>(&self, record: &'a SequenceRecord<'a>) -> Cow<'a, [u8]> {
        let norm_seq = record.normalize(false);
        match (self.circular_rotation, self.circular_kmers) {
            (true, _) => Cow::Owned(self.lmsr_rotation(norm_seq.as_ref())),
            (false, true) => Cow::Owned(self.adjust_for_circular_kmers(norm_seq.as_ref())),
            _ => {
                let norm_seq_rc = norm_seq.as_ref().reverse_complement();
                if norm_seq.as_ref() < norm_seq_rc.as_slice() {
                    norm_seq
                } else {
                    Cow::Owned(norm_seq_rc)
                }
            }
        }
    }

    fn adjust_for_circular_kmers(&self, seq: &[u8]) -> Vec<u8> {
        let wrap_len = usize::from(self.k.get()) - 1;
        let mut adjusted_seq = Vec::with_capacity(seq.len() + wrap_len);
        adjusted_seq.extend_from_slice(seq);
        adjusted_seq.extend_from_slice(&seq[..wrap_len]);
        adjusted_seq
    }

    #[inline]
    fn circular_index(index: usize, len: usize) -> usize {
        debug_assert!(len > 0);
        if index < len {
            index
        } else {
            debug_assert!(index - len < len);
            index - len
        }
    }

    fn minimal_rotation_index(seq: &[u8]) -> usize {
        let seq_len = seq.len();
        if seq_len <= 1 {
            return 0;
        }

        let (mut left, mut right, mut offset) = (0, 1, 0);
        while left < seq_len && right < seq_len && offset < seq_len {
            let left_base = seq[Self::circular_index(left + offset, seq_len)];
            let right_base = seq[Self::circular_index(right + offset, seq_len)];
            match left_base.cmp(&right_base) {
                Ordering::Equal => offset += 1,
                Ordering::Greater => {
                    left += offset + 1;
                    if left <= right {
                        left = right + 1;
                    }
                    offset = 0;
                }
                Ordering::Less => {
                    right += offset + 1;
                    if right <= left {
                        right = left + 1;
                    }
                    offset = 0;
                }
            }
        }

        usize::min(left, right)
    }

    fn compare_rotations(
        seq: &[u8],
        seq_start: usize,
        other: &[u8],
        other_start: usize,
    ) -> Ordering {
        let len = seq.len();
        debug_assert_eq!(len, other.len());
        for offset in 0..len {
            let seq_base = seq[Self::circular_index(seq_start + offset, len)];
            let other_base = other[Self::circular_index(other_start + offset, len)];
            let ordering = seq_base.cmp(&other_base);
            if ordering != Ordering::Equal {
                return ordering;
            }
        }
        Ordering::Equal
    }

    fn build_rotation(seq: &[u8], start: usize) -> Vec<u8> {
        let mut rotation = Vec::with_capacity(seq.len());
        rotation.extend_from_slice(&seq[start..]);
        rotation.extend_from_slice(&seq[..start]);
        rotation
    }

    fn lmsr_rotation(&self, seq: &[u8]) -> Vec<u8> {
        let rc = seq.reverse_complement();
        let seq_start = Self::minimal_rotation_index(seq);
        let rc_start = Self::minimal_rotation_index(&rc);
        if Self::compare_rotations(seq, seq_start, &rc, rc_start) == Ordering::Less {
            Self::build_rotation(seq, seq_start)
        } else {
            Self::build_rotation(&rc, rc_start)
        }
    }
}
