use needletail::{parser::SequenceRecord, Sequence};

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
    pub k: u8,
}

impl SequenceProcessor {
    pub fn new(circular_rotation: bool, circular_kmers: bool, k: u8) -> Self {
        SequenceProcessor {
            circular_rotation,
            circular_kmers,
            k,
        }
    }

    pub fn process_sequence(&self, record: SequenceRecord) -> Vec<u8> {
        let norm_seq = record.normalize(false).to_vec();
        match (self.circular_rotation, self.circular_kmers) {
            (true, _) => self.lmsr_rotation(&norm_seq),
            (false, true) => self.adjust_for_circular_kmers(&norm_seq),
            _ => {
                let norm_seq_rc = norm_seq.reverse_complement();
                if norm_seq < norm_seq_rc {
                    norm_seq
                } else {
                    norm_seq_rc
                }
            }
        }
    }

    fn adjust_for_circular_kmers(&self, seq: &[u8]) -> Vec<u8> {
        let mut adjusted_seq = Vec::with_capacity(seq.len() + self.k as usize - 1);
        adjusted_seq.extend_from_slice(seq);
        adjusted_seq.extend_from_slice(&seq[..(self.k as usize - 1)]);
        adjusted_seq
    }

    fn compute_lmsr(&self, seq: &[u8]) -> Vec<u8> {
        let seq_len = seq.len();
        let doubled_seq: Vec<u8> = seq.iter().chain(seq.iter()).copied().collect();
        let (mut start_idx, mut min_rotation_idx) = (0, 0);
        while start_idx < seq_len {
            min_rotation_idx = start_idx;
            let (mut compare_idx, mut current_idx) = (start_idx + 1, start_idx);
            while compare_idx < 2 * seq_len && doubled_seq[current_idx] <= doubled_seq[compare_idx]
            {
                current_idx = if doubled_seq[current_idx] < doubled_seq[compare_idx] {
                    start_idx
                } else {
                    current_idx + 1
                };
                compare_idx += 1;
            }
            start_idx += compare_idx - current_idx;
        }
        doubled_seq[min_rotation_idx..min_rotation_idx + seq_len].to_vec()
    }

    fn lmsr_rotation(&self, seq: &[u8]) -> Vec<u8> {
        let rc = seq.reverse_complement();
        let lmsr_seq = self.compute_lmsr(seq);
        let lmsr_seq_rc = self.compute_lmsr(&rc);
        if lmsr_seq < lmsr_seq_rc {
            lmsr_seq
        } else {
            lmsr_seq_rc
        }
    }
}
