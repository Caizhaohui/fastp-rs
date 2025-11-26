use crate::fastq::FastqRecord;

pub struct OverlapResult {
    pub overlapped: bool,
    pub offset: i32,
    pub overlap_len: usize,
    pub diff: usize,
}

pub struct OverlapAnalyzer;

impl OverlapAnalyzer {
    pub fn analyze(r1: &FastqRecord, r2: &FastqRecord) -> OverlapResult {
        Self::analyze_with_params(r1, r2, 10, 5, 0.2)
    }

    pub fn analyze_with_params(r1: &FastqRecord, r2: &FastqRecord, min_overlap: usize, diff_limit: usize, diff_percent_limit: f32) -> OverlapResult {
        let len1 = r1.seq.len();
        let len2 = r2.seq.len();
        let seq1 = r1.seq.as_bytes();
        let r2_rc = reverse_complement(r2.seq.as_bytes());
        let s2 = r2_rc.as_bytes();

        let mut best_offset = 0;
        let mut best_diff = usize::MAX;
        let mut best_overlap_len = 0;
        let mut found = false;

        // Direction 1: R1 starts before S2 (offset >= 0)
        // R1: [offset] matches S2: [0]
        for offset in 0..len1 {
             let overlap_len = std::cmp::min(len1 - offset, len2);
             if overlap_len < min_overlap { continue; }
             
             let diff = count_diff(&seq1[offset..], &s2[0..overlap_len], overlap_len);
             let limit = std::cmp::min(diff_limit, (overlap_len as f32 * diff_percent_limit) as usize);
             
             if diff <= limit {
                 // Prefer larger overlap length? Or smaller diff ratio?
                 // Here we simply check if it's better than current best
                 if diff < best_diff || (diff == best_diff && overlap_len > best_overlap_len) {
                     best_diff = diff;
                     best_offset = offset as i32;
                     best_overlap_len = overlap_len;
                     found = true;
                 }
             }
        }

        // Direction 2: S2 starts before R1 (offset < 0)
        // S2: [offset_pos] matches R1: [0]
        // offset = -offset_pos
        for offset_pos in 1..len2 {
            let overlap_len = std::cmp::min(len2 - offset_pos, len1);
            if overlap_len < min_overlap { continue; }

            let diff = count_diff(&s2[offset_pos..], &seq1[0..overlap_len], overlap_len);
            let limit = std::cmp::min(diff_limit, (overlap_len as f32 * diff_percent_limit) as usize);

            if diff <= limit {
                 if diff < best_diff || (diff == best_diff && overlap_len > best_overlap_len) {
                     best_diff = diff;
                     best_offset = -(offset_pos as i32);
                     best_overlap_len = overlap_len;
                     found = true;
                 }
            }
        }

        OverlapResult {
            overlapped: found,
            offset: best_offset,
            overlap_len: best_overlap_len,
            diff: best_diff,
        }
    }
}

fn count_diff(s1: &[u8], s2: &[u8], len: usize) -> usize {
    let mut diff = 0;
    for i in 0..len {
        if s1[i] != s2[i] {
            diff += 1;
        }
    }
    diff
}


fn reverse_complement(seq: &[u8]) -> String {
    let mut res = String::with_capacity(seq.len());
    for b in seq.iter().rev() {
        res.push(match b {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'N' => 'N',
            _ => 'N',
        });
    }
    res
}
