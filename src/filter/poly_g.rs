use crate::fastq::FastqRecord;

pub struct PolyGTrimmer;

impl PolyGTrimmer {
    pub fn trim_poly_g(rec: &mut FastqRecord, min_len: usize, report: &mut super::Report) {
        const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
        const MAX_MISMATCH: usize = 5;

        let seq = rec.seq.as_bytes();
        let rlen = seq.len();

        let mut mismatch = 0;
        let mut i = 0;
        let mut first_g_pos = rlen; // If no poly G, no trim

        // Scan from tail to head
        // i is number of bases checked from tail
        for check_idx in 0..rlen {
            i = check_idx;
            // index from tail: rlen - 1 - check_idx
            if seq[rlen - 1 - check_idx] != b'G' {
                mismatch += 1;
            } else {
                first_g_pos = rlen - 1 - check_idx;
            }

            let allowed_mismatch = (check_idx + 1) / ALLOW_ONE_MISMATCH_FOR_EACH;
            if mismatch > MAX_MISMATCH || (mismatch > allowed_mismatch && check_idx >= min_len - 1) {
                break;
            }
        }

        // if the polyG length (i) >= min_len, trim it
        if i >= min_len {
            if first_g_pos < rlen {
                 let trimmed_len = rlen - first_g_pos;
                 // Update report (assuming we track polyG trimmed reads/bases?)
                 // Currently report struct doesn't have polyG fields, but we can reuse or add them.
                 // For now, let's just trim.
                 // To match C++, we need to add fields to Report?
                 // User asked to implement the function.
                 
                 // Note: C++ uses `firstGPos` which is the index of the first G in the valid PolyG tail.
                 // So we truncate to `firstGPos`.
                 report.poly_g_trimmed_reads += 1;
                 report.poly_g_trimmed_bases += trimmed_len as u64;
                 rec.seq.truncate(first_g_pos);
                 rec.qual.truncate(first_g_pos);
            }
        }
    }
}
