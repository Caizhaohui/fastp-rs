use crate::fastq::FastqRecord;

pub struct PolyXTrimmer;

impl PolyXTrimmer {
    pub fn trim_poly_x(rec: &mut FastqRecord, min_len: usize, report: &mut super::Report) {
        const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
        const MAX_MISMATCH: usize = 5;

        let seq = rec.seq.as_bytes();
        let rlen = seq.len();
        if rlen == 0 { return; }

        let mut mismatch = 0usize;
        let mut i = 0usize;
        let mut first_pos = rlen;
        let mut base: u8 = seq[rlen - 1];

        for check_idx in 0..rlen {
            i = check_idx;
            let cur = seq[rlen - 1 - check_idx];
            if check_idx == 0 { base = cur; }
            if cur != base { mismatch += 1; } else { first_pos = rlen - 1 - check_idx; }

            let allowed_mismatch = (check_idx + 1) / ALLOW_ONE_MISMATCH_FOR_EACH;
            if mismatch > MAX_MISMATCH || (mismatch > allowed_mismatch && check_idx >= min_len - 1) { break; }
        }

        if i >= min_len && first_pos < rlen {
            let trimmed_len = rlen - first_pos;
            report.poly_x_trimmed_reads += 1;
            report.poly_x_trimmed_bases += trimmed_len as u64;
            rec.seq.truncate(first_pos);
            rec.qual.truncate(first_pos);
        }
    }
}
