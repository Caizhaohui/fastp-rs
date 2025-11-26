use crate::fastq::FastqRecord;

pub struct BaseCorrector;

impl BaseCorrector {
    pub fn correct(r1: &mut FastqRecord, r2: &mut FastqRecord, offset: i32, overlap_len: usize) {
        let len2 = r2.seq.len();
        let mut s1 = r1.seq.as_bytes().to_vec();
        let mut q1 = r1.qual.as_bytes().to_vec();
        let mut s2 = r2.seq.as_bytes().to_vec();
        let mut q2 = r2.qual.as_bytes().to_vec();

        if offset >= 0 {
            let off = offset as usize;
            for i in 0..overlap_len {
                let i1 = off + i;
                let j = len2 - 1 - i;
                let b1 = s1[i1];
                let b2 = s2[j];
                let rc_b2 = complement(b2);
                if b1 != rc_b2 {
                    let q1v = q1[i1].saturating_sub(33);
                    let q2v = q2[j].saturating_sub(33);
                    if q1v >= q2v {
                        s2[j] = complement(b1);
                        q2[j] = q1[i1];
                    } else {
                        s1[i1] = rc_b2;
                        q1[i1] = q2[j];
                    }
                }
            }
        } else {
            let k = (-offset) as usize;
            for i in 0..overlap_len {
                let i1 = i;
                let j = len2 - 1 - (k + i);
                let b1 = s1[i1];
                let b2 = s2[j];
                let rc_b2 = complement(b2);
                if b1 != rc_b2 {
                    let q1v = q1[i1].saturating_sub(33);
                    let q2v = q2[j].saturating_sub(33);
                    if q1v >= q2v {
                        s2[j] = complement(b1);
                        q2[j] = q1[i1];
                    } else {
                        s1[i1] = rc_b2;
                        q1[i1] = q2[j];
                    }
                }
            }
        }

        r1.seq = String::from_utf8(s1).unwrap_or_default();
        r1.qual = String::from_utf8(q1).unwrap_or_default();
        r2.seq = String::from_utf8(s2).unwrap_or_default();
        r2.qual = String::from_utf8(q2).unwrap_or_default();
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

