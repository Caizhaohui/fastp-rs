use std::cmp;
use crate::fastq::FastqRecord;
use crate::filter::Report;
use crate::filter::matcher::Matcher;

pub struct AdapterTrimmer;

impl AdapterTrimmer {
    pub fn trim_by_sequence(
        rec: &mut FastqRecord,
        adapter_seq: Option<&str>,
        report: &mut Report,
        _is_r2: bool,
    ) -> bool {
        let adapter = match adapter_seq {
            Some(s) if !s.is_empty() => s,
            _ => return false,
        };

        // C++ defaults: matchReq depends on adapter list size, but for single sequence:
        // default matchReq = 4
        // allowOneMismatchForEach = 8
        let match_req = 4;
        let allow_one_mismatch_for_each = 8;

        let rlen = rec.seq.len();
        let alen = adapter.len();
        
        if alen < match_req {
            return false;
        }

        let rdata = rec.seq.as_bytes();
        let adata = adapter.as_bytes();

        let mut start = 0isize;
        if alen >= 16 {
            start = -4;
        } else if alen >= 12 {
            start = -3;
        } else if alen >= 8 {
            start = -2;
        }

        let mut pos = 0isize;
        let mut found = false;

        // 1. Exact match with hamming distance
        let end_pos = (rlen as isize) - (match_req as isize);
        
        for p in start..end_pos {
            pos = p;
            let cmplen = cmp::min((rlen as isize - pos) as usize, alen);
            let allowed_mismatch = cmplen / allow_one_mismatch_for_each;
            let mut mismatch = 0;
            let mut matched = true;

            let loop_start = cmp::max(0, -pos) as usize;
            for i in loop_start..cmplen {
                if adata[i] != rdata[(i as isize + pos) as usize] {
                    mismatch += 1;
                    if mismatch > allowed_mismatch {
                        matched = false;
                        break;
                    }
                }
            }
            if matched {
                found = true;
                break;
            }
        }

        // 2. Try one gap (insertion in sequence)
        if !found {
             for p in 0..((rlen as isize - match_req as isize - 1) as usize) {
                pos = p as isize;
                let cmplen = cmp::min(rlen - p - 1, alen);
                let allowed_mismatch = if cmplen / allow_one_mismatch_for_each > 0 {
                    (cmplen / allow_one_mismatch_for_each) - 1
                } else { 0 };
                
                // matchWithOneInsertion(rdata, adata, ...) -> insertion in rdata
                // C++: matchWithOneInsertion(rdata, adata, cmplen, allowedMismatch)
                // Rust Matcher expects (ins_data, normal_data)
                // Here rdata has insertion relative to adata? 
                // C++: Matcher::matchWithOneInsertion(rdata, adata, cmplen, allowedMismatch);
                // "rdata" is passed as insData. So read has insertion.
                let matched = Matcher::match_with_one_insertion(
                    &rdata[p..], 
                    adata, 
                    cmplen, 
                    allowed_mismatch
                );
                
                if matched {
                    found = true;
                    // hasInsertion = true;
                    break;
                }
            }
        }

        // 3. Try deletion in sequence (insertion in adapter)
        if !found {
            for p in 0..((rlen as isize - match_req as isize) as usize) {
                pos = p as isize;
                let cmplen = cmp::min(rlen - p, alen - 1);
                let allowed_mismatch = if cmplen / allow_one_mismatch_for_each > 0 {
                    (cmplen / allow_one_mismatch_for_each) - 1
                } else { 0 };

                // C++: Matcher::matchWithOneInsertion(adata, rdata, cmplen, allowedMismatch);
                // "adata" passed as insData. Adapter has insertion (meaning read has deletion).
                let matched = Matcher::match_with_one_insertion(
                    adata,
                    &rdata[p..],
                    cmplen,
                    allowed_mismatch
                );

                if matched {
                    found = true;
                    // hasDeletion = true;
                    break;
                }
            }
        }

        if found {
            if pos < 0 {
                // Adapter starts before read start (partial match at beginning?)
                // logic in C++:
                // string adapter = adapterseq.substr(0, alen+pos);
                // r->mSeq->resize(0);
                // r->mQuality->resize(0);
                // Actually if pos < 0, it means the adapter match starts BEFORE the read. 
                // Wait, if pos is negative, rlen-pos is > rlen. 
                // In C++, if pos < 0, it empties the read?
                // "Illumina adapter dimer usually have the first A skipped"
                // If the adapter is found at -4, it means the read STARTS with the adapter (minus first 4 bases).
                // So the whole read is adapter?
                // Yes, resize(0).
                
                // Capture trimmed bases for reporting
                // In C++: adapter = adapterseq.substr(0, alen+pos)
                // But we just track counts here for now
                let trimmed_len = rec.seq.len();
                rec.seq.clear();
                rec.qual.clear();
                
                report.adapter_trimmed_reads += 1;
                report.adapter_trimmed_bases += trimmed_len as u64;
            } else {
                let p = pos as usize;
                let trimmed_len = rec.seq.len() - p;
                rec.seq.truncate(p);
                rec.qual.truncate(p);
                
                report.adapter_trimmed_reads += 1;
                report.adapter_trimmed_bases += trimmed_len as u64;
            }
            return true;
        }

        false
    }
}
