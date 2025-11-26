mod matcher;
mod adapter_trimmer;
mod sliding_window;
mod overlap;
mod poly_g;
mod poly_x;
mod base_correction;

use serde::Serialize;
use crate::fastq::FastqRecord;
use crate::config::Cli;
use self::adapter_trimmer::AdapterTrimmer;
use self::sliding_window::SlidingWindow;
use self::overlap::OverlapAnalyzer;
use self::poly_g::PolyGTrimmer;
use self::poly_x::PolyXTrimmer;
use self::base_correction::BaseCorrector;

#[derive(Default, Serialize, Clone)]
pub struct Report {
    pub total_reads: u64,
    pub passed_reads: u64,
    pub failed_too_short: u64,
    pub failed_low_quality: u64,
    pub failed_n_excess: u64,
    pub failed_low_average_qual: u64,
    pub adapter_trimmed_reads: u64,
    pub adapter_trimmed_bases: u64,
    pub poly_g_trimmed_reads: u64,
    pub poly_g_trimmed_bases: u64,
    pub poly_x_trimmed_reads: u64,
    pub poly_x_trimmed_bases: u64,
    pub pe_overlap_avg_diff: f32,
    pub pe_overlap_count: u64,
}

impl Report {
    pub fn merge(&mut self, other: &Report) {
        self.total_reads += other.total_reads;
        self.passed_reads += other.passed_reads;
        self.failed_too_short += other.failed_too_short;
        self.failed_low_quality += other.failed_low_quality;
        self.failed_n_excess += other.failed_n_excess;
        self.failed_low_average_qual += other.failed_low_average_qual;
        self.adapter_trimmed_reads += other.adapter_trimmed_reads;
        self.adapter_trimmed_bases += other.adapter_trimmed_bases;
        self.poly_g_trimmed_reads += other.poly_g_trimmed_reads;
        self.poly_g_trimmed_bases += other.poly_g_trimmed_bases;
        self.poly_x_trimmed_reads += other.poly_x_trimmed_reads;
        self.poly_x_trimmed_bases += other.poly_x_trimmed_bases;
        // weighted average for overlap diff
        let total = self.pe_overlap_count + other.pe_overlap_count;
        if total > 0 {
            let sum = self.pe_overlap_avg_diff * (self.pe_overlap_count as f32)
                    + other.pe_overlap_avg_diff * (other.pe_overlap_count as f32);
            self.pe_overlap_avg_diff = sum / (total as f32);
            self.pe_overlap_count = total;
        }
    }
}

pub struct Filter {
    config: Cli,
}

impl Filter {
    pub fn new(config: Cli) -> Self {
        Self { config }
    }

    pub fn trim_pair(&self, mut r1: FastqRecord, mut r2: FastqRecord, report: &mut Report) -> (FastqRecord, FastqRecord) {
        let res = OverlapAnalyzer::analyze_with_params(&r1, &r2, self.config.overlap_len_require, self.config.overlap_diff_limit, (self.config.overlap_diff_percent_limit as f32) / 100.0);
        if res.overlapped {
            // accumulate average diff
            let prev = report.pe_overlap_avg_diff * report.pe_overlap_count as f32;
            report.pe_overlap_count += 1;
            report.pe_overlap_avg_diff = (prev + res.diff as f32) / report.pe_overlap_count as f32;
        }
        if self.config.correction && res.overlapped && res.overlap_len >= self.config.overlap_len_require {
            BaseCorrector::correct(&mut r1, &mut r2, res.offset, res.overlap_len);
        }
        if !self.config.disable_adapter_trimming 
           && self.config.adapter_sequence.is_none() 
           && self.config.adapter_sequence_r2.is_none() {
            if res.overlapped {
                let offset = res.offset;
                let overlap_len = res.overlap_len;
                if offset >= 0 {
                    let off = offset as usize;
                    if r1.seq.len() > off + overlap_len {
                        let trimmed_len = r1.seq.len() - (off + overlap_len);
                        report.adapter_trimmed_bases += trimmed_len as u64;
                        r1.seq.truncate(off + overlap_len);
                        r1.qual.truncate(off + overlap_len);
                        report.adapter_trimmed_reads += 1;
                    }
                } else {
                    let k = (-offset) as usize;
                    if r1.seq.len() > overlap_len {
                        let trimmed_len = r1.seq.len() - overlap_len;
                        report.adapter_trimmed_bases += trimmed_len as u64;
                        r1.seq.truncate(overlap_len);
                        r1.qual.truncate(overlap_len);
                        report.adapter_trimmed_reads += 1;
                    }
                    if r2.seq.len() > k {
                        let new_len = r2.seq.len().saturating_sub(k);
                        let trimmed_len = r2.seq.len() - new_len;
                        if trimmed_len > 0 {
                            report.adapter_trimmed_bases += trimmed_len as u64;
                            r2.seq.truncate(new_len);
                            r2.qual.truncate(new_len);
                            report.adapter_trimmed_reads += 1;
                        }
                    }
                }
            }
        }
        
        // Continue with individual record trimming (Quality, etc.)
        let r1 = self.trim_record(r1, false, report);
        let r2 = self.trim_record(r2, true, report);
        
        (r1, r2)
    }

    pub fn trim_record(&self, mut rec: FastqRecord, is_r2: bool, report: &mut Report) -> FastqRecord {
        let front = if is_r2 { self.config.trim_front2 } else { self.config.trim_front1 };
        let tail = if is_r2 { self.config.trim_tail2 } else { self.config.trim_tail1 };
        let max_len = if is_r2 { self.config.max_len2 } else { self.config.max_len1 };

        // Adapter Trimming
        if !self.config.disable_adapter_trimming {
            let adapter_seq = if is_r2 {
                self.config.adapter_sequence_r2.as_deref()
            } else {
                self.config.adapter_sequence.as_deref()
            };
            
            // If explicit adapter sequence is provided, use it
            if let Some(_) = adapter_seq {
                AdapterTrimmer::trim_by_sequence(&mut rec, adapter_seq, report, is_r2);
            }
        }
        
        if self.config.trim_poly_x {
            PolyXTrimmer::trim_poly_x(&mut rec, self.config.poly_x_min_len, report);
        }
        if self.config.trim_poly_g && !self.config.disable_trim_poly_g {
            PolyGTrimmer::trim_poly_g(&mut rec, self.config.poly_g_min_len, report);
        }

        // Sliding Window Quality Cutting
        let mut front_trimmed_count = 0;
        SlidingWindow::trim_and_cut(&mut rec, &self.config, &mut front_trimmed_count);

        let start = front.min(rec.seq.len());
        let mut end = rec.seq.len().saturating_sub(tail);
        
        if max_len > 0 {
            let max_end = start.saturating_add(max_len);
            end = end.min(max_end);
        }
        
        if start >= end { 
            rec.seq.clear(); 
            rec.qual.clear(); 
            return rec; 
        }
        
        let seq = rec.seq[start..end].to_string();
        let qual = rec.qual[start..end].to_string();
        rec.seq = seq;
        rec.qual = qual;
        rec
    }

    pub fn pass_filters(&self, rec: &FastqRecord, rep: &mut Report) -> bool {
        let qmin = self.config.qualified_quality_phred;
        let unq_limit = self.config.unqualified_percent_limit;
        let len_req = self.config.length_required;
        let n_limit = self.config.n_base_limit;
        let avg_req = self.config.average_qual;

        if rec.seq.len() < len_req { 
            rep.failed_too_short += 1; 
            return false; 
        }
        
        let n_count = rec.seq.bytes().filter(|&b| b == b'N' || b == b'n').count();
        if n_count > n_limit { 
            rep.failed_n_excess += 1; 
            return false; 
        }
        
        if avg_req > 0 && self.avg_phred(&rec.qual) < avg_req as f32 { 
            rep.failed_low_average_qual += 1; 
            return false; 
        }
        
        let mut low = 0usize;
        for b in rec.qual.bytes() { 
            if (b.saturating_sub(33)) < qmin { 
                low += 1; 
            } 
        }
        
        let pct = if rec.qual.is_empty() { 
            100.0 
        } else { 
            (low as f32) * 100.0 / rec.qual.len() as f32 
        };
        
        if pct > unq_limit as f32 { 
            rep.failed_low_quality += 1; 
            return false; 
        }
        
        true
    }

    fn avg_phred(&self, q: &str) -> f32 { 
        if q.is_empty() {
            0.0
        } else { 
            q.bytes().map(|b| (b.saturating_sub(33)) as f32).sum::<f32>() / q.len() as f32 
        } 
    }
}
