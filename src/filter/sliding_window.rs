use crate::fastq::FastqRecord;
use crate::config::Cli;

pub struct SlidingWindow;

impl SlidingWindow {
    /// Sliding window trimming logic
    /// Ported from fastp C++ Filter::trimAndCut
    pub fn trim_and_cut(
        rec: &mut FastqRecord,
        config: &Cli,
        front_trimmed: &mut usize
    ) -> bool {
        let mut front = if *front_trimmed > 0 { *front_trimmed } else { 
            // This usually starts with static trim options which we handled before calling this?
            // In C++, trimAndCut takes 'front' and 'tail' arguments which are static trims.
            // But we already did static trimming in trim_record.
            // However, trim_record calls cut logic?
            // Let's assume we pass the already trimmed record, but we need to track cumulative trimming?
            0 
        };
        // We will operate on rec.seq directly, so "front" is relative to current rec.seq?
        // No, C++ logic:
        // int rlen = r->length() - front - tail;
        // The record passed to trimAndCut in C++ is the original record?
        // Let's see call site in C++:
        // mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed1);
        // So it combines static trim + quality cut.
        
        // In our Rust implementation, we already did static trimming in `trim_record`.
        // So `rec` passed here already has front/tail removed.
        // So we should treat front=0, tail=0 relative to `rec`.
        // But wait, if we cut more from front, we need to update front_trimmed for reporting/adapters?
        // Yes.
        
        let mut tail = 0;
        let l = rec.seq.len();
        
        // 1. Quality Cut Front (5')
        if config.cut_front {
            let w = config.cut_front_window_size;
            let mean_qual = config.cut_front_mean_quality as f32 + 33.0; // Phred+33
            
            if l >= w {
                let qual_bytes = rec.qual.as_bytes();
                let seq_bytes = rec.seq.as_bytes();
                
                let mut s = 0;
                let mut total_qual = 0;
                
                // Initialize window
                for i in 0..w-1 {
                    total_qual += qual_bytes[i] as i32;
                }
                
                let mut found = false;
                // Sliding
                for i in 0..(l - w + 1) {
                    s = i;
                    // Add new char at right of window
                    total_qual += qual_bytes[s + w - 1] as i32;
                    // Remove old char at left of window (if moved)
                    if s > 0 {
                        total_qual -= qual_bytes[s - 1] as i32;
                    }
                    
                    if (total_qual as f32 / w as f32) >= mean_qual {
                        found = true;
                        break;
                    }
                }
                
                if found {
                    // s is the start of the "good" window.
                    // C++: if(s > 0) s = s + w - 1;
                    // Why s + w - 1? 
                    // Example: w=4. Bad, Bad, Bad, Good(start at 3).
                    // Window [3,4,5,6] is good.
                    // It says we should cut up to... where?
                    // C++ comment: "the trimming in front is forwarded"
                    // If we find a good window at s, it means [s, s+w) is good.
                    // The code says: s = s + w - 1.
                    // This moves s to the LAST base of the good window.
                    // So we drop everything before that? That seems aggressive.
                    // Wait, let's re-read C++.
                    /*
                    if(s > 0)
                        s = s+w-1;
                    while(s<l && seq[s] == 'N')
                        s++;
                    front = s;
                    */
                    // If window [s, s+w) is the FIRST good window.
                    // It means [s-1, s+w-1) was bad (or s=0).
                    // If s=0, window [0, w) is good. s remains 0. front = 0. Keep all.
                    // If s=1, window [1, 1+w) is good. [0, w) was bad.
                    // s becomes 1 + w - 1 = w.
                    // So we cut [0..w].
                    // This implies if we slide until we find a good window, we cut everything BEFORE the END of that window?
                    // That sounds like we assume the "bad" region extends into the window?
                    // "fastp uses a sliding window ... if the mean quality is lower than required, the bases in the window are dropped."
                    // If we scan from front, we drop windows until we find a good one.
                    // But here we are looking for the FIRST good window.
                    // So all previous windows were bad.
                    // If [0,4) is bad, [1,5) is bad... [10,14) is good.
                    // Then we cut everything up to... 10? or 14?
                    // C++ says s = s + w - 1. So 10+4-1 = 13.
                    // So we cut 0..13. Start at 13.
                    // This means we cut the window ITSELF too? except the last base?
                    // If s=0 (first window good), s remains 0.
                    // If s > 0, s jumps to end of window.
                    
                    if s > 0 {
                        s = s + w - 1;
                    }
                    
                    // Also skip Ns
                    while s < l && seq_bytes[s] == b'N' {
                        s += 1;
                    }
                    
                    front = s;
                } else {
                    // No good window found, cut everything?
                    // C++: loop finishes, s will be l-w+1?
                    // Actually if break not hit, s goes up to l-tail-w.
                    // If not found, we probably drop everything?
                    // C++ doesn't explicitly handle "not found" inside the if(mOptions->qualityCut.enabledFront).
                    // But if loop finishes without break, s will be large.
                    // But wait, the loop variable s is local to the for loop scope in C++?
                    // "int s = front;" declared before loop.
                    // Loop: "for(s=front; s+w<l-tail; s++)"
                    // If loop completes, s = l - tail - w.
                    // Then "if(s > 0) s = s+w-1" -> s = l - tail - 1.
                    // So we trim almost everything.
                    if !found {
                        front = l; // Trim all
                    }
                }
            } else {
                // Too short for window, trim all? or keep all?
                // C++: if(l - front - tail - w <= 0) return NULL; -> Discard read
                front = l;
            }
        }
        
        // 2. Quality Cut Right (Sliding from Front, but cut tail if bad window found)
        // Wait, "cut_right" in fastp means:
        // "Scan from beginning, if we find a bad window, we cut everything from there to the end."
        // This is different from cut_tail (which scans from end).
        let mut right_cut_len = 0; // Amount to cut from right
        if config.cut_right {
             let w = config.cut_right_window_size;
             let mean_qual = config.cut_right_mean_quality as f32 + 33.0;
             let current_l = l.saturating_sub(front); // Length after front trim
             
             if current_l >= w {
                 let qual_bytes = &rec.qual.as_bytes()[front..];
                 // let seq_bytes = &rec.seq.as_bytes()[front..]; // Not needed for right cut?
                 
                 let mut s;
                 let mut total_qual = 0;
                 
                 // Init window
                 for i in 0..w-1 {
                     total_qual += qual_bytes[i] as i32;
                 }
                 
                 let mut found_bad = false;
                 let mut bad_pos = 0;
                 
                 for i in 0..(current_l - w + 1) {
                     s = i;
                     total_qual += qual_bytes[s + w - 1] as i32;
                     if s > 0 {
                         total_qual -= qual_bytes[s - 1] as i32;
                     }
                     
                     if (total_qual as f32 / w as f32) < mean_qual {
                         found_bad = true;
                         bad_pos = s;
                         break;
                     }
                 }
                 
                 if found_bad {
                     // Found a bad window starting at bad_pos.
                     // We want to keep bases BEFORE this window.
                     // But inside the bad window, maybe some bases are good?
                     // C++:
                     /*
                     if(foundLowQualWindow ) {
                        // keep the good bases in the window
                        while(s<l-1 && qualstr[s]>=33 + mOptions->qualityCut.qualityRight)
                            s++;
                        rlen = s - front;
                    }
                     */
                     // s is relative to qualstr (which includes front trim in C++? No, C++ uses index s starting from front).
                     // Here we used slice starting at front. So s is relative to slice.
                     // We check bases starting at bad_pos.
                     
                     let mut check_pos = bad_pos;
                     while check_pos < current_l && qual_bytes[check_pos] >= mean_qual as u8 {
                         check_pos += 1;
                     }
                     // Cut from check_pos to end.
                     // New length is check_pos.
                     // So we cut (current_l - check_pos) from right.
                     right_cut_len = current_l - check_pos;
                 }
             } else {
                 // Too short, discard?
                 right_cut_len = current_l;
             }
        }
        
        // 3. Quality Cut Tail (3') - Scan from tail to front
        // Only if cut_right is NOT enabled (C++ logic prefers cut_right over cut_tail if both? No, "if(!mOptions->qualityCut.enabledRight && mOptions->qualityCut.enabledTail)")
        // So they are mutually exclusive in C++ implementation for the tail side?
        // Yes.
        if !config.cut_right && config.cut_tail {
             let w = config.cut_tail_window_size;
             let mean_qual = config.cut_tail_mean_quality as f32 + 33.0;
             let current_l = l.saturating_sub(front);
             
             if current_l >= w {
                 let qual_bytes = &rec.qual.as_bytes()[front..];
                 // We scan from end.
                 // C++: t = l - tail - 1;
                 // loop t down to front.
                 
                 let mut total_qual = 0;
                 let mut t = current_l - 1;
                 
                 // Init window at very end
                 // Window: [t-w+1, t]
                 // Init: [current_l-w, current_l-1]
                 // Note: C++ loop seems to slide window backwards.
                 
                 for i in 0..w-1 {
                     // totalQual += qualstr[t-i]; -> sum of last w-1 bases
                     total_qual += qual_bytes[t - i] as i32;
                 }
                 
                 let mut found_good = false;
                 // Loop t from end down to w-1
                 // Window ends at t, starts at t-w+1
                 // We want to find the first GOOD window from the back (which means closest to end that is good? Or closest to front?)
                 // C++: "for(t=l-tail-1; t-w>=front; t--)"
                 // It scans from right to left.
                 // Finds FIRST window (from right) that is GOOD.
                 // Break.
                 
                 // Re-implementation:
                 // t is the END index of the window (inclusive)
                 let mut best_t = 0; 
                 
                 // We need to iterate such that we cover all windows.
                 // Start t at current_l - 1.
                 // End when t - w + 1 == 0.
                 
                 for i in (w-1..current_l).rev() {
                     t = i;
                     // Add left-most char of new window (which is t-w+1)
                     total_qual += qual_bytes[t - w + 1] as i32;
                     
                     // Remove char that fell off the right (t+1)
                     if t < current_l - 1 {
                         total_qual -= qual_bytes[t + 1] as i32;
                     }
                     
                     if (total_qual as f32 / w as f32) >= mean_qual {
                         found_good = true;
                         best_t = t;
                         break;
                     }
                 }
                 
                 if found_good {
                     // Found good window ending at best_t.
                     // C++: if(t < l-1) t = t-w+1;
                     // "t = t - w + 1" moves position to START of the good window.
                     // Then "while(t>=0 && seq[t] == 'N') t--;"
                     // Then "rlen = t - front + 1;" (so t is the last kept index)
                     // So we keep up to start of good window?
                     // If [10, 14) is good. t=13. t becomes 10.
                     // We keep 0..10. (Length 11).
                     // Wait, if [10, 14) is good, why cut 11,12,13?
                     // "quality cutting backward" -> we drop bases from tail until we find a good window.
                     // So if we find a good window, we keep it?
                     // C++ logic: "t = t - w + 1". This sets t to the START of the window.
                     // Then rlen = t ...
                     // So we KEEP the good window? No, t is the LAST index.
                     // So if t=10, we keep index 10.
                     // So we keep the first base of the good window, and drop the rest of the window?
                     // That seems weird.
                     // Let's trace carefully.
                     // C++: "if(t < l-1) t = t-w+1;"
                     // If we found the good window at the very end (t = l-1), then t doesn't change. We keep everything.
                     // If we found it earlier (t < l-1), we move t to start of window.
                     // This implies that if we had to trim something (found good window "inside"), 
                     // we assume the good window marks the END of the valid read?
                     // Wait, we are scanning from TAIL.
                     // Bad, Bad, Bad, Good.
                     // We discard the Bads. We keep the Good.
                     // So we should cut everything AFTER the Good window.
                     // The Good window [start, end].
                     // If we cut after Good, we keep up to 'end'.
                     // But C++ sets t = start.
                     // This means we keep up to 'start'.
                     // So we discard the Good window too (except first base)?
                     // Maybe because sliding window average is high, but individual bases towards end might be low?
                     // Or maybe it's just heuristic.
                     // Let's stick to C++ logic: t = t - w + 1.
                     
                     if best_t < current_l - 1 {
                         best_t = best_t + 1 - w;
                     }
                     
                     // Check Ns backwards
                     let seq_bytes = &rec.seq.as_bytes()[front..];
                     while best_t > 0 && seq_bytes[best_t] == b'N' {
                         best_t -= 1;
                     }
                     
                     // Calculate cut length
                     // Kept length is best_t + 1.
                     tail = current_l - (best_t + 1);
                 } else {
                     // No good window found, trim all
                     tail = current_l;
                 }
             } else {
                 tail = current_l;
             }
        }
        
        // Apply Trimming
        if front > 0 || tail > 0 || right_cut_len > 0 {
            let total_cut_right = tail.max(right_cut_len);
            let new_len = l.saturating_sub(front).saturating_sub(total_cut_right);
            
            if new_len == 0 {
                rec.seq.clear();
                rec.qual.clear();
                return false; // filtered out (effectively empty)
            }
            
            // Do the substring
            // rec.seq = rec.seq[front .. l - total_cut_right]
            let end = l - total_cut_right;
            if front >= end {
                 rec.seq.clear();
                 rec.qual.clear();
                 return false;
            }
            
            rec.seq = rec.seq[front..end].to_string();
            rec.qual = rec.qual[front..end].to_string();
            
            *front_trimmed += front;
            return true;
        }
        
        true
    }
}
