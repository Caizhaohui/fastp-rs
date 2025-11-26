// use std::cmp;

pub struct Matcher;

impl Matcher {
    /// Checks if a sequence matches another with exactly one insertion in the first sequence (ins_data).
    /// 
    /// # Arguments
    /// * `ins_data` - The sequence with potential insertion (length should be at least cmplen + 1)
    /// * `normal_data` - The reference sequence (length should be at least cmplen)
    /// * `cmplen` - The length of the comparison window in normal_data
    /// * `diff_limit` - The maximum allowed mismatches (excluding the insertion itself)
    /// 
    /// # Logic
    /// We try all possible positions for the insertion.
    /// If an insertion is at index `i` (1-based relative to comparison window),
    /// then we compare:
    /// - `ins_data[0..i]` with `normal_data[0..i]` (left part)
    /// - `ins_data[i+1..cmplen+1]` with `normal_data[i..cmplen]` (right part)
    pub fn match_with_one_insertion(ins_data: &[u8], normal_data: &[u8], cmplen: usize, diff_limit: usize) -> bool {
        if ins_data.len() < cmplen + 1 || normal_data.len() < cmplen {
            return false;
        }

        let mut acc_mismatch_from_left = vec![0; cmplen];
        let mut acc_mismatch_from_right = vec![0; cmplen];

        // acc_mismatch_from_left[0]: head vs. head
        acc_mismatch_from_left[0] = if ins_data[0] == normal_data[0] { 0 } else { 1 };
        
        // acc_mismatch_from_right[cmplen-1]: tail vs. tail
        acc_mismatch_from_right[cmplen - 1] = if ins_data[cmplen] == normal_data[cmplen - 1] { 0 } else { 1 };

        // Calculate accumulated mismatches from left
        for i in 1..cmplen {
            let mismatch = if ins_data[i] == normal_data[i] { 0 } else { 1 };
            acc_mismatch_from_left[i] = acc_mismatch_from_left[i - 1] + mismatch;

            // Early exit optimization? (from C++ code: if accumulated > diffLimit, break)
            // But we need to compute full array for later check or can we?
            // The C++ code breaks early, which means subsequent values are not computed/valid.
            // However, the C++ code uses separate loops for left and right.
        }

        // Calculate accumulated mismatches from right
        for i in (0..cmplen - 1).rev() {
            let mismatch = if ins_data[i + 1] == normal_data[i] { 0 } else { 1 };
            acc_mismatch_from_right[i] = acc_mismatch_from_right[i + 1] + mismatch;
        }

        // Check for possible insertion points
        // Insertion can be from pos = 1 to cmplen - 1
        for i in 1..cmplen {
            // Check if left part alone exceeds limit (optimization from C++)
            // In C++: if(accMismatchFromLeft[i-1] + accMismatchFromRight[cmplen-1]> diffLimit) return false;
            // This seems to check if even with "perfect" right match (which is not what accMismatchFromRight[cmplen-1] is?), 
            // wait, accMismatchFromRight[cmplen-1] is just the last char comparison.
            // Let's stick to the core logic:
            // diff = mismatches before insertion + mismatches after insertion
            
            let diff = acc_mismatch_from_left[i - 1] + acc_mismatch_from_right[i];
            if diff <= diff_limit {
                return true;
            }
        }

        false
    }
}
