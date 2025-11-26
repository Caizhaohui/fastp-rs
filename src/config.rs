use clap::{Parser, ArgAction};
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(name = "fastp-rs", version = "0.1.0", about = "FASTQ preprocessor (Rust)")]
pub struct Cli {
    #[arg(short='i', long="in1")]
    pub in1: Option<String>,
    #[arg(short='I', long="in2")]
    pub in2: Option<String>,
    #[arg(short='o', long="out1")]
    pub out1: Option<String>,
    #[arg(short='O', long="out2")]
    pub out2: Option<String>,
    #[arg(long="stdin", action=ArgAction::SetTrue)]
    pub stdin: bool,
    #[arg(long="stdout", action=ArgAction::SetTrue)]
    pub stdout: bool,
    
    // Trimming Options
    #[arg(short='f', long="trim_front1", default_value_t=0)]
    pub trim_front1: usize,
    #[arg(short='t', long="trim_tail1", default_value_t=0)]
    pub trim_tail1: usize,
    #[arg(short='b', long="max_len1", default_value_t=0)]
    pub max_len1: usize,
    
    #[arg(short='F', long="trim_front2", default_value_t=0)]
    pub trim_front2: usize,
    #[arg(short='T', long="trim_tail2", default_value_t=0)]
    pub trim_tail2: usize,
    #[arg(short='B', long="max_len2", default_value_t=0)]
    pub max_len2: usize,
    
    // Filtering Options
    #[arg(short='l', long="length_required", default_value_t=15)]
    pub length_required: usize,
    #[arg(short='q', long="qualified_quality_phred", default_value_t=15)]
    pub qualified_quality_phred: u8,
    #[arg(short='u', long="unqualified_percent_limit", default_value_t=40)]
    pub unqualified_percent_limit: u8,
    #[arg(short='e', long="average_qual", default_value_t=0)]
    pub average_qual: u8,
    #[arg(short='n', long="n_base_limit", default_value_t=5)]
    pub n_base_limit: usize,
    
    // Reporting
    #[arg(short='j', long="json", default_value = "fastp.json")]
    pub json: String,
    #[arg(long="html", default_value = "fastp.html")]
    pub html: String,
    #[arg(short='R', long="report_title", default_value = "fastp report")]
    pub report_title: String,

    // Sliding Window Quality Cutting
    #[arg(short='5', long="cut_front", action=ArgAction::SetTrue)]
    pub cut_front: bool,
    #[arg(short='3', long="cut_tail", action=ArgAction::SetTrue)]
    pub cut_tail: bool,
    #[arg(short='r', long="cut_right", action=ArgAction::SetTrue)]
    pub cut_right: bool,
    #[arg(short='W', long="cut_window_size", default_value_t=4)]
    pub cut_window_size: usize,
    #[arg(short='M', long="cut_mean_quality", default_value_t=20)]
    pub cut_mean_quality: u8,
    #[arg(long="cut_front_window_size", default_value_t=4)]
    pub cut_front_window_size: usize,
    #[arg(long="cut_front_mean_quality", default_value_t=20)]
    pub cut_front_mean_quality: u8,
    #[arg(long="cut_tail_window_size", default_value_t=4)]
    pub cut_tail_window_size: usize,
    #[arg(long="cut_tail_mean_quality", default_value_t=20)]
    pub cut_tail_mean_quality: u8,
    #[arg(long="cut_right_window_size", default_value_t=4)]
    pub cut_right_window_size: usize,
    #[arg(long="cut_right_mean_quality", default_value_t=20)]
    pub cut_right_mean_quality: u8,

    // Adapter Trimming
    #[arg(short='A', long="disable_adapter_trimming", action=ArgAction::SetTrue)]
    pub disable_adapter_trimming: bool,
    #[arg(short='a', long="adapter_sequence")]
    pub adapter_sequence: Option<String>,
    #[arg(long="adapter_sequence_r2")]
    pub adapter_sequence_r2: Option<String>,
    
    // PolyG Trimming
    #[arg(long="trim_poly_g", action=ArgAction::SetTrue)]
    pub trim_poly_g: bool,
    #[arg(long="poly_g_min_len", default_value_t=10)]
    pub poly_g_min_len: usize,
    #[arg(short='G', long="disable_trim_poly_g", action=ArgAction::SetTrue)]
    pub disable_trim_poly_g: bool,

    // PolyX Trimming
    #[arg(short='x', long="trim_poly_x", action=ArgAction::SetTrue)]
    pub trim_poly_x: bool,
    #[arg(long="poly_x_min_len", default_value_t=10)]
    pub poly_x_min_len: usize,

    // Overlap analysis and correction (PE)
    #[arg(short='c', long="correction", action=ArgAction::SetTrue)]
    pub correction: bool,
    #[arg(long="overlap_len_require", default_value_t=30)]
    pub overlap_len_require: usize,
    #[arg(long="overlap_diff_limit", default_value_t=5)]
    pub overlap_diff_limit: usize,
    #[arg(long="overlap_diff_percent_limit", default_value_t=20)]
    pub overlap_diff_percent_limit: u8,

    // Threading
    #[arg(short='w', long="thread", default_value_t=2)]
    pub thread: usize,

    // Performance tuning
    #[arg(long="pack_size", default_value_t=1000)]
    pub pack_size: usize,
    #[arg(long="queue_depth", default_value_t=0)]
    pub queue_depth: usize,
    #[arg(short='z', long="compression", default_value_t=4)]
    pub compression: u32,

    // External compressor (pigz)
    #[arg(long="pigz", action=ArgAction::SetTrue)]
    pub pigz: bool,
    #[arg(long="pigz_threads", default_value_t=0)]
    pub pigz_threads: usize,
}
