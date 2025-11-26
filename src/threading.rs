use crate::fastq::FastqRecord;
use crate::filter::Report;

pub struct Pack {
    pub id: u64,
    pub data: Vec<(FastqRecord, Option<FastqRecord>)>, // (R1, R2) - R2 is None for SE
}

pub struct ProcessedPack {
    pub id: u64,
    pub data: Vec<(FastqRecord, Option<FastqRecord>)>,
    pub report: Report,
}

// To maintain order, we can use a MinHeap or just a simple BTreeMap buffer in the writer
// But since packs come from workers, they might be out of order.
