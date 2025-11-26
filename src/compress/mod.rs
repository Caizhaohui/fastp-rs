use crossbeam::channel::{unbounded, Sender, Receiver};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::Write;

pub struct CompressTask {
    pub id: u64,
    pub which: u8, // 1 for R1, 2 for R2
    pub data: Vec<u8>,
}

pub struct CompressedResult {
    pub id: u64,
    pub which: u8,
    pub data: Vec<u8>,
}

pub struct CompressionPool {
    tx: Sender<CompressTask>,
    pub rx: Receiver<CompressedResult>,
    level: u32,
}

impl CompressionPool {
    pub fn new(threads: usize, level: u32) -> Self {
        let (tx, worker_rx) = unbounded::<CompressTask>();
        let (worker_tx, rx) = unbounded::<CompressedResult>();
        for _ in 0..threads.max(1) {
            let rx = worker_rx.clone();
            let tx = worker_tx.clone();
            let lvl = level;
            std::thread::spawn(move || {
                while let Ok(task) = rx.recv() {
                    let mut enc = GzEncoder::new(Vec::new(), Compression::new(lvl));
                    // write all bytes
                    let _ = enc.write_all(&task.data);
                    let out = enc.finish().unwrap_or_default();
                    let _ = tx.send(CompressedResult { id: task.id, which: task.which, data: out });
                }
            });
        }
        Self { tx, rx, level }
    }

    pub fn submit(&self, id: u64, which: u8, data: Vec<u8>) {
        let _ = self.tx.send(CompressTask { id, which, data });
    }
}
