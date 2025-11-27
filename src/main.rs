mod config;
mod fastq;
mod filter;
mod threading;
mod html_report;
mod compress;

use clap::Parser;
use std::io;
use std::fs::File;
use std::thread;
use std::sync::{Arc, Mutex};
use std::collections::BinaryHeap;
use std::cmp::Ordering;
use crossbeam::channel::{bounded, Sender, Receiver};

use crate::config::Cli;
use crate::fastq::{Reader, Writer, FastqRecord};
use crate::filter::{Filter, Report};
use crate::threading::{Pack, ProcessedPack};
use crate::html_report::write_html_report;
use crate::compress::CompressionPool;
use std::io::Write as IoWrite;
use clap::Subcommand;
// use serde::Serialize;

// Helper for ordering ProcessedPack in BinaryHeap (MinHeap)
struct OrderedPack(ProcessedPack);

impl PartialEq for OrderedPack {
    fn eq(&self, other: &Self) -> bool {
        self.0.id == other.0.id
    }
}
impl Eq for OrderedPack {}
impl PartialOrd for OrderedPack {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Reverse order for MinHeap
        other.0.id.partial_cmp(&self.0.id)
    }
}
impl Ord for OrderedPack {
    fn cmp(&self, other: &Self) -> Ordering {
        other.0.id.cmp(&self.0.id)
    }
}

#[derive(Subcommand, Debug, Clone)]
enum Cmd {
    EmitSbatch {
        #[arg(long="partition", default_value = "<your_partition>")]
        partition: String,
        #[arg(long="cpus", default_value_t = 24)]
        cpus: usize,
        #[arg(long="mem", default_value = "64G")]
        mem: String,
    },
    Run,
}

fn main() -> io::Result<()> {
    let mut cli = Cli::parse();
    // simple subcommand via env var FASTP_RS_CMD, to avoid extra clap changes to Cli
    if let Ok(cmd) = std::env::var("FASTP_RS_CMD") {
        match cmd.as_str() {
            "emit_sbatch" => {
                let partition = std::env::var("FASTP_RS_PARTITION").unwrap_or_else(|_| "<your_partition>".to_string());
                let cpus: usize = std::env::var("FASTP_RS_CPUS").ok().and_then(|v| v.parse().ok()).unwrap_or(24);
                let mem = std::env::var("FASTP_RS_MEM").unwrap_or_else(|_| "64G".to_string());
                let script = format!("#!/usr/bin/env bash\n#SBATCH -p {partition}\n#SBATCH -c {cpus}\n#SBATCH --mem={mem}\n#SBATCH -J fastp_rs_pool\n#SBATCH -o fastp_rs.%j.out\n#SBATCH -e fastp_rs.%j.err\ncd \"$SLURM_SUBMIT_DIR\"\n\nP={cpus}\nPACK={pack}\nQD={qd}\nZ={z}\n\ncargo build --release\n/usr/bin/time -v ./target/release/fastp_rs \\\n  -i R1.fq.gz -I R2.fq.gz \\\n  -o out1.fq.gz -O out2.fq.gz \\\n  --json fastp.json --html fastp.html \\\n  -w $P --pack_size $PACK --queue_depth $QD -z $Z \\\n  -x --poly_x_min_len 10 \\\n  --trim_poly_g --poly_g_min_len 10 \\\n  -c --overlap_len_require 30 \\\n     --overlap_diff_limit 5 \\\n     --overlap_diff_percent_limit 20\n",
                    partition=partition,
                    cpus=cpus,
                    mem=mem,
                    pack=cli.pack_size,
                    qd=if cli.queue_depth == 0 { (if cli.thread==0 { num_cpus::get() } else { cli.thread })*4 } else { cli.queue_depth },
                    z=cli.compression,
                );
                println!("{}", script);
                return Ok(());
            },
            _ => {}
        }
    }
    let thread_num = if cli.thread == 0 { num_cpus::get() } else { cli.thread };
    let pack_size = if cli.pack_size == 0 { 1000 } else { cli.pack_size };

    let qd = if cli.queue_depth == 0 { thread_num * 2 } else { cli.queue_depth };
    let (tx_pack, rx_pack): (Sender<Pack>, Receiver<Pack>) = bounded(qd);
    let (tx_out, rx_out): (Sender<ProcessedPack>, Receiver<ProcessedPack>) = bounded(qd);

    let filter = Arc::new(Filter::new(cli.clone()));
    
    // 1. Workers
    let mut workers = Vec::new();
    for _ in 0..thread_num {
        let rx = rx_pack.clone();
        let tx = tx_out.clone();
        let filter = filter.clone();
        
        let handle = thread::spawn(move || {
            while let Ok(pack) = rx.recv() {
                let mut processed_data = Vec::with_capacity(pack.data.len());
                let mut local_report = Report::default();
                
                for (r1, r2_opt) in pack.data {
                    local_report.total_reads += 1;
                    
                    if let Some(r2) = r2_opt {
                        // PE Processing
                        let (rec1, rec2) = filter.trim_pair(r1, r2, &mut local_report);
                        
                        if filter.pass_filters(&rec1, &mut local_report) && filter.pass_filters(&rec2, &mut local_report) {
                            local_report.passed_reads += 1;
                            processed_data.push((rec1, Some(rec2)));
                        }
                    } else {
                        // SE Processing
                        let rec1 = filter.trim_record(r1, false, &mut local_report);
                        if filter.pass_filters(&rec1, &mut local_report) {
                            local_report.passed_reads += 1;
                            processed_data.push((rec1, None));
                        }
                    }
                }
                
                tx.send(ProcessedPack {
                    id: pack.id,
                    data: processed_data,
                    report: local_report,
                }).unwrap();
            }
        });
        workers.push(handle);
    }
    
    // Drop original tx_out so receiver closes when workers finish
    drop(tx_out);

    // 2. Writer Thread
    let cli_writer = cli.clone();
    let final_report = Arc::new(Mutex::new(Report::default()));
    let final_report_clone = final_report.clone();
    
    let writer_handle = thread::spawn(move || -> io::Result<()> {
        // compression pool for .gz outputs when not using external pigz
        let use_pool = (!cli_writer.pigz) && (cli_writer.out1.as_deref().map(|p| p.ends_with(".gz")).unwrap_or(false)
                       || cli_writer.out2.as_deref().map(|p| p.ends_with(".gz")).unwrap_or(false));
        let pool = if use_pool { Some(CompressionPool::new(thread_num, cli_writer.compression)) } else { None };
        let mut w1 = if let Some(p) = &cli_writer.out1 {
            if cli_writer.pigz && p.ends_with(".gz") {
                Some(Writer::new(Some("/dev/stdout"), cli_writer.stdout, cli_writer.compression)?)
            } else {
                Some(Writer::new(Some(p), cli_writer.stdout, cli_writer.compression)?)
            }
        } else if cli_writer.stdout {
             Some(Writer::new(None, true, cli_writer.compression)?)
        } else {
            None
        };
        
        let mut w2 = if let Some(p) = &cli_writer.out2 {
            if cli_writer.pigz && p.ends_with(".gz") {
                Some(Writer::new(Some("/dev/stdout"), false, cli_writer.compression)?)
            } else {
                Some(Writer::new(Some(p), false, cli_writer.compression)?)
            }
        } else {
            None
        };

        let mut next_id = 0;
        let mut buffer = BinaryHeap::new();

        for pack in rx_out {
            buffer.push(OrderedPack(pack));
            
            while let Some(top) = buffer.peek() {
                if top.0.id == next_id {
                    let OrderedPack(p) = buffer.pop().unwrap();
                    
                    // Aggregate Report
                    {
                        let mut rep = final_report_clone.lock().unwrap();
                        rep.merge(&p.report);
                    }
                    
                    // Write Output
                    for (r1, r2_opt) in p.data {
                        // write or submit to compression pool
                        if let Some(pool) = &pool {
                            if let Some(path) = &cli_writer.out1 { if path.ends_with(".gz") {
                                let mut buf = Vec::new();
                                use std::fmt::Write as FmtWrite;
                                write!(&mut buf, "{}\n{}\n{}\n{}\n", r1.name, r1.seq, r1.plus, r1.qual).unwrap();
                                pool.submit(p.id, 1, buf);
                            } else if let Some(w) = &mut w1 { w.write_record(&r1)?; } }
                            if let Some(r2) = r2_opt {
                                if let Some(path) = &cli_writer.out2 { if path.ends_with(".gz") {
                                    let mut buf = Vec::new();
                                    use std::fmt::Write as FmtWrite;
                                    write!(&mut buf, "{}\n{}\n{}\n{}\n", r2.name, r2.seq, r2.plus, r2.qual).unwrap();
                                    pool.submit(p.id, 2, buf);
                                } else if let Some(w) = &mut w2 { w.write_record(&r2)?; } }
                            }
                        } else {
                            if let Some(w) = &mut w1 { w.write_record(&r1)?; }
                            if let Some(r2) = r2_opt { if let Some(w) = &mut w2 { w.write_record(&r2)?; } }
                        }
                    }
                    // drain compressed results in order
                    if let Some(pool) = &pool {
                        use std::collections::BTreeMap;
                        let mut bucket: BTreeMap<(u64,u8), Vec<u8>> = BTreeMap::new();
                        // drain up to current pack id
                        while let Ok(res) = pool.rx.try_recv() {
                            bucket.insert((res.id, res.which), res.data);
                        }
                        // write sorted by (id, which)
                        for ((pid, which), data) in bucket {
                            if pid == p.id {
                                if which == 1 { if let Some(w) = &mut w1 { IoWrite::write_all(w, &data)?; } }
                                else { if let Some(w) = &mut w2 { IoWrite::write_all(w, &data)?; } }
                            }
                        }
                    }
                    next_id += 1;
                } else {
                    break;
                }
            }
        }
        Ok(())
    });

    // 3. Reader (Main Thread)
    let mut pack_data = Vec::with_capacity(pack_size);
    let mut pack_id = 0;

    if cli.in1.is_some() && cli.in2.is_some() {
        // PE
        let mut r1 = Reader::new(cli.in1.as_deref(), cli.stdin)?;
        let mut r2 = Reader::new(cli.in2.as_deref(), false)?;
        
        loop {
            let rec1_opt = r1.next_record()?;
            let rec2_opt = r2.next_record()?;
            
            match (rec1_opt, rec2_opt) {
                (Some(rec1), Some(rec2)) => {
                    pack_data.push((rec1, Some(rec2)));
                    if pack_data.len() >= pack_size {
                        tx_pack.send(Pack { id: pack_id, data: pack_data }).unwrap();
                        pack_data = Vec::with_capacity(pack_size);
                        pack_id += 1;
                    }
                },
                (None, None) => break,
                _ => {
                    eprintln!("Error: PE input files have different number of reads");
                    break;
                }
            }
        }
    } else {
        // SE
        let mut r1 = Reader::new(cli.in1.as_deref(), cli.stdin)?;
        loop {
            match r1.next_record()? {
                Some(rec1) => {
                    pack_data.push((rec1, None));
                    if pack_data.len() >= pack_size {
                        tx_pack.send(Pack { id: pack_id, data: pack_data }).unwrap();
                        pack_data = Vec::with_capacity(pack_size);
                        pack_id += 1;
                    }
                },
                None => break,
            }
        }
    }
    
    // Send remaining data
    if !pack_data.is_empty() {
        tx_pack.send(Pack { id: pack_id, data: pack_data }).unwrap();
    }
    
    // Close pack channel to notify workers
    drop(tx_pack);
    
    // Wait for workers
    for w in workers {
        w.join().unwrap();
    }
    
    // Wait for writer
    writer_handle.join().unwrap()?;

    // Generate JSON Report
    let jf = cli.json;
    let mut jf_w = File::create(jf)?;
    let rep = final_report.lock().unwrap();
    serde_json::to_writer_pretty(&mut jf_w, &*rep)?;

    // Generate HTML Report
    let hf = cli.html;
    if !hf.is_empty() {
        write_html_report(&hf, &*rep, &cli.report_title)?;
    }
    
    Ok(())
}
