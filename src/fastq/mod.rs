use std::io::{self, Read, Write, BufRead, BufReader};
use std::fs::File;
use flate2::write::GzEncoder;
use flate2::Compression;
use flate2::read::MultiGzDecoder;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub name: String,
    pub seq: String,
    pub plus: String,
    pub qual: String,
}

impl FastqRecord {
    pub fn new(name: String, seq: String, plus: String, qual: String) -> Self {
        Self { name, seq, plus, qual }
    }
}

pub struct Reader {
    reader: Box<dyn BufRead>,
}

impl Reader {
    pub fn new(path: Option<&str>, stdin: bool) -> io::Result<Self> {
        let reader: Box<dyn BufRead> = if stdin || path == Some("/dev/stdin") || path.is_none() {
            Box::new(BufReader::new(io::stdin()))
        } else {
            let p = path.unwrap();
            let f = File::open(p)?;
            if p.ends_with(".gz") {
                Box::new(BufReader::new(MultiGzDecoder::new(f)))
            } else {
                Box::new(BufReader::new(f))
            }
        };
        Ok(Self { reader })
    }

    pub fn next_record(&mut self) -> io::Result<Option<FastqRecord>> {
        // block buffer to reduce allocations
        let mut name = String::new();
        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();

        if self.reader.read_line(&mut name)? == 0 { return Ok(None); }
        if self.reader.read_line(&mut seq)? == 0 { return Ok(None); }
        if self.reader.read_line(&mut plus)? == 0 { return Ok(None); }
        if self.reader.read_line(&mut qual)? == 0 { return Ok(None); }

        if name.ends_with('\n') { name.pop(); }
        if name.ends_with('\r') { name.pop(); }
        if seq.ends_with('\n') { seq.pop(); }
        if seq.ends_with('\r') { seq.pop(); }
        if plus.ends_with('\n') { plus.pop(); }
        if plus.ends_with('\r') { plus.pop(); }
        if qual.ends_with('\n') { qual.pop(); }
        if qual.ends_with('\r') { qual.pop(); }

        Ok(Some(FastqRecord { name, seq, plus, qual }))
    }
}

pub struct Writer {
    writer: Box<dyn Write>,
}

impl Writer {
    pub fn new(path: Option<&str>, stdout: bool, compression_level: u32) -> io::Result<Self> {
        let writer: Box<dyn Write> = if stdout || path == Some("/dev/stdout") || path.is_none() {
            Box::new(io::stdout())
        } else {
            let p = path.unwrap();
            let f = File::create(p)?;
            if p.ends_with(".gz") {
                let enc = GzEncoder::new(f, Compression::new(compression_level));
                Box::new(enc)
            } else {
                Box::new(f)
            }
        };
        Ok(Self { writer })
    }

    pub fn write_record(&mut self, rec: &FastqRecord) -> io::Result<()> {
        use std::io::Write;
        self.writer.write_all(rec.name.as_bytes())?; self.writer.write_all(b"\n")?;
        self.writer.write_all(rec.seq.as_bytes())?;  self.writer.write_all(b"\n")?;
        self.writer.write_all(rec.plus.as_bytes())?; self.writer.write_all(b"\n")?;
        self.writer.write_all(rec.qual.as_bytes())?; self.writer.write_all(b"\n")?;
        Ok(())
    }
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> { self.writer.write(buf) }
    fn flush(&mut self) -> io::Result<()> { self.writer.flush() }
}
