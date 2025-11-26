use std::fs::File;
use std::io::{self, Write};
use crate::filter::Report;

pub fn write_html_report(path: &str, report: &Report, title: &str) -> io::Result<()> {
    let mut f = File::create(path)?;
    
    writeln!(f, "<!DOCTYPE html>")?;
    writeln!(f, "<html>")?;
    writeln!(f, "<head>")?;
    writeln!(f, "<title>{}</title>", title)?;
    writeln!(f, "<style>")?;
    writeln!(f, "body {{ font-family: Arial, sans-serif; margin: 20px; }}")?;
    writeln!(f, "table {{ border-collapse: collapse; width: 100%; max-width: 800px; }}")?;
    writeln!(f, "th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}")?;
    writeln!(f, "th {{ background-color: #f2f2f2; }}")?;
    writeln!(f, "h1 {{ color: #333; }}")?;
    writeln!(f, "</style>")?;
    writeln!(f, "</head>")?;
    writeln!(f, "<body>")?;
    
    writeln!(f, "<h1>{}</h1>", title)?;
    
    writeln!(f, "<h2>General Statistics</h2>")?;
    writeln!(f, "<table>")?;
    writeln!(f, "<tr><th>Metric</th><th>Value</th></tr>")?;
    writeln!(f, "<tr><td>Total Reads</td><td>{}</td></tr>", report.total_reads)?;
    writeln!(f, "<tr><td>Passed Reads</td><td>{}</td></tr>", report.passed_reads)?;
    writeln!(f, "<tr><td>Failed (Too Short)</td><td>{}</td></tr>", report.failed_too_short)?;
    writeln!(f, "<tr><td>Failed (Low Quality)</td><td>{}</td></tr>", report.failed_low_quality)?;
    writeln!(f, "<tr><td>Failed (Too many N)</td><td>{}</td></tr>", report.failed_n_excess)?;
    writeln!(f, "<tr><td>Failed (Low Avg Qual)</td><td>{}</td></tr>", report.failed_low_average_qual)?;
    writeln!(f, "</table>")?;

    writeln!(f, "<h2>Adapter Trimming</h2>")?;
    writeln!(f, "<table>")?;
    writeln!(f, "<tr><th>Metric</th><th>Value</th></tr>")?;
    writeln!(f, "<tr><td>Trimmed Reads</td><td>{}</td></tr>", report.adapter_trimmed_reads)?;
    writeln!(f, "<tr><td>Trimmed Bases</td><td>{}</td></tr>", report.adapter_trimmed_bases)?;
    writeln!(f, "</table>")?;

    writeln!(f, "<h2>PolyG Trimming</h2>")?;
    writeln!(f, "<table>")?;
    writeln!(f, "<tr><th>Metric</th><th>Value</th></tr>")?;
    writeln!(f, "<tr><td>Trimmed Reads</td><td>{}</td></tr>", report.poly_g_trimmed_reads)?;
    writeln!(f, "<tr><td>Trimmed Bases</td><td>{}</td></tr>", report.poly_g_trimmed_bases)?;
    writeln!(f, "</table>")?;

    writeln!(f, "<h2>PolyX Trimming</h2>")?;
    writeln!(f, "<table>")?;
    writeln!(f, "<tr><th>Metric</th><th>Value</th></tr>")?;
    writeln!(f, "<tr><td>Trimmed Reads</td><td>{}</td></tr>", report.poly_x_trimmed_reads)?;
    writeln!(f, "<tr><td>Trimmed Bases</td><td>{}</td></tr>", report.poly_x_trimmed_bases)?;
    writeln!(f, "</table>")?;

    writeln!(f, "<h2>PE Overlap Stats</h2>")?;
    writeln!(f, "<table>")?;
    writeln!(f, "<tr><th>Metric</th><th>Value</th></tr>")?;
    writeln!(f, "<tr><td>Overlap Pairs</td><td>{}</td></tr>", report.pe_overlap_count)?;
    writeln!(f, "<tr><td>Average Diff</td><td>{:.3}</td></tr>", report.pe_overlap_avg_diff)?;
    writeln!(f, "</table>")?;
    
    writeln!(f, "</body>")?;
    writeln!(f, "</html>")?;
    
    Ok(())
}
