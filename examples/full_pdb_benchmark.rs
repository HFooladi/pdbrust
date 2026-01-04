//! Comprehensive PDB Archive Benchmark
//!
//! Benchmarks PDBRust against the full PDB archive of gzip-compressed files.
//! This example tests parsing correctness, measures performance, and generates
//! detailed reports for quality assurance before public release.
//!
//! # Usage
//!
//! ```bash
//! cargo run --release --example full_pdb_benchmark \
//!     --features "gzip,parallel,descriptors,quality,summary" \
//!     -- /path/to/pdb/archive --output-dir ./benchmark_results
//! ```
//!
//! # Output Files
//!
//! - `benchmark_report.txt` - Comprehensive statistics and timing
//! - `failures.tsv` - Detailed log of all parsing failures
//! - `timing_histogram.txt` - Parse time distribution

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use pdbrust::{PdbError, PdbStructure, parse_gzip_pdb_file};

// ============== Error Classification ==============

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub enum ErrorCategory {
    IoError,
    InvalidRecord,
    ParseError,
    EmptyStructure,
}

impl std::fmt::Display for ErrorCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ErrorCategory::IoError => write!(f, "IoError"),
            ErrorCategory::InvalidRecord => write!(f, "InvalidRecord"),
            ErrorCategory::ParseError => write!(f, "ParseError"),
            ErrorCategory::EmptyStructure => write!(f, "EmptyStructure"),
        }
    }
}

impl From<&PdbError> for ErrorCategory {
    fn from(err: &PdbError) -> Self {
        match err {
            PdbError::IoError(_) => ErrorCategory::IoError,
            PdbError::InvalidRecord(_) => ErrorCategory::InvalidRecord,
            PdbError::ParseError(_) => ErrorCategory::ParseError,
        }
    }
}

// ============== Statistics Tracking ==============

#[derive(Debug, Default)]
pub struct BenchmarkStats {
    pub total_files: usize,
    pub successful_parses: usize,
    pub failed_parses: usize,

    // Timing statistics (in microseconds for precision)
    pub parse_times_us: Vec<u64>,

    // Error breakdown
    pub error_counts: HashMap<ErrorCategory, usize>,
    pub error_examples: HashMap<ErrorCategory, Vec<String>>,

    // Structure statistics
    pub total_atoms: u64,
    pub total_residues: u64,
    pub atom_counts: Vec<usize>,
    pub largest_structure: (String, usize), // (pdb_id, num_atoms)
    pub smallest_structure: (String, usize),

    // Analysis statistics (if enabled)
    #[cfg(feature = "summary")]
    pub analysis_times_us: Vec<u64>,
    #[cfg(feature = "summary")]
    pub successful_analyses: usize,
}

#[derive(Debug, Default)]
pub struct TimingStats {
    pub count: usize,
    pub mean_ms: f64,
    pub median_ms: f64,
    pub p95_ms: f64,
    pub p99_ms: f64,
    pub min_ms: f64,
    pub max_ms: f64,
    pub std_dev_ms: f64,
    pub total_seconds: f64,
}

impl TimingStats {
    fn from_microseconds(times_us: &[u64]) -> Self {
        if times_us.is_empty() {
            return Self::default();
        }

        let mut sorted: Vec<u64> = times_us.to_vec();
        sorted.sort_unstable();

        let count = sorted.len();
        let sum: u64 = sorted.iter().sum();
        let mean = sum as f64 / count as f64;

        let median_idx = count / 2;
        let median = sorted[median_idx] as f64;

        let p95_idx = ((count as f64) * 0.95) as usize;
        let p95 = sorted[p95_idx.min(count - 1)] as f64;

        let p99_idx = ((count as f64) * 0.99) as usize;
        let p99 = sorted[p99_idx.min(count - 1)] as f64;

        let min = sorted[0] as f64;
        let max = sorted[count - 1] as f64;

        let variance: f64 = sorted
            .iter()
            .map(|&t| (t as f64 - mean).powi(2))
            .sum::<f64>()
            / count as f64;
        let std_dev = variance.sqrt();

        Self {
            count,
            mean_ms: mean / 1000.0,
            median_ms: median / 1000.0,
            p95_ms: p95 / 1000.0,
            p99_ms: p99 / 1000.0,
            min_ms: min / 1000.0,
            max_ms: max / 1000.0,
            std_dev_ms: std_dev / 1000.0,
            total_seconds: sum as f64 / 1_000_000.0,
        }
    }
}

// ============== Benchmark Runner ==============

pub struct BenchmarkRunner {
    input_dir: PathBuf,
    output_dir: PathBuf,
    progress_interval: usize,
    run_analysis: bool,
}

impl BenchmarkRunner {
    pub fn new(input_dir: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            input_dir,
            output_dir,
            progress_interval: 10000,
            run_analysis: cfg!(feature = "summary"),
        }
    }

    /// Discover all .ent.gz files in the PDB archive directory structure
    pub fn discover_files(&self) -> Vec<PathBuf> {
        let mut files = Vec::new();

        // PDB archive structure: pdb/{two-char-code}/pdb{code}.ent.gz
        if let Ok(entries) = fs::read_dir(&self.input_dir) {
            for entry in entries.filter_map(|e| e.ok()) {
                let subdir = entry.path();
                if subdir.is_dir() {
                    if let Ok(subentries) = fs::read_dir(&subdir) {
                        for subentry in subentries.filter_map(|e| e.ok()) {
                            let path = subentry.path();
                            if let Some(ext) = path.extension() {
                                if ext == "gz" {
                                    files.push(path);
                                }
                            }
                        }
                    }
                }
            }
        }

        files.sort();
        files
    }

    /// Extract PDB ID from filename (e.g., pdb1ubq.ent.gz -> 1ubq)
    fn extract_pdb_id(path: &Path) -> String {
        path.file_stem()
            .and_then(|s| s.to_str())
            .map(|s| s.trim_start_matches("pdb").trim_end_matches(".ent"))
            .unwrap_or("unknown")
            .to_string()
    }

    /// Run the full benchmark
    pub fn run(&self) -> BenchmarkStats {
        println!("╔══════════════════════════════════════════════════════════════╗");
        println!("║           PDBRust Full Archive Benchmark                     ║");
        println!("╚══════════════════════════════════════════════════════════════╝\n");

        // Phase 1: File Discovery
        println!("Phase 1: Discovering files...");
        let files = self.discover_files();
        let total_files = files.len();
        println!("  Found {} files to process\n", total_files);

        if total_files == 0 {
            eprintln!(
                "ERROR: No .ent.gz files found in {}",
                self.input_dir.display()
            );
            return BenchmarkStats::default();
        }

        // Create output directory
        fs::create_dir_all(&self.output_dir).expect("Failed to create output directory");

        // Initialize thread-safe counters
        let processed = AtomicUsize::new(0);
        let success_count = AtomicUsize::new(0);
        let fail_count = AtomicUsize::new(0);
        let start_time = Instant::now();

        // Initialize results collectors
        let stats = Mutex::new(BenchmarkStats {
            total_files,
            smallest_structure: (String::new(), usize::MAX),
            ..Default::default()
        });
        let failures: Mutex<Vec<(String, String, String)>> = Mutex::new(Vec::new());

        // Phase 2: Parallel Processing
        let num_threads = rayon::current_num_threads();
        println!(
            "Phase 2: Processing {} files with {} threads...\n",
            total_files, num_threads
        );

        files.par_iter().for_each(|path| {
            let pdb_id = Self::extract_pdb_id(path);

            // Parse
            let parse_start = Instant::now();
            let parse_result = parse_gzip_pdb_file(path);
            let parse_duration = parse_start.elapsed();

            // Update progress atomically
            let current = processed.fetch_add(1, Ordering::Relaxed) + 1;
            if current % self.progress_interval == 0 || current == total_files {
                let elapsed = start_time.elapsed().as_secs_f64();
                let rate = current as f64 / elapsed;
                let eta = if rate > 0.0 {
                    (total_files - current) as f64 / rate
                } else {
                    0.0
                };
                let success = success_count.load(Ordering::Relaxed);
                let failed = fail_count.load(Ordering::Relaxed);
                println!(
                    "  Progress: {:>7}/{} ({:>5.1}%) | Success: {:>6} | Failed: {:>5} | Rate: {:>6.0}/s | ETA: {:>5.0}s",
                    current,
                    total_files,
                    100.0 * current as f64 / total_files as f64,
                    success,
                    failed,
                    rate,
                    eta
                );
            }

            match parse_result {
                Ok(structure) => {
                    let num_atoms = structure.atoms.len();

                    if num_atoms == 0 {
                        // Empty structure - record as failure
                        fail_count.fetch_add(1, Ordering::Relaxed);
                        let mut s = stats.lock().unwrap();
                        s.failed_parses += 1;
                        *s.error_counts
                            .entry(ErrorCategory::EmptyStructure)
                            .or_insert(0) += 1;
                        let examples = s
                            .error_examples
                            .entry(ErrorCategory::EmptyStructure)
                            .or_default();
                        if examples.len() < 10 {
                            examples.push(pdb_id.clone());
                        }

                        let mut f = failures.lock().unwrap();
                        f.push((pdb_id, "EmptyStructure".to_string(), "No atoms parsed".to_string()));
                        return;
                    }

                    success_count.fetch_add(1, Ordering::Relaxed);

                    // Get residue count
                    let num_residues = Self::count_residues(&structure);

                    // Update stats
                    {
                        let mut s = stats.lock().unwrap();
                        s.successful_parses += 1;
                        s.parse_times_us.push(parse_duration.as_micros() as u64);
                        s.total_atoms += num_atoms as u64;
                        s.total_residues += num_residues as u64;
                        s.atom_counts.push(num_atoms);

                        // Track largest/smallest
                        if num_atoms > s.largest_structure.1 {
                            s.largest_structure = (pdb_id.clone(), num_atoms);
                        }
                        if num_atoms < s.smallest_structure.1 {
                            s.smallest_structure = (pdb_id.clone(), num_atoms);
                        }
                    }

                    // Run analysis if enabled
                    #[cfg(feature = "summary")]
                    if self.run_analysis {
                        let analysis_start = Instant::now();
                        let _summary = structure.summary();
                        let analysis_duration = analysis_start.elapsed();

                        let mut s = stats.lock().unwrap();
                        s.successful_analyses += 1;
                        s.analysis_times_us.push(analysis_duration.as_micros() as u64);
                    }
                }
                Err(e) => {
                    fail_count.fetch_add(1, Ordering::Relaxed);
                    let category = ErrorCategory::from(&e);
                    let message = format!("{}", e);

                    let mut s = stats.lock().unwrap();
                    s.failed_parses += 1;
                    *s.error_counts.entry(category.clone()).or_insert(0) += 1;

                    let examples = s
                        .error_examples
                        .entry(category.clone())
                        .or_default();
                    if examples.len() < 10 {
                        examples.push(format!("{}: {}", pdb_id, message.lines().next().unwrap_or("")));
                    }

                    let mut f = failures.lock().unwrap();
                    f.push((pdb_id, category.to_string(), message));
                }
            }
        });

        let total_time = start_time.elapsed();

        // Extract results
        let stats = stats.into_inner().unwrap();
        let failures = failures.into_inner().unwrap();

        // Phase 3: Write reports
        println!("\nPhase 3: Writing results...");
        self.write_failure_log(&failures);
        self.write_summary_report(&stats, total_time);
        self.write_timing_histogram(&stats);

        println!("\n  Results written to: {}", self.output_dir.display());

        stats
    }

    fn count_residues(structure: &PdbStructure) -> usize {
        use std::collections::HashSet;
        let mut seen = HashSet::new();
        for atom in &structure.atoms {
            seen.insert((atom.chain_id.clone(), atom.residue_seq, atom.ins_code));
        }
        seen.len()
    }

    fn write_failure_log(&self, failures: &[(String, String, String)]) {
        let path = self.output_dir.join("failures.tsv");
        let file = File::create(&path).expect("Failed to create failures.tsv");
        let mut writer = BufWriter::new(file);

        writeln!(writer, "pdb_id\terror_category\terror_message").unwrap();
        for (pdb_id, category, message) in failures {
            // Escape tabs and newlines in error message
            let clean_message = message.replace(['\t', '\n'], " ");
            writeln!(writer, "{}\t{}\t{}", pdb_id, category, clean_message).unwrap();
        }

        println!("  - Failure log: {}", path.display());
    }

    fn write_summary_report(&self, stats: &BenchmarkStats, total_time: Duration) {
        let path = self.output_dir.join("benchmark_report.txt");
        let file = File::create(&path).expect("Failed to create report");
        let mut w = BufWriter::new(file);

        let parse_timing = TimingStats::from_microseconds(&stats.parse_times_us);

        writeln!(
            w,
            "╔══════════════════════════════════════════════════════════════╗"
        )
        .unwrap();
        writeln!(
            w,
            "║           PDBRust Full Archive Benchmark Report              ║"
        )
        .unwrap();
        writeln!(
            w,
            "╚══════════════════════════════════════════════════════════════╝"
        )
        .unwrap();
        writeln!(w).unwrap();
        writeln!(w, "Total Execution Time: {:.2}s", total_time.as_secs_f64()).unwrap();
        writeln!(w).unwrap();

        // File Processing Summary
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "FILE PROCESSING SUMMARY").unwrap();
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "Total Files:        {:>12}", stats.total_files).unwrap();
        writeln!(
            w,
            "Successful Parses:  {:>12} ({:.4}%)",
            stats.successful_parses,
            100.0 * stats.successful_parses as f64 / stats.total_files as f64
        )
        .unwrap();
        writeln!(
            w,
            "Failed Parses:      {:>12} ({:.4}%)",
            stats.failed_parses,
            100.0 * stats.failed_parses as f64 / stats.total_files as f64
        )
        .unwrap();
        writeln!(w).unwrap();

        // Parsing Timing
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "PARSING TIMING STATISTICS").unwrap();
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(
            w,
            "Total Parse Time:   {:>12.2}s",
            parse_timing.total_seconds
        )
        .unwrap();
        writeln!(w, "Mean:               {:>12.3}ms", parse_timing.mean_ms).unwrap();
        writeln!(w, "Median:             {:>12.3}ms", parse_timing.median_ms).unwrap();
        writeln!(w, "P95:                {:>12.3}ms", parse_timing.p95_ms).unwrap();
        writeln!(w, "P99:                {:>12.3}ms", parse_timing.p99_ms).unwrap();
        writeln!(w, "Min:                {:>12.3}ms", parse_timing.min_ms).unwrap();
        writeln!(w, "Max:                {:>12.3}ms", parse_timing.max_ms).unwrap();
        writeln!(w, "Std Dev:            {:>12.3}ms", parse_timing.std_dev_ms).unwrap();
        if parse_timing.total_seconds > 0.0 {
            writeln!(
                w,
                "Throughput:         {:>12.0} files/sec",
                stats.successful_parses as f64 / parse_timing.total_seconds
            )
            .unwrap();
        }
        writeln!(w).unwrap();

        // Analysis Timing (if enabled)
        #[cfg(feature = "summary")]
        {
            let analysis_timing = TimingStats::from_microseconds(&stats.analysis_times_us);
            if analysis_timing.count > 0 {
                writeln!(
                    w,
                    "══════════════════════════════════════════════════════════════"
                )
                .unwrap();
                writeln!(w, "ANALYSIS TIMING STATISTICS").unwrap();
                writeln!(
                    w,
                    "══════════════════════════════════════════════════════════════"
                )
                .unwrap();
                writeln!(w, "Structures Analyzed: {:>11}", stats.successful_analyses).unwrap();
                writeln!(
                    w,
                    "Total Analysis Time: {:>11.2}s",
                    analysis_timing.total_seconds
                )
                .unwrap();
                writeln!(
                    w,
                    "Mean:                {:>11.3}ms",
                    analysis_timing.mean_ms
                )
                .unwrap();
                writeln!(
                    w,
                    "Median:              {:>11.3}ms",
                    analysis_timing.median_ms
                )
                .unwrap();
                writeln!(w, "P99:                 {:>11.3}ms", analysis_timing.p99_ms).unwrap();
                writeln!(w).unwrap();
            }
        }

        // Error Breakdown
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "ERROR BREAKDOWN").unwrap();
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        if stats.error_counts.is_empty() {
            writeln!(w, "No errors encountered!").unwrap();
        } else {
            let mut errors: Vec<_> = stats.error_counts.iter().collect();
            errors.sort_by(|a, b| b.1.cmp(a.1));

            for (category, count) in errors {
                writeln!(
                    w,
                    "{:20} {:>8} ({:.4}%)",
                    format!("{:?}:", category),
                    count,
                    100.0 * *count as f64 / stats.total_files as f64
                )
                .unwrap();
                if let Some(examples) = stats.error_examples.get(category) {
                    for ex in examples.iter().take(3) {
                        writeln!(w, "  └─ {}", ex).unwrap();
                    }
                }
            }
        }
        writeln!(w).unwrap();

        // Structure Statistics
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "STRUCTURE STATISTICS").unwrap();
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "Total Atoms Parsed:     {:>15}", stats.total_atoms).unwrap();
        writeln!(w, "Total Residues Parsed:  {:>15}", stats.total_residues).unwrap();
        if stats.successful_parses > 0 {
            writeln!(
                w,
                "Average Atoms/Structure:    {:>11.0}",
                stats.total_atoms as f64 / stats.successful_parses as f64
            )
            .unwrap();
            writeln!(
                w,
                "Average Residues/Structure: {:>11.0}",
                stats.total_residues as f64 / stats.successful_parses as f64
            )
            .unwrap();
        }
        writeln!(
            w,
            "Largest Structure:  {} ({} atoms)",
            stats.largest_structure.0, stats.largest_structure.1
        )
        .unwrap();
        if stats.smallest_structure.1 < usize::MAX {
            writeln!(
                w,
                "Smallest Structure: {} ({} atoms)",
                stats.smallest_structure.0, stats.smallest_structure.1
            )
            .unwrap();
        }
        writeln!(w).unwrap();

        // Memory Estimates
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "MEMORY ESTIMATES").unwrap();
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        let bytes_per_atom = 200; // Rough estimate for Atom struct + overhead
        let peak_memory_mb = (stats.largest_structure.1 * bytes_per_atom) as f64 / 1e6;
        writeln!(w, "Est. Peak Memory (largest): {:>8.1} MB", peak_memory_mb).unwrap();
        writeln!(w).unwrap();

        // Version Info
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "ENVIRONMENT").unwrap();
        writeln!(
            w,
            "══════════════════════════════════════════════════════════════"
        )
        .unwrap();
        writeln!(w, "PDBRust Version: {}", env!("CARGO_PKG_VERSION")).unwrap();
        writeln!(w, "Rayon Threads:   {}", rayon::current_num_threads()).unwrap();

        println!("  - Benchmark report: {}", path.display());
    }

    fn write_timing_histogram(&self, stats: &BenchmarkStats) {
        if stats.parse_times_us.is_empty() {
            return;
        }

        let path = self.output_dir.join("timing_histogram.txt");
        let file = File::create(&path).expect("Failed to create histogram");
        let mut w = BufWriter::new(file);

        writeln!(w, "Parse Time Distribution (milliseconds)").unwrap();
        writeln!(w, "========================================").unwrap();

        // Create histogram buckets (log scale)
        let buckets = [
            0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0,
        ];
        let mut counts = vec![0usize; buckets.len() + 1];

        for &time_us in &stats.parse_times_us {
            let time_ms = time_us as f64 / 1000.0;
            let bucket_idx = buckets
                .iter()
                .position(|&b| time_ms < b)
                .unwrap_or(buckets.len());
            counts[bucket_idx] += 1;
        }

        let max_count = *counts.iter().max().unwrap_or(&1);
        let bar_width = 50;

        // Print histogram
        let mut prev_bound = 0.0;
        for (i, &count) in counts.iter().enumerate() {
            let bound = if i < buckets.len() {
                buckets[i]
            } else {
                f64::INFINITY
            };

            let bar_len = if max_count > 0 {
                count * bar_width / max_count
            } else {
                0
            };
            let bar: String = "█".repeat(bar_len);

            if bound.is_infinite() {
                writeln!(w, "{:>8.1}ms+   {:>8} │{}", prev_bound, count, bar).unwrap();
            } else {
                writeln!(
                    w,
                    "{:>6.1}-{:<6.1}ms {:>8} │{}",
                    prev_bound, bound, count, bar
                )
                .unwrap();
            }
            prev_bound = bound;
        }

        println!("  - Timing histogram: {}", path.display());
    }
}

// ============== Main ==============

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 || args.iter().any(|a| a == "--help" || a == "-h") {
        println!("PDBRust Full Archive Benchmark");
        println!();
        println!("Usage: {} <pdb_archive_dir> [OPTIONS]", args[0]);
        println!();
        println!("Arguments:");
        println!("  <pdb_archive_dir>    Path to PDB archive directory containing .ent.gz files");
        println!();
        println!("Options:");
        println!(
            "  --output-dir <dir>   Output directory for results (default: ./benchmark_results)"
        );
        println!("  --help, -h           Show this help message");
        println!();
        println!("Example:");
        println!("  cargo run --release --example full_pdb_benchmark \\");
        println!("      --features \"gzip,parallel,descriptors,quality,summary\" \\");
        println!("      -- /data/shared/comp3d/databases/pdb_data/pdb/ --output-dir ./results");
        std::process::exit(0);
    }

    let input_dir = PathBuf::from(&args[1]);
    let output_dir = args
        .iter()
        .position(|a| a == "--output-dir")
        .and_then(|i| args.get(i + 1))
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("./benchmark_results"));

    if !input_dir.exists() {
        eprintln!(
            "ERROR: Input directory does not exist: {}",
            input_dir.display()
        );
        std::process::exit(1);
    }

    let runner = BenchmarkRunner::new(input_dir, output_dir);
    let stats = runner.run();

    // Print final summary
    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                      FINAL SUMMARY                           ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!(
        "Success Rate: {}/{} ({:.4}%)",
        stats.successful_parses,
        stats.total_files,
        if stats.total_files > 0 {
            100.0 * stats.successful_parses as f64 / stats.total_files as f64
        } else {
            0.0
        }
    );
    println!("Failed: {}", stats.failed_parses);

    let timing = TimingStats::from_microseconds(&stats.parse_times_us);
    if timing.count > 0 {
        println!(
            "Parse Timing: mean={:.3}ms, median={:.3}ms, p99={:.3}ms",
            timing.mean_ms, timing.median_ms, timing.p99_ms
        );
        println!(
            "Throughput: {:.0} files/sec",
            if timing.total_seconds > 0.0 {
                stats.successful_parses as f64 / timing.total_seconds
            } else {
                0.0
            }
        );
    }

    // Exit with error code if too many failures
    if stats.failed_parses > 0 {
        let failure_rate = stats.failed_parses as f64 / stats.total_files as f64;
        if failure_rate > 0.05 {
            eprintln!(
                "\nWARNING: High failure rate ({:.2}%) - check failures.tsv for details",
                failure_rate * 100.0
            );
            std::process::exit(1);
        }
    }
}
