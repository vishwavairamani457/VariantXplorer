# FastQC Backend Code
import os
import json
import tempfile
import shutil
import gzip
import re
from collections import Counter, defaultdict
from statistics import mean
import os
import sys
import matplotlib.pyplot as plt
try:
    import pyfastx
    HAS_PYFASTX = True
except Exception:
    HAS_PYFASTX = False
from Bio import SeqIO

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# Create a temporary directory
def get_fixed_outdir():
    base_temp = tempfile.gettempdir()
    outdir = os.path.join(base_temp, "pyfastqc_out")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)
    return outdir

# Main Backend of FastQC
class FastQCBackend:
    def __init__(self):
        pass
    def run_fastqc(self, fastq_files, params, end_type="Single End", max_reads=500000000, progress_callback=None):
        if isinstance(fastq_files, str):
            fastq_files = [fastq_files]
        results, outdirs = {}, []
        if end_type == "Single End":
            result, outdir = self._run_single_file(fastq_files[0], params, max_reads, progress_callback)
            results = result
            outdirs.append(outdir)
        elif end_type == "Paired End":
            # run both R1 and R2 sequentially (could be parallelized later)
            result_r1, outdir_r1 = self._run_single_file(fastq_files[0], params, max_reads, progress_callback)
            result_r2, outdir_r2 = self._run_single_file(fastq_files[1], params, max_reads, progress_callback)
            results = {**{f"R1_{k}": v for k, v in result_r1.items()},
                       **{f"R2_{k}": v for k, v in result_r2.items()}}
            outdirs.extend([outdir_r1, outdir_r2])
        else:
            raise ValueError("Invalid end_type. Choose 'Single End' or 'Paired End'.")
        return results, outdirs

    #  Single File Run
    def _run_single_file(self, fastq_file, params, max_reads, progress_callback):
        outdir = get_fixed_outdir()
        qc_params = params.get("Quality Control", {})
        # Parameters
        phred_cutoff = int(qc_params.get("Phred Score", 0))
        min_len = int(qc_params.get("Min Read Length", 10))
        max_len = int(qc_params.get("Max Read Length", 1e9))
        # process all reads by default; set Sample Reads in params to limit
        sample_reads = qc_params.get("Sample Reads", None)
        if sample_reads in [None, "", "None", 0, "0"]:
            sample_reads = None
        else:
            sample_reads = int(sample_reads)
        # Load adapters/contaminants
        valid_bases = set("ACGTUN")
        user_adapters = self._split_sequences(qc_params.get("Input Adapter Sequence", ""))
        user_contaminants = self._split_sequences(qc_params.get("Contaminants/Primers Sequence", ""))
        default_adapters, default_contaminants = [], []
        default_adapters_file = resource_path("adapters.txt")
        default_contam_file = resource_path("contaminants.txt")
        if os.path.exists(default_adapters_file):
            with open(default_adapters_file) as f:
                default_adapters = [a.strip().upper() for a in f if a.strip() and set(a.strip().upper()) <= valid_bases]
        if os.path.exists(default_contam_file):
            with open(default_contam_file) as f:
                default_contaminants = [c.strip().upper() for c in f if c.strip() and set(c.strip().upper()) <= valid_bases]
        adapters = list(dict.fromkeys(default_adapters + user_adapters))
        contaminants = list(dict.fromkeys(default_contaminants + user_contaminants))
        # compile regex patterns (escape sequences to avoid regex meta-issues)
        adapter_patterns = {a: re.compile(re.escape(a)) for a in adapters if a}
        contaminant_patterns = {c: re.compile(re.escape(c)) for c in contaminants if c}
        # Containers
        read_lengths = []
        base_sums, base_counts = [], []
        gc_percentages = []
        adapter_hits, contaminant_hits = 0, 0
        adapter_freq, contaminant_freq = Counter(), Counter()
        adapter_examples, contaminant_examples = [], []
        #  Fastq Reading
        if progress_callback:
            progress_callback(15, f"Reading sequences from {os.path.basename(fastq_file)}...")
        def parse_reads():
            if HAS_PYFASTX:
                for name, seq, qual in pyfastx.Fastq(fastq_file, build_index=False):
                    yield seq, [ord(c) - 33 for c in qual]
            else:
                handle = gzip.open(fastq_file, "rt") if fastq_file.endswith(".gz") else open(fastq_file, "rt")
                for rec in SeqIO.parse(handle, "fastq-sanger"):
                    yield str(rec.seq), rec.letter_annotations.get("phred_quality", [])
                handle.close()
        # Processing Reads
        total_reads = sample_reads if sample_reads else 0
        if total_reads == 0 and HAS_PYFASTX:
            total_reads = len(pyfastx.Fastq(fastq_file, build_index=True))
        if total_reads == 0:
            total_reads = 1
        step = max(1, total_reads // 100)
        for i, (seq, quals) in enumerate(parse_reads()):
            if sample_reads and i >= sample_reads:
                break
            if not seq:
                continue
            seq = seq.upper()
            if len(seq) < min_len or len(seq) > max_len:
                continue
            if quals and mean(quals) < phred_cutoff:
                continue
            # basic stats
            read_lengths.append(len(seq))
            gc_pct = (seq.count("G") + seq.count("C")) / len(seq) * 100
            gc_percentages.append(gc_pct)
            # incremental per-base sums (avoid storing per-read arrays)
            for pos, q in enumerate(quals):
                if pos >= len(base_sums):
                    base_sums.append(q)
                    base_counts.append(1)
                else:
                    base_sums[pos] += q
                    base_counts[pos] += 1
            # adapters: check each adapter pattern, record freq & example
            for a_name, pattern in adapter_patterns.items():
                if pattern.search(seq):
                    adapter_hits += 1
                    adapter_freq[a_name] += 1
                    if len(adapter_examples) < 5:
                        adapter_examples.append(seq)
                    break
            # contaminants
            for c_name, pattern in contaminant_patterns.items():
                if pattern.search(seq):
                    contaminant_hits += 1
                    contaminant_freq[c_name] += 1
                    if len(contaminant_examples) < 5:
                        contaminant_examples.append(seq)
                    break
            # progress updates
            if progress_callback and (i % step == 0 or i == total_reads - 1):
                percent = 15 + int((i + 1) / total_reads * 75)
                progress_callback(percent, f"Processing reads {i+1}/{total_reads}")
        # Compile Results
        if not read_lengths:
            result = {"error": "No valid reads found!"}
            with open(os.path.join(outdir, "qc_results.json"), "w") as f:
                json.dump(result, f, indent=4)
            return result, outdir
        if progress_callback:
            progress_callback(90, "Compiling results...")
        per_base_quality = [{"Base": str(i + 1), "Mean": base_sums[i] / base_counts[i]} for i in range(len(base_sums))]
        mean_quality_scores = [p["Mean"] for p in per_base_quality]
        gc_mean = mean(gc_percentages) if gc_percentages else 0.0
        read_mean = mean(read_lengths)
        stats = {
            "Total Reads": len(read_lengths),
            "Mean Read Length": read_mean,
            "Mean GC Content": gc_mean,
            "Adapter Hits": adapter_hits,
            "Contaminant Hits": contaminant_hits,
        }
        result = {
            "basic_statistics": stats,
            "per_base_quality": per_base_quality,
            "mean_quality_scores": mean_quality_scores,
            "gc_mean": gc_mean,
            "length_distribution": read_lengths,
            "sequence_length_distribution": [{"Length": k, "Count": v} for k, v in sorted(Counter(read_lengths).items())],
            "per_sequence_gc_content": [{"GC": k, "Count": v} for k, v in sorted(Counter(int(round(g)) for g in gc_percentages).items())],
            "adapter_counts": dict(adapter_freq),
            "contaminant_counts": dict(contaminant_freq),
            "adapter_examples": adapter_examples,
            "contaminant_examples": contaminant_examples,
            "module_status": {
                "Basic Statistics": "pass" if len(read_lengths) > 0 else "fail",
                "Adapter Content": "warn" if adapter_hits > 0 else "pass",
                "Contaminants": "warn" if contaminant_hits > 0 else "pass",
            },
        }
        with open(os.path.join(outdir, "qc_results.json"), "w") as f:
            json.dump(result, f, indent=4)
        # Figures
        result["_figures"] = {
            "per_base_quality": self.plot_per_base_quality(per_base_quality),
            "gc_content": self.plot_gc_content(result["per_sequence_gc_content"]),
            "length_distribution": self.plot_length_distribution(result["sequence_length_distribution"]),
            "adapters_contaminants": self.plot_adapter_contaminants(stats),
        }
        if progress_callback:
            progress_callback(100, "QC Analysis Completed ✅")
        return result, outdir     
    # Helpers
    def _split_sequences(self, s):
        if isinstance(s, str):
            return [a.strip().upper() for a in s.split(",") if a.strip()]
        elif isinstance(s, list):
            return [str(a).upper().strip() for a in s if a]
        else:
            return []
    # Plots
    def plot_per_base_quality(self, per_base_quality):
        fig, ax = plt.subplots()
        if per_base_quality:
            bases = [int(p["Base"]) for p in per_base_quality]
            means = [p["Mean"] for p in per_base_quality]
            ax.axhspan(30, 60, facecolor="green", alpha=0.2)
            ax.axhspan(20, 30, facecolor="orange", alpha=0.2)
            ax.axhspan(0, 20, facecolor="red", alpha=0.2)
            ax.plot(bases, means, linewidth=1.2, label="Mean")
            ax.set_ylim(0, 60)
        ax.set_title("Per Base Sequence Quality")
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Phred Score")
        ax.legend(loc="upper right")
        return fig
    def plot_gc_content(self, gc_data):
        fig, ax = plt.subplots()
        if gc_data:
            x = [d["GC"] for d in gc_data]
            y = [d["Count"] for d in gc_data]
            ax.bar(x, y, color="skyblue", edgecolor="black")
        ax.set_title("GC Content Distribution")
        ax.set_xlabel("GC %")
        ax.set_ylabel("Count")
        return fig
    def plot_length_distribution(self, len_data):
        fig, ax = plt.subplots()
        if len_data:
            x = [d["Length"] for d in len_data]
            y = [d["Count"] for d in len_data]
            ax.bar(x, y, color="orange", edgecolor="black")
        ax.set_title("Sequence Length Distribution")
        ax.set_xlabel("Length")
        ax.set_ylabel("Count")
        return fig
    def plot_adapter_contaminants(self, stats):
        fig, ax = plt.subplots()
        labels = ["Adapter Hits", "Contaminant Hits"]
        values = [stats.get("Adapter Hits", 0), stats.get("Contaminant Hits", 0)]
        ax.bar(labels, values, color=["purple", "brown"], edgecolor="black")
        ax.set_title("Adapter and Contaminant Hits")
        return fig

# Interpretation
class FastQCInterpreter:
    def __init__(self, qc_results):
        self.qc_results = qc_results or {}
    def interpret_quality(self):
        scores = self.qc_results.get("mean_quality_scores", [])
        if not scores:
            return "No quality score data available."
        mean_q = sum(scores) / len(scores)
        if mean_q >= 30:
            return f"Per-base quality is good (average Q{mean_q:.2f} > Q30)."
        elif mean_q >= 20:
            return f"Per-base quality is moderate (average Q{mean_q:.2f}). Some trimming may help."
        else:
            return f"Per-base quality is poor (average Q{mean_q:.2f}). Trimming recommended."
    def interpret_gc(self):
        gc_mean = self.qc_results.get("gc_mean", None)
        if gc_mean is None:
            return "No GC content data available."
        if 40 <= gc_mean <= 60:
            return f"GC content ({gc_mean:.2f}%) is within expected range."
        else:
            return f"GC content ({gc_mean:.2f}%) is outside 40–60%. Possible contamination or bias."
    def interpret_adapters(self):
        adapter_counts = self.qc_results.get("adapter_counts", {}) or {}
        total = sum(adapter_counts.values())
        if total == 0:
            total = int(self.qc_results.get("basic_statistics", {}).get("Adapter Hits", 0))
        if total == 0:
            return "No adapter detected."
        details = ", ".join([f"{a}: {c} hits" for a, c in adapter_counts.items() if c > 0]) if adapter_counts else "No detailed adapter list available."
        if total < 10:
            return f"Low adapter presence ({total} hits). Details: {details}"
        else:
            return f"High adapter presence ({total} hits). Details: {details}"
    def interpret_contaminants(self):
        contaminant_counts = self.qc_results.get("contaminant_counts", {}) or {}
        total = sum(contaminant_counts.values())
        if total == 0:
            total = int(self.qc_results.get("basic_statistics", {}).get("Contaminant Hits", 0))
        if total == 0:
            return "No contaminants detected."
        details = ", ".join([f"{c}: {n} hits" for c, n in contaminant_counts.items() if n > 0]) if contaminant_counts else "No detailed contaminant list available."
        if total < 10:
            return f"Low contaminant presence ({total} hits). Details: {details}"
        else:
            return f"High contaminant presence ({total} hits). Details: {details}"
    def interpret_length_distribution(self):
        lengths = self.qc_results.get("length_distribution", [])
        if not lengths:
            return "No sequence length distribution data available."
        if len(set(lengths)) == 1:
            return f"Reads are uniform length ({lengths[0]} bp)."
        else:
            return f"Reads have variable lengths (range {min(lengths)}–{max(lengths)} bp)."
    def section_reports(self):
        return {
            "per_base_quality": self.interpret_quality(),
            "gc_content": self.interpret_gc(),
            "adapters": self.interpret_adapters(),
            "contaminants": self.interpret_contaminants(),
            "length_distribution": self.interpret_length_distribution(),
        }
