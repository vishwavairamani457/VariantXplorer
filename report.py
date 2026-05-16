import os
import sys
import subprocess
import platform
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # ✅ Fix for PyInstaller - no display needed
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from datetime import datetime
import json
import tempfile

# ─────────────────────────────────────────
# DETECT IF RUNNING INSIDE PYWEBVIEW APP
# ─────────────────────────────────────────
def is_packaged():
    return getattr(sys, 'frozen', False)

def get_downloads_folder():
    """Get the user's Downloads folder cross-platform"""
    if platform.system() == "Windows":
        return os.path.join(os.environ["USERPROFILE"], "Downloads")
    elif platform.system() == "Darwin":
        return os.path.join(os.path.expanduser("~"), "Downloads")
    else:
        return os.path.join(os.path.expanduser("~"), "Downloads")

def open_file_in_os(filepath):
    """Open the file using the OS default viewer"""
    try:
        if platform.system() == "Windows":
            os.startfile(filepath)
        elif platform.system() == "Darwin":
            subprocess.Popen(["open", filepath])
        else:
            subprocess.Popen(["xdg-open", filepath])
    except Exception as e:
        print(f"Could not open file: {e}")

# ─────────────────────────────────────────
# SAFE IMPORTS
# ─────────────────────────────────────────
def safe_import(func_name, module_name):
    try:
        module = __import__(module_name)
        return getattr(module, func_name)
    except:
        return None

get_fastqc_summary  = safe_import("get_fastqc_summary",  "fastqc_backend")
get_trimming_stats  = safe_import("get_trimming_stats",   "trimming")
get_alignment_stats = safe_import("get_alignment_stats",  "bwa_backend")
get_variant_count   = safe_import("get_variant_count",    "freebayes")
get_annotation_summary = safe_import("get_annotation_summary", "annotation")
get_filtering_results  = safe_import("get_filtering_results",  "filtering")

# ─────────────────────────────────────────
# SAFE WRAPPERS (unchanged from your code)
# ─────────────────────────────────────────
def get_fastqc_summary_from_json(qc_json_path):
    if not qc_json_path or not os.path.exists(qc_json_path):
        return {
            "Total Reads": "NA", "Mean Read Length": "NA",
            "Mean GC Content": "NA", "GC Interpretation": "NA",
            "Adapter Hits": "NA", "Contaminants": "NA", "Per Base Quality": "NA"
        }
    try:
        with open(qc_json_path) as f:
            data = json.load(f)
        stats = data.get("basic_statistics", {})
        gc_mean = stats.get("Mean GC Content", 0)
        gc_interp = f"GC content ({gc_mean:.2f}%) is normal" if 40 <= gc_mean <= 60 else f"GC content ({gc_mean:.2f}%) deviates"
        mean_q = sum(data.get("mean_quality_scores", [])) / max(len(data.get("mean_quality_scores", [])), 1)
        quality_interp = "High quality (Q30+)" if mean_q >= 30 else ("Moderate quality" if mean_q >= 20 else "Poor quality")
        return {
            "Total Reads": stats.get("Total Reads", "NA"),
            "Mean Read Length": round(stats.get("Mean Read Length", 0), 2),
            "Mean GC Content": f"{round(gc_mean, 2)}%",
            "GC Interpretation": gc_interp,
            "Adapter Hits": stats.get("Adapter Hits", 0),
            "Contaminants": stats.get("Contaminant Hits", 0),
            "Per Base Quality": quality_interp
        }
    except:
        return None

def get_trimming_summary():
    try:
        import streamlit as st
        if "trimming_stats" in st.session_state:
            stats = st.session_state["trimming_stats"]
            return {
                "Total Reads": stats.get("total_reads", "NA"),
                "Trimmed Reads": stats.get("trimmed_reads", "NA"),
                "Dropped Reads": stats.get("dropped_reads", "NA"),
                "Adapter Trimmed": stats.get("adapter_hits_trimmed", "NA"),
                "Contaminant Trimmed": stats.get("contaminant_hits_trimmed", "NA"),
                "Mean Read Length": stats.get("mean_read_length", "NA"),
                "Mean GC Content": f"{stats.get('mean_gc_content', 'NA')}%"
            }
    except:
        pass
    return None

def get_alignment_summary_from_sam(temp_dir):
    try:
        import glob
        sam_files = glob.glob(os.path.join(temp_dir, "*.sam"))
        if not sam_files:
            raise FileNotFoundError("No SAM file found")
        total, mapped = 0, 0
        with open(sam_files[0]) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                total += 1
                if not (int(line.split("\t")[1]) & 4):
                    mapped += 1
        return {
            "Total Reads": total,
            "Mapped Reads": mapped,
            "Mapping Rate": f"{(mapped / total) * 100:.2f}%" if total > 0 else "0%"
        }
    except:
        return {"Total Reads": "NA", "Mapped Reads": "NA", "Mapping Rate": "NA"}

def safe_alignment(temp_dir):
    try:
        if get_alignment_stats:
            data = get_alignment_stats(temp_dir)
            if data and isinstance(data, dict):
                if "Total Reads" in data and "Mapped Reads" in data:
                    total = int(data["Total Reads"])
                    mapped = int(data["Mapped Reads"])
                    data["Mapping Rate"] = f"{(mapped / total) * 100:.2f}%"
                return data
        return get_alignment_summary_from_sam(temp_dir)
    except:
        return {"Total Reads": "NA", "Mapped Reads": "NA", "Mapping Rate": "NA"}

def get_annotation_summary_from_df(temp_dir):
    try:
        import streamlit as st
        if "annotation_df" in st.session_state:
            df = st.session_state["annotation_df"]
        else:
            import glob
            vcf_files = glob.glob(os.path.join(temp_dir, "*annotated*.vcf"))
            if not vcf_files:
                raise FileNotFoundError("No annotated VCF found")
            impacts = []
            with open(vcf_files[0]) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    if "ANN=" in line:
                        ann = line.split("ANN=")[1].split(";")[0]
                        parts = ann.split("|")
                        if len(parts) >= 2:
                            impacts.append(parts[1])
            if not impacts:
                raise ValueError("No annotation data")
            from collections import Counter
            counts = Counter(impacts)
            return {
                "High Impact": counts.get("HIGH", 0),
                "Moderate Impact": counts.get("MODERATE", 0),
                "Low Impact": counts.get("LOW", 0)
            }
        if "Impact" not in df.columns:
            raise ValueError("Impact column missing")
        return {
            "High Impact": int((df["Impact"] == "HIGH").sum()),
            "Moderate Impact": int((df["Impact"] == "MODERATE").sum()),
            "Low Impact": int((df["Impact"] == "LOW").sum())
        }
    except:
        return {"High Impact": "NA", "Moderate Impact": "NA", "Low Impact": "NA"}

def safe_annotation(temp_dir):
    try:
        if get_annotation_summary:
            data = get_annotation_summary(temp_dir)
            if data and isinstance(data, dict):
                return data
        return get_annotation_summary_from_df(temp_dir)
    except:
        return {"High Impact": "NA", "Moderate Impact": "NA", "Low Impact": "NA"}

def get_filtering_summary(filtered_df):
    if filtered_df is None or filtered_df.empty:
        return {"Total Passed Variants": 0, "High Impact": 0, "Moderate Impact": 0, "Low Impact": 0}
    impact_col = next((col for col in filtered_df.columns if col.lower() == "impact"), None)
    if impact_col is None:
        return {"Total Passed Variants": len(filtered_df), "High Impact": 0, "Moderate Impact": 0, "Low Impact": len(filtered_df)}
    return {
        "Total Passed Variants": len(filtered_df),
        "High Impact": int((filtered_df[impact_col] == "HIGH").sum()),
        "Moderate Impact": int((filtered_df[impact_col] == "MODERATE").sum()),
        "Low Impact": int((filtered_df[impact_col] == "LOW").sum()),
    }

# ─────────────────────────────────────────
# VCF PARSER
# ─────────────────────────────────────────
def parse_vcf(vcf_file):
    variants = []
    if not os.path.exists(vcf_file):
        return pd.DataFrame()
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt, _, _, info = parts[:8]
            dp, af = "NA", "NA"
            for field in info.split(";"):
                if field.startswith("DP="):
                    dp = field.split("=")[1]
                elif field.startswith("AF="):
                    af = field.split("=")[1]
            vtype = "SNP" if len(ref) == 1 and len(alt) == 1 else ("INS" if len(ref) < len(alt) else "DEL")
            variants.append([chrom, pos, ref, alt, vtype, dp, af])
    return pd.DataFrame(variants, columns=["CHROM", "POS", "REF", "ALT", "TYPE", "DP", "AF"])

# ─────────────────────────────────────────
# VARIANT PLOT
# ─────────────────────────────────────────
def generate_plot(df, temp_dir):
    if df.empty:
        return None
    counts = df["TYPE"].value_counts()
    plt.figure()
    counts.plot(kind="bar")
    plt.title("Variant Distribution")
    plot_path = os.path.join(temp_dir, "variant_plot.png")
    plt.savefig(plot_path)
    plt.close()
    return plot_path

def find_fastqc_json(temp_dir):
    path1 = os.path.join(temp_dir, "qc_results.json")
    path2 = os.path.join(tempfile.gettempdir(), "pyfastqc_out", "qc_results.json")
    return path1 if os.path.exists(path1) else (path2 if os.path.exists(path2) else None)

# ─────────────────────────────────────────
# MAIN DATA COLLECTION
# ─────────────────────────────────────────
def collect_pipeline_data(temp_dir, vcf_file):
    qc_json_path = find_fastqc_json(temp_dir)
    fastqc    = get_fastqc_summary_from_json(qc_json_path) if qc_json_path else None
    trimming  = get_trimming_summary()
    alignment = safe_alignment(temp_dir)
    annotation = safe_annotation(temp_dir)
    total_variants = get_variant_count(vcf_file) if get_variant_count else len(parse_vcf(vcf_file))
    filtered_df = None
    if get_filtering_results:
        filtered_df = get_filtering_results(temp_dir)
    if filtered_df is None:
        try:
            import streamlit as st
            filtered_df = st.session_state.get("filtered_df", None)
        except:
            filtered_df = None
    return fastqc, trimming, alignment, total_variants, annotation, filtered_df

# ─────────────────────────────────────────
# REPORT GENERATION
# ─────────────────────────────────────────
def generate_report(output_pdf, temp_dir, vcf_file, logo_path=None):
    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(output_pdf, pagesize=A4)
    elements = []
    fastqc, trimming, alignment, total_variants, annotation, filtered_df = collect_pipeline_data(temp_dir, vcf_file)
    variant_df = parse_vcf(vcf_file)
    plot_path = generate_plot(variant_df, temp_dir)

    # HEADER
    if logo_path and os.path.exists(logo_path):
        elements.append(Image(logo_path, width=80, height=50))
    elements.append(Paragraph("<b>VARIANT ANALYSIS REPORT</b>", styles['Title']))
    elements.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles['Normal']))
    elements.append(Spacer(1, 20))

    # FASTQC
    elements.append(Paragraph("<b>FastQC Summary</b>", styles['Heading2']))
    if fastqc:
        elements.append(Table(list(fastqc.items())))
    else:
        elements.append(Paragraph("No FastQC data available", styles['Normal']))
    elements.append(Spacer(1, 15))

    # TRIMMING
    elements.append(Paragraph("<b>Trimming Summary</b>", styles['Heading2']))
    if trimming:
        elements.append(Table(list(trimming.items())))
    else:
        elements.append(Paragraph("No trimming performed", styles['Normal']))
    elements.append(Spacer(1, 15))

    # ALIGNMENT
    elements.append(Paragraph("<b>Alignment Summary (BWA)</b>", styles['Heading2']))
    elements.append(Table(list(alignment.items())))
    elements.append(Spacer(1, 15))

    # VARIANT COUNT
    elements.append(Paragraph("<b>Variant Calling Summary</b>", styles['Heading2']))
    elements.append(Table([["Total Variants", total_variants]]))
    elements.append(Spacer(1, 15))

    # ANNOTATION
    elements.append(Paragraph("<b>Annotation Summary</b>", styles['Heading2']))
    elements.append(Table(list(annotation.items())))
    elements.append(Spacer(1, 15))

    # FILTERED VARIANTS
    elements.append(Paragraph("<b>Filtered Important Variants</b>", styles['Heading2']))
    try:
        import streamlit as st
        filtered_df = st.session_state.get("annotated_df")
        if isinstance(filtered_df, pd.DataFrame) and not filtered_df.empty:
            df = filtered_df.copy()
            df.columns = [c.strip() for c in df.columns]
            if "Impact" in df.columns:
                df = df[df["Impact"].isin(["HIGH", "MODERATE"])]
            selected_cols = ["Start", "Ref", "Alt", "Gene", "HGVS.c", "HGVS.p", "VariantName"]
            available_cols = [c for c in selected_cols if c in df.columns]
            clean_df = df[available_cols].head(20) if available_cols else df.head(20)
            if not clean_df.empty:
                rename_map = {
                    "Start": "Position", "Ref": "Ref", "Alt": "Alt",
                    "Gene": "Gene", "HGVS.c": "cDNA Change",
                    "HGVS.p": "Protein Change", "VariantName": "Variant"
                }
                clean_df.rename(columns=rename_map, inplace=True)
                table_data = [clean_df.columns.tolist()]
                for row in clean_df.values:
                    table_data.append([str(x)[:40] for x in row])
                table = Table(table_data, repeatRows=1)
                table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), colors.darkgreen),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.black),
                    ("FONTSIZE", (0, 0), (-1, -1), 8),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ]))
                elements.append(table)
            else:
                elements.append(Paragraph("No high/moderate impact variants found", styles['Normal']))
        else:
            elements.append(Paragraph("No filtering results found", styles['Normal']))
    except:
        elements.append(Paragraph("No filtering results found", styles['Normal']))

    # PLOT
    if plot_path and os.path.exists(plot_path):
        elements.append(Paragraph("<b>Variant Distribution</b>", styles['Heading2']))
        elements.append(Image(plot_path, width=400, height=300))

    # FOOTER
    elements.append(Spacer(1, 30))
    elements.append(Paragraph("Generated by VariantXplorer", styles['Italic']))
    doc.build(elements)

# ─────────────────────────────────────────
# ✅ MAIN RUN FUNCTION — FIXED FOR APP
# ─────────────────────────────────────────
def run_report(temp_dir, vcf_file):
    # ✅ Always save to Downloads folder so user can find it
    downloads = get_downloads_folder()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_pdf = os.path.join(downloads, f"VariantXplorer_Report_{timestamp}.pdf")

    generate_report(
        output_pdf=output_pdf,
        temp_dir=temp_dir,
        vcf_file=vcf_file,
        logo_path="logo.png"
    )

    # ✅ Auto-open the PDF after generation
    open_file_in_os(output_pdf)

    return output_pdf