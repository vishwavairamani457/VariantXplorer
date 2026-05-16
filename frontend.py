# Main Frontend Code
import numpy as np
import streamlit as st
import base64
from fastqc_backend import FastQCBackend, FastQCInterpreter
import matplotlib.pyplot as plt
import tempfile, os, shutil, json, time



import os
import sys

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))

    return os.path.join(base_path, relative_path)
# Set the font as Times New Roman
st.markdown(
    """
    <style>
    /* Global font override to Times New Roman */
    html, body, [class*="css"] {
        font-family: "Times New Roman", Times, serif !important;
        color: white !important;
    }
    h1, h2, h3, h4, h5, h6, p, div, label, input, textarea, span, button {
        font-family: "Times New Roman", Times, serif !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)
# page config
st.set_page_config(page_title="VariantXplorer", layout="wide")
# convert image to base64
def get_base64_image(image_path):
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode()
    except Exception:
        return None
bg_image = get_base64_image(resource_path("geneimage.jpg"))
# Custom CSS for background
if bg_image:
    st.markdown(
        f"""
        <style>
        .stApp {{
            height: 100%;
            margin: 0;
            padding: 0;
            background-image: url("data:image/jpeg;base64,{bg_image}");
            background-size: cover;
            background-repeat: no-repeat;
            background-position: center center;
            background-attachment: fixed;
            image-rendering: auto;
        }}
        .content {{
            color: white;
            text-align: center;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )
# Display Title
st.markdown('<div class="content"><h1> VariantXplorer - An Integrated Variant Analysis Tool', unsafe_allow_html=True)
# Upload box commands
st.markdown("""
    <style>
    section[data-testid="stFileUploader"] label {
        color: white !important;
    }
    /* This targets the label of file uploader */
    .stFileUploader > label {
        color: white !important;
    }
    </style>
""", unsafe_allow_html=True)
# Create columns to position the uploader on the right
col1, col2, col3 = st.columns([1, 1, 2])
with col3:
    fastq_file = st.file_uploader(
        "Choose FASTQ file(s) (single-end: 1 file, paired-end: 2 files)",
        type=["fastq"],
        accept_multiple_files=True
    )
if fastq_file:
    with col3:
        for f in fastq_file:
            st.success(f"FASTQ uploaded: {f.name}")
with col3:
    reference_file = st.file_uploader(
        "Upload Reference Genome (FASTA) and GTF/GFF file together",
        type=["fasta", "gtf", "gff", "gff3"],
        accept_multiple_files=True
    )
if reference_file:
    with col3:
        for f in reference_file:
            st.success(f"Reference/Annotation uploaded: {f.name}")
# Split reference_files into FASTA and GFF/GTF
fasta_file = None
gff_file = None
if reference_file:
    for f in reference_file:
        fname = f.name.lower()
        if fname.endswith((".fa", ".fasta", ".fna")):
            fasta_file = f
        elif fname.endswith((".gtf", ".gff", ".gff3")):
            gff_file = f
    if fasta_file is None:
        st.error("Reference FASTA file is required")
        st.stop()
    if gff_file is None:
        st.error("No GTF/GFF found — intron/exon annotations disabled")
        st.stop()
# Session state initialization / safer button handling
if "run_clicked" not in st.session_state:
    st.session_state.run_clicked = False
if "uploaded_paths" not in st.session_state:
    st.session_state.uploaded_paths = None
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = None
if "reference_path" not in st.session_state:
    st.session_state.reference_path = None
if "fastqc_history" not in st.session_state:
    st.session_state.fastqc_history = []  # list of dicts: {"results":..., "outdirs":..., "timestamp":...}
if "trimming_history" not in st.session_state:
    st.session_state.trimming_history = []  # list of dicts: {"paths":..., "stats":..., "timestamp":...}
if "fastqc_on_trim_history" not in st.session_state:
    st.session_state.fastqc_on_trim_history = []  # list of dicts similar to fastqc_history
if "bwa_history" not in st.session_state:
    st.session_state.bwa_history = []  # list of dicts: {"sam_path":..., "df":..., "timestamp":...}
if "fastqc_done" not in st.session_state:
    st.session_state.fastqc_done = False
if "do_trimming" not in st.session_state:
    st.session_state.do_trimming = None 
if "trimming_done" not in st.session_state:
    st.session_state.trimming_done = False
if "fastqc_on_trim_done" not in st.session_state:
    st.session_state.fastqc_on_trim_done = False
if "bwa_done" not in st.session_state:
    st.session_state.bwa_done = False
if "bwa_counter" not in st.session_state:
    st.session_state.bwa_counter = 0
if "end_type" not in st.session_state:
    st.session_state.end_type = "Single End"
# persistent cached display data
if "cached_sam_df" not in st.session_state:
    st.session_state.cached_sam_df = None          # DataFrame of SAM records
if "cached_vcf_df" not in st.session_state:
    st.session_state.cached_vcf_df = None          # DataFrame of VCF records
if "cached_annotated_df" not in st.session_state:
    st.session_state.cached_annotated_df = None    # annotated variants DataFrame
if "cached_filtered_df" not in st.session_state:
    st.session_state.cached_filtered_df = None     # filtered variants DataFrame
if "cached_trimming_stats" not in st.session_state:
    st.session_state.cached_trimming_stats = {}    # {label: stats}
# Create Run Analysis button that sets session state rather than local variable
col1, col2, col3 = st.columns([1, 1, 2])
with col3:
    if st.button("Run Analysis"):
        st.session_state.run_clicked = True
    st.markdown(
        "<p style='color:white;'>Click here to start analysis with default parameters or change the parameter based on your requirements !!!</p>",
        unsafe_allow_html=True
    )
# Add readend type
col1, col2, col3, col4 = st.columns(4)
with col3:
    st.markdown(
        "<p style='color:white;'>Select the readend type (Paired end / Single end)</p>",
        unsafe_allow_html=True)  
if "end_type" not in st.session_state:
    st.session_state.end_type = "Single End"
col1, col2, col3, col4 = st.columns(4)
with col3:
    if st.button("Paired End", type="primary" if st.session_state.end_type == "Paired End" else "secondary"):
        st.session_state.end_type = "Paired End"
with col4:
    if st.button("Single End", type="primary" if st.session_state.end_type == "Single End" else "secondary"):
        st.session_state.end_type = "Single End"
end_type = st.session_state.end_type
if "tool_params" in st.session_state:
    st.session_state.tool_params["Trimming"]["Paired or Single End"] = end_type
# Initialize tool params and defaults
tools = ["Quality Control", "Trimming", "Reference Mapping", "Variant Calling", "Filter"]
if "default_params" not in st.session_state:
    st.session_state.default_params = {
        "Quality Control": {
            "Phred Score": 0,
            "Min Read Length": 10,
            "Max Read Length": 999999999999,
            "Input Adapter Sequence": "",
            "Contaminants/Primers Sequence": ""
        },
        "Trimming": {
            "Window Size": 5,
            "Window Qual": 10,
            "Leading Quality": 3,
            "Trailing Quality": 3,
            "Min Read Length": 10,
            "Read Crop Length": 0,
            "Read Head Crop Length": 0,
            "Min GC Content": 40,
            "Max GC Content": 60,
            "ns_max_p  (Percentage of Ambiguous Nucleotide)": 1,
            "ns_max_n  (No.of Ambiguous Nucleotide)": 3,
            "lc_threshold (Low Complexity Sequence)": 10
        },
        "Reference Mapping": {
            "Mismatch Penalty": 4,
            "Threads": 8,
            "Seed length": 20,
            "Band Width": 100,
            "Z-dropoff": 100,
            "Seed Split ratio": 1.5,
            "Max Occurence": 10000,
            "Match Score": 1,
            "Gap Open Penalty":6,
            "Gap Extension Penalty":1,
            "Clipping Penalty":5,
            "Min Map Score": 1,
            "LowComplexityThreshold": 0.8,
        },
        "Variant Calling": {            
            "ploidy": 1,
            "min_mapping_quality": 0,
            "min_base_quality": 0,
            "reference_mapping_quality": 100,
            "reference_base_quality": 60,
            "min_coverage": 1,
            "min_alt_count": 1,
            "min_alt_fraction": 0.001,
            "min_depth": 0,
            "max_haplotypes": 12,
            "max_haplotype_length": 500,
            "min_haplotype_support": 1,
            "theta": 0.001,
            "pvar": 0.0001,
            "read_dependence_factor": 0.9,
            "genotyping_max_iter": 500,
            "site_selection_max_iter": 5,
            "indel_exclusion_window": 5,
            "debug": False,
        },
        "Filter": {
            "Min Depth": 0,
            "Max Allele Frequency": 1.0,
            "Allowed Impacts": ["HIGH", "MODERATE"]            
        }
    }
if "tool_params" not in st.session_state or not isinstance(st.session_state.tool_params, dict):
    st.session_state.tool_params = st.session_state.default_params.copy()
else:
    # Update defaults only for missing tools, don't overwrite user-edited values
    for tool, defaults in st.session_state.default_params.items():
        if tool not in st.session_state.tool_params:
            st.session_state.tool_params[tool] = defaults.copy()
if "show_dialog" not in st.session_state:
    st.session_state.show_dialog = None
st.markdown(
    """
    <h2 id="tools-parameter">Tools Parameter</h2>
    """,
    unsafe_allow_html=True
)
st.markdown(
    """
    <style>
    #tools-parameter {
        color: white !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)
if "show_dialog" not in st.session_state:
    st.session_state.show_dialog = None
if "tool_params" not in st.session_state:
    st.session_state.tool_params = st.session_state.default_params.copy()
cols = st.columns(len(tools))
for i, tool in enumerate(tools):
    if cols[i].button(tool, key=f"btn_{tool}"):
        st.session_state.show_dialog = tool
        st.rerun()     
# OPEN PARAMETER DIALOG
def open_param_dialog(tool_name):
    @st.dialog(f"{tool_name} Parameters")
    def dialog_content():
        st.write(f"Edit parameters for **{tool_name}**")
        # Render editable parameters
        for key, val in st.session_state.tool_params[tool_name].items():
            unique_key = f"{tool_name}_{key}"
            if isinstance(val, bool):
                new_val = st.checkbox(key, value=val, key=unique_key)
            elif isinstance(val, int):
                new_val = st.number_input(key, value=val, key=unique_key)
            elif isinstance(val, float):
                new_val = st.number_input(key, value=val, format="%.5f", key=unique_key)
            else:
                new_val = st.text_input(key, value=str(val), key=unique_key)
            st.session_state.tool_params[tool_name][key] = new_val
        # CLOSE BUTTON
        if st.button("Close", key=f"close_{tool_name}"):
            st.session_state.show_dialog = None
            st.rerun()
    dialog_content()
if st.session_state.show_dialog:
    open_param_dialog(st.session_state.show_dialog)
# Prepare uploads on initial Run click
if st.session_state.run_clicked and st.session_state.uploaded_paths is None:
    if not fastq_file:
        st.error("Please upload at least one FASTQ file first.")
        st.session_state.run_clicked = False
    else:
        tmp_dir = tempfile.mkdtemp()
        st.session_state.temp_dir = tmp_dir
        uploaded_paths = []
        for f in fastq_file:
            f.seek(0)
            dst = os.path.join(tmp_dir, f.name)
            with open(dst, "wb") as out_f:
                out_f.write(f.read())
            uploaded_paths.append(dst)
        st.session_state.uploaded_paths = uploaded_paths
        if reference_file:
            ref_dst = os.path.join(tmp_dir, "reference.fasta")
            fasta_file.seek(0)
            if gff_file:
                gff_file.seek(0)            
            with open(ref_dst, "wb") as rf:
                fasta_file.seek(0)
                rf.write(fasta_file.read())                
            st.session_state.reference_path = ref_dst
        if not st.session_state.end_type:
            st.session_state.end_type = "Paired End" if len(uploaded_paths) == 2 else "Single End"
# HELPER: render a single FastQC result dict
def _render_fastqc_single(result, label_prefix=""):
    st.markdown(f"<h4 style='color:white;'>{label_prefix}Basic Statistics</h4>", unsafe_allow_html=True)
    st.json(result["basic_statistics"])
    st.markdown(f"<h4 style='color:white;'>{label_prefix}Module Status</h4>", unsafe_allow_html=True)
    st.json(result["module_status"])
    interpreter = FastQCInterpreter({
        "mean_quality_scores": [p["Mean"] for p in result["per_base_quality"]],
        "gc_mean": result["basic_statistics"]["Mean GC Content"],
        "adapter_counts": result.get("adapter_counts", {}),
        "contaminant_counts": result.get("contaminant_counts", {}),
        "length_distribution": [d["Length"] for d in result["sequence_length_distribution"]],
    })
    interpretations = interpreter.section_reports()
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.pyplot(result["_figures"]["per_base_quality"])
        st.markdown(f"<p style='color:white; font-size:14px'>{interpretations['per_base_quality']}</p>", unsafe_allow_html=True)
    with col2:
        st.pyplot(result["_figures"]["gc_content"])
        st.markdown(f"<p style='color:white; font-size:14px'>{interpretations['gc_content']}</p>", unsafe_allow_html=True)
    with col3:
        st.pyplot(result["_figures"]["length_distribution"])
        st.markdown(f"<p style='color:white; font-size:14px'>{interpretations['length_distribution']}</p>", unsafe_allow_html=True)
    with col4:
        st.pyplot(result["_figures"]["adapters_contaminants"])
        st.markdown(f"<p style='color:white; font-size:14px'>{interpretations['adapters']}</p>", unsafe_allow_html=True)
        st.markdown(f"<p style='color:white; font-size:14px'>{interpretations['contaminants']}</p>", unsafe_allow_html=True)
TABS_CSS = """
<style>
.stTabs [role="tab"] p { color: white !important; }
.stTabs [role="tab"][aria-selected="true"] p { color: white !important; font-weight: bold; }
.stTabs [role="tab"][aria-selected="true"] { border-bottom: 3px solid white !important; }
</style>
"""
# Backend Call Code for FASTQC 
if st.session_state.run_clicked and not st.session_state.fastqc_done and st.session_state.uploaded_paths:
    backend = FastQCBackend()
    progress_text = "🔄 Running QC Analysis..."
    my_bar = st.progress(0, text="")
    progress_placeholder = st.empty()
    status_placeholder = st.empty()
    progress_placeholder.markdown(
        f"<h4 style='color:white;'>{progress_text} 0%</h4>",
        unsafe_allow_html=True
    )
    status_placeholder.markdown(
        f"<p style='color:white;'>Step: Initializing...</p>",
        unsafe_allow_html=True
    )
    def update_progress(percent, step_name):
        my_bar.progress(percent)
        progress_placeholder.markdown(
            f"<h4 style='color:white;'>{progress_text} {percent}%</h4>",
            unsafe_allow_html=True
        )
        status_placeholder.markdown(
            f"<p style='color:white;'>⚙️ {step_name}... {percent}%</p>",
            unsafe_allow_html=True
        )
    try:
        results, outdirs = backend.run_fastqc(
            st.session_state.uploaded_paths,
            st.session_state.tool_params,
            end_type=st.session_state.end_type,
            progress_callback=update_progress
        )
    except Exception as e:
        my_bar.empty()
        status_placeholder.empty()
        progress_placeholder.empty()
        st.error(f"FastQC failed: {e}")
        st.session_state.run_clicked = False
    else:
        my_bar.empty()
        status_placeholder.empty()
        progress_placeholder.markdown("<h3 style='color:white;'>✅ QC Analysis Completed</h3>", unsafe_allow_html=True)
        st.session_state.fastqc_history.append({
            "results": results,
            "outdirs": outdirs,
            "timestamp": time.strftime("%Y%m%d-%H%M%S")
        })
        st.session_state.fastqc_done = True
# Always display every FastQC run from history
if st.session_state.fastqc_history:
    for idx, entry in enumerate(st.session_state.fastqc_history):
        results = entry["results"]
        st.markdown(f"<h3 style='color:white;'>✅ FastQC Result</h3>", unsafe_allow_html=True)
        if st.session_state.end_type == "Single End":
            _render_fastqc_single(results)
        else:
            r1_result = {k[len("R1_"):]: v for k, v in results.items() if k.startswith("R1_")}
            r2_result = {k[len("R2_"):]: v for k, v in results.items() if k.startswith("R2_")}
            tab_r1, tab_r2 = st.tabs(["Read 1 (R1)", "Read 2 (R2)"])
            st.markdown(TABS_CSS, unsafe_allow_html=True)
            for tab, result, label in zip([tab_r1, tab_r2], [r1_result, r2_result], ["R1", "R2"]):
                with tab:
                    _render_fastqc_single(result, label_prefix=f"{label} - ")

# Backend Call Code for Trimming
if st.session_state.fastqc_done and st.session_state.do_trimming is None:
    st.markdown("<h3 style='color:white;text-align:center;'>FastQC Completed: Do you want to run trimming before alignment?</h3>", unsafe_allow_html=True)
    col1,col2= st.columns([1, 1])
    with col1:
        if st.button("Yes, Perform Trimming", key="trim_yes", use_container_width=True):
            st.session_state.do_trimming = True
    with col2:
        if st.button("No, Skip Trimming", key="trim_no", use_container_width=True):
            st.session_state.do_trimming = False


# Backend Call Code for Trimming
from trimming import trim_reads, read_fastq, get_trimming_params
def load_list_from_file(filename):
    if os.path.exists(filename):
        with open(filename) as f:
            return [line.strip().upper() for line in f if line.strip()]
    return []
if st.session_state.do_trimming and not st.session_state.trimming_done:
    if not st.session_state.uploaded_paths:
        st.error("No uploaded FASTQ paths available for trimming.")
    else:
        tmp_dir = st.session_state.temp_dir or tempfile.mkdtemp()
        st.session_state.temp_dir = tmp_dir
        uploaded_paths = st.session_state.uploaded_paths
        trimming_params_input = st.session_state.tool_params
        trimming_params = get_trimming_params(trimming_params_input)
        qc_params = trimming_params_input.get("Quality Control", {})
        user_adapter = qc_params.get("Input Adapter Sequence", "").upper().strip()
        user_adapters = [user_adapter] if user_adapter else []
        user_contaminant = qc_params.get("Contaminants/Primers Sequence", "").upper().strip()
        user_contaminants = [user_contaminant] if user_contaminant else []
        default_adapters = load_list_from_file("adapters.txt")
        default_contaminants = load_list_from_file("contaminants.txt")
        adapters = list(dict.fromkeys(default_adapters + user_adapters))
        contaminants = list(dict.fromkeys(default_contaminants + user_contaminants))
        trimmed_files = {}
        labels = ["SE"] if st.session_state.end_type == "Single End" else ["R1", "R2"]
        trim_progress_placeholder = st.empty()
        trim_status_placeholder = st.empty()
        trim_progress_bar = trim_progress_placeholder.progress(0)
        total_reads_est = 0
        for path in uploaded_paths:
            total_reads_est += sum(1 for _ in read_fastq(path))
        current_trim_step = [0]
        def update_trim_progress(increment=1, step_name="Trimming"):
            current_trim_step[0] += increment
            pct = int((current_trim_step[0] / max(1, total_reads_est)) * 100)
            pct = min(pct, 100)
            trim_progress_bar.progress(pct)
            trim_status_placeholder.markdown(
                f"<h3 style='color:white;'>⚙️ {step_name} ({pct}%)</h3>",
                unsafe_allow_html=True
            )
        for idx, file_path in enumerate(uploaded_paths):
            label = labels[idx] if idx < len(labels) else f"R{idx+1}"
            def trimming_progress_callback(pct, step_name="Trimming"):
                update_trim_progress(increment=1, step_name=f"{label} Trimming")
            out_name = os.path.join(tmp_dir, f"trimmed_{label}_{int(time.time())}.fastq")
            output_file, stats = trim_reads(
                file_path,
                out_name,
                adapters=adapters,
                contaminants=contaminants,
                user_params=trimming_params,
                progress_callback=trimming_progress_callback
            )
            trimmed_files[label] = (output_file, stats)
        trim_progress_bar.progress(100)
        trim_progress_bar.empty()
        trim_status_placeholder.markdown("<h3 style='color:white;'>✅ Trimming Completed</h3>", unsafe_allow_html=True)
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        trimming_entry = {"paths": {}, "stats": {}, "timestamp": timestamp}
        for label, (path, stats) in trimmed_files.items():
            st.session_state[f"trimmed_{label}_path"] = path
            st.session_state[f"trimmed_{label}_stats"] = stats
            trimming_entry["paths"][label] = path
            trimming_entry["stats"][label] = stats
            st.session_state.cached_trimming_stats[label] = stats
        st.session_state.trimming_history.append(trimming_entry)
        st.session_state.trimming_done = True
        st.session_state.fastqc_on_trim_done = False
# Always display trimming stats from cache
if st.session_state.trimming_done and st.session_state.cached_trimming_stats:
    st.markdown("<h3 style='color:white;'>✅ Trimming Statistics</h3>", unsafe_allow_html=True)
    if st.session_state.end_type == "Single End":
        label = "SE"
        if label in st.session_state.cached_trimming_stats:
            st.json({f"{label} - Trimming Statistics": st.session_state.cached_trimming_stats[label]})
    else:
        tab_r1, tab_r2 = st.tabs(["Read 1 (R1)", "Read 2 (R2)"])
        st.markdown(TABS_CSS, unsafe_allow_html=True)
        for tab, label in zip([tab_r1, tab_r2], ["R1", "R2"]):
            with tab:
                if label in st.session_state.cached_trimming_stats:
                    st.json({f"{label} - Trimming Statistics": st.session_state.cached_trimming_stats[label]})


# Backend Call Code for FASTQC
if st.session_state.trimming_done and not st.session_state.fastqc_on_trim_done and st.session_state.do_trimming:
    fastqc = FastQCBackend()
    tmp_dir = st.session_state.temp_dir
    trimmed_paths = []
    if st.session_state.end_type == "Single End":
        trimmed_paths = [st.session_state.get("trimmed_SE_path")]
    else:
        trimmed_paths = [st.session_state.get("trimmed_R1_path"), st.session_state.get("trimmed_R2_path")]
    trimmed_paths = [p for p in trimmed_paths if p is not None]
    if not trimmed_paths:
        st.error("Trimmed files not found for FastQC stage.")
    else:
        fastqc_progress_placeholder = st.empty()
        fastqc_status_placeholder = st.empty()
        fastqc_progress_bar = fastqc_progress_placeholder.progress(0)
        fastqc_files = len(trimmed_paths)
        current_fastqc_step = [0]
        def fastqc_progress_callback(pct, step_name="FastQC"):
            overall_pct = int((current_fastqc_step[0] + pct / 100) / max(1, fastqc_files) * 100)
            overall_pct = min(overall_pct, 100)
            fastqc_progress_bar.progress(overall_pct)
            fastqc_status_placeholder.markdown(
                f"<h3 style='color:white;'>⚙️ {step_name} ({overall_pct}%)</h3>",
                unsafe_allow_html=True
            )
        try:
            if st.session_state.end_type == "Single End":
                for label, path in zip(["SE"], trimmed_paths):
                    results, outdirs = fastqc.run_fastqc(
                        path, params=st.session_state.tool_params, end_type="Single End",
                        progress_callback=lambda pct, step_name=f"{label} FastQC": fastqc_progress_callback(pct, step_name)
                    )
                    current_fastqc_step[0] += 1
                    st.session_state.fastqc_on_trim_history.append({
                        "results": results,
                        "outdirs": outdirs,
                        "timestamp": time.strftime("%Y%m%d-%H%M%S")
                    })
            else:
                results, outdirs = fastqc.run_fastqc(
                    trimmed_paths, params=st.session_state.tool_params, end_type="Paired End",
                    progress_callback=lambda pct, step_name="FastQC": fastqc_progress_callback(pct, step_name)
                )
                current_fastqc_step[0] += 2
                st.session_state.fastqc_on_trim_history.append({
                    "results": results,
                    "outdirs": outdirs,
                    "timestamp": time.strftime("%Y%m%d-%H%M%S")
                })
        except Exception as e:
            st.error(f"FastQC on trimmed reads failed: {e}")
        fastqc_progress_bar.progress(100)
        fastqc_progress_bar.empty()
        fastqc_status_placeholder.empty()
        st.session_state.fastqc_on_trim_done = True
# Always display FastQC-on-trim results from history
if st.session_state.fastqc_on_trim_history:
    st.markdown("<h3 style='color:white;'>✅ FastQC Analysis on Trimmed Reads</h3>", unsafe_allow_html=True)
    for entry in st.session_state.fastqc_on_trim_history:
        results = entry["results"]
        if st.session_state.end_type == "Single End":
            _render_fastqc_single(results)
        else:
            r1_result = {k[len("R1_"):]: v for k, v in results.items() if k.startswith("R1_")}
            r2_result = {k[len("R2_"):]: v for k, v in results.items() if k.startswith("R2_")}
            tab_r1, tab_r2 = st.tabs(["Read 1 (R1)", "Read 2 (R2)"])
            st.markdown(TABS_CSS, unsafe_allow_html=True)
            for tab, result, label in zip([tab_r1, tab_r2], [r1_result, r2_result], ["R1", "R2"]):
                with tab:
                    _render_fastqc_single(result, label_prefix=f"{label} - ")

# Backend Call Code for BWA Alignment
from bwa_backend import SimpleBWAligner
import pandas as pd
should_run_bwa = False
if (st.session_state.get("do_trimming") is False and st.session_state.get("fastqc_done") and not st.session_state.get("bwa_done")) or \
   (st.session_state.get("do_trimming") is True and st.session_state.get("trimming_done") and not st.session_state.get("bwa_done")):
    should_run_bwa = True
if should_run_bwa:
    if not st.session_state.get("reference_path") or not st.session_state.get("uploaded_paths"):
        st.error("Please upload reference FASTA and FASTQ file(s).")
    else:
        temp_dir = st.session_state.get("temp_dir") or tempfile.mkdtemp()
        ref_path = st.session_state["reference_path"]
        fq_path_list = []
        if st.session_state.get("do_trimming") and st.session_state.get("trimming_done"):
            if st.session_state.get("end_type") == "Single End":
                trimmed = st.session_state.get("trimmed_SE_path")
                fq_path_list = [trimmed] if trimmed else [st.session_state["uploaded_paths"][0]]
            else:
                r1 = st.session_state.get("trimmed_R1_path")
                r2 = st.session_state.get("trimmed_R2_path")
                fq_path_list = [r1, r2] if r1 and r2 else st.session_state["uploaded_paths"][:2]
        else:
            fq_path_list = st.session_state["uploaded_paths"][:2] if len(st.session_state["uploaded_paths"]) >= 2 else [st.session_state["uploaded_paths"][0]]
        params = {"Reference Mapping": st.session_state.get("reference_mapping_params", st.session_state.default_params["Reference Mapping"])}
        bwa = SimpleBWAligner(ref_path, params)
        output_sam = os.path.join(temp_dir, "alignment.sam")
        progress_bar = st.progress(0)
        progress_text = st.empty()
        def progress_callback(percent, step_name):
            percent = max(0, min(100, int(percent)))
            progress_bar.progress(percent)
            progress_text.markdown(
                f"<h4 style='color:white; font-weight:bold;'>🔄 {step_name} - {percent}%</h4>",
                unsafe_allow_html=True
            )
        try:
            progress_callback(0, "Initializing alignment")
            if len(fq_path_list) == 1:
                output_sam = os.path.join(temp_dir, "output_single.sam")
                aligned_reads = bwa.align_reads_single(
                    fq_path_list[0], output_sam, progress_callback=progress_callback
                )
            else:
                output_sam = os.path.join(temp_dir, "output_paired.sam")
                bwa.align_reads_paired(fq_path_list[0], fq_path_list[1], output_sam, progress_callback=progress_callback)
                aligned_reads = []
            progress_callback(100, "Finalizing")
            progress_bar.empty()
            progress_text.empty()
            if os.path.exists(output_sam):
                try:
                    bam_path = output_sam.replace(".sam", ".bam")
                    bwa.convert_sam_to_bam(output_sam)
                except Exception as e:
                    st.error(f"BAM conversion error: {e}")
               # Build and CACHE SAM DataFrame
                sam_data = []
                with open(output_sam, "r") as sam:
                    for line in sam:
                        if line.startswith("@"):
                            continue
                        fields = line.strip().split("\t")
                        if len(fields) < 11:
                            continue
                        sam_data.append({
                            "QNAME": fields[0], "FLAG": fields[1], "RNAME": fields[2],
                            "POS": fields[3], "MAPQ": fields[4], "CIGAR": fields[5],
                            "RNEXT": fields[6], "PNEXT": fields[7], "TLEN": fields[8],
                            "SEQ": fields[9], "QUAL": fields[10]
                        })
                st.session_state.cached_sam_df = pd.DataFrame(sam_data) if sam_data else pd.DataFrame()
            st.session_state["bwa_done"] = True
        except Exception as e:
            progress_bar.empty()
            progress_text.empty()
            st.error(f"Alignment failed: {e}")
# Always display BWA results from cache
if st.session_state.get("bwa_done"):
    st.markdown("<h3 style='color:white;'>✅ Alignment Completed</h3>", unsafe_allow_html=True)
    if st.session_state.cached_sam_df is not None:
        if not st.session_state.cached_sam_df.empty:
            st.dataframe(st.session_state.cached_sam_df, width=1500, height=400)
        else:
            st.warning("No alignments to display in SAM file.")

# Backend Call Code for Freebayes
from freebayes import FreeBayesParallelCaller as SimpleFreeBayesCaller
import glob
if st.session_state.get("bwa_done") and not st.session_state.get("freebayes_done"):
    temp_dir = st.session_state.get("temp_dir") or tempfile.mkdtemp()
    sample_name = st.session_state.get("sample_name", "Sample1")
    ref_path = st.session_state.get("reference_path")
    bam_candidates = glob.glob(os.path.join(temp_dir, "*.bam"))
    if not bam_candidates:
        st.error(f"❌ No BAM file found in {temp_dir}. Please ensure BWA step completed successfully.")
    else:
        bam_file = bam_candidates[0]
        output_vcf = os.path.join(temp_dir, f"{sample_name}_variants.vcf")
        if not os.path.exists(ref_path):
            st.error("❌ Reference FASTA not found. Please upload a valid reference genome.")
        else:
            progress_bar = st.progress(0, text="Initializing FreeBayes...")
            progress_text = st.empty()
            freebayes = SimpleFreeBayesCaller(
                reference_fasta=ref_path,
                bam_path=bam_file,
                params=st.session_state.tool_params,
                progress_callback=progress_callback
            )
            try:
                progress_callback(0, "Counting Haplotypes")
                freebayes.run(output_vcf=output_vcf)
                progress_callback(100, "Finalizing and writing VCF")
                progress_bar.empty()
                progress_text.empty()
                if os.path.exists(output_vcf):
                    st.session_state["freebayes_done"] = True
                    # Build and CACHE VCF DataFrame
                    vcf_records = []
                    with open(output_vcf, "r") as vcf:
                        for line in vcf:
                            if line.startswith("#"):
                                continue
                            fields = line.strip().split("\t")
                            if len(fields) < 8:
                                continue
                            vcf_records.append({
                                "CHROM": fields[0], "POS": fields[1], "ID": fields[2],
                                "REF": fields[3], "ALT": fields[4],
                                "FILTER": fields[6], "INFO": fields[7]
                            })
                    st.session_state.cached_vcf_df = pd.DataFrame(vcf_records) if vcf_records else pd.DataFrame()
                else:
                    st.error("FreeBayes did not produce a VCF file.")
            except Exception as e:
                progress_bar.empty()
                progress_text.empty()
                st.error(f"FreeBayes execution failed: {e}")
# Always display FreeBayes results from cache
if st.session_state.get("freebayes_done"):
    st.markdown("<h3 style='color:white;'>✅ Variant Analysis Completed</h3>", unsafe_allow_html=True)
    if st.session_state.cached_vcf_df is not None:
        if not st.session_state.cached_vcf_df.empty:
            st.dataframe(st.session_state.cached_vcf_df, width=1500, height=400)
        else:
            st.warning("No variants detected in VCF file.")


# Backend Call Code for Variant Annotation
from annotation import annotate
if st.session_state.get("freebayes_done") and not st.session_state.get("vep_done"):
    temp_dir = st.session_state.get("temp_dir")
    vcf_candidates = glob.glob(os.path.join(temp_dir, "*_variants.vcf"))
    if not vcf_candidates:
        st.error("❌ No VCF file found from FreeBayes.")
        st.stop()
    input_vcf = vcf_candidates[0]
    if reference_file is None:
        st.warning("⚠️ Please upload Reference FASTA and GTF/GFF file.")
        st.stop()
    fasta_file_ann = None
    gff_file_ann = None
    for f in reference_file:
        fname = f.name.lower()
        if fname.endswith((".fa", ".fasta")):
            fasta_file_ann = f
        elif fname.endswith((".gff", ".gff3", ".gtf")):
            gff_file_ann = f
    if fasta_file_ann is None:
        st.error("❌ Reference FASTA file not found")
        st.stop()
    if gff_file_ann is None:
        st.error("⚠️ No GFF/GTF found → Ensembl REST will be used")
        st.stop()
    fasta_path = os.path.join(temp_dir, "reference.fa")
    with open(fasta_path, "wb") as f:
        f.write(fasta_file_ann.getbuffer())
    gtf_path = os.path.join(temp_dir, "annotation.gtf")
    with open(gtf_path, "wb") as f:
        f.write(gff_file_ann.getbuffer())
    fasta_files = glob.glob(os.path.join(temp_dir, "*.fa")) + glob.glob(os.path.join(temp_dir, "*.fasta"))
    gtf_files = glob.glob(os.path.join(temp_dir, "*.gtf")) + glob.glob(os.path.join(temp_dir, "*.gff")) + glob.glob(os.path.join(temp_dir, "*.gff3"))
    if not fasta_files or not gtf_files:
        st.error("❌ Both Reference FASTA and GTF/GFF must be uploaded.")
        st.stop()
    fasta_path = fasta_files[0]
    gtf_path = gtf_files[0]
    progress_bar = st.progress(0)
    progress_text = st.empty()
    with open(input_vcf) as f:
        total_variants = sum(1 for l in f if not l.startswith("#"))
    def update_progress(current):
        percent = min(int((current / total_variants) * 100), 100)
        progress_bar.progress(percent)
        progress_text.markdown(
            f"<h4 style='color:white;'>🔄 Annotating variants {current}/{total_variants} ({percent}%)</h4>",
            unsafe_allow_html=True
        )
        time.sleep(0.002)
    try:
        output_vcf_path = os.path.join(temp_dir, "annotated.vcf")
        df = annotate(vcf_path=input_vcf, fasta=fasta_path, gtf=gtf_path, output_vcf=output_vcf_path)
        for i in range(len(df)):
            update_progress(i + 1)
        progress_bar.progress(100)
        progress_bar.empty()
        progress_text.empty()
        st.session_state["vep_done"] = True
        st.session_state.cached_annotated_df = df
        st.session_state["annotated_df"] = df
        st.session_state["vep_full_df"] = df
    except Exception as e:
        progress_bar.empty()
        progress_text.empty()
        st.error(f"❌ Variant Annotation failed: {e}")
# Always display annotation results from cache
if st.session_state.get("vep_done"):
    st.markdown("<h3 style='color:white;'>✅ Variant Annotation Completed</h3>", unsafe_allow_html=True)
    if st.session_state.cached_annotated_df is not None and not st.session_state.cached_annotated_df.empty:
        st.dataframe(st.session_state.cached_annotated_df, width=1600, height=450)
    else:
        st.warning("⚠️ No variants were annotated")


# Backend Call Code for Variant Filtering
from filtering import VariantFilter
if st.session_state.get("vep_done") and not st.session_state.get("filter_done"):
    temp_dir = st.session_state.get("temp_dir")
    annotated_vcf_path = os.path.join(temp_dir, "annotated.vcf")
    filter_params = st.session_state.default_params["Filter"]
    progress_bar = st.progress(0)
    progress_text = st.empty()
    def update_progress(current, total):
        percent = min(int((current / total) * 100), 100)
        progress_bar.progress(percent)
        progress_text.markdown(
            f"<h4 style='color:white;'>🔄 Filtering variants {current}/{total} ({percent}%)</h4>",
            unsafe_allow_html=True
        )
        time.sleep(0.001)
    try:
        vf = VariantFilter(
            annotated_vcf_path=annotated_vcf_path,
            filter_params=filter_params,
            progress_callback=update_progress
        )
        vf.load_vcf()
        filtered_df = vf.pass_filters()
        if isinstance(filtered_df, pd.DataFrame):
            if "Impact" not in filtered_df.columns:
                filtered_df["Impact"] = "LOW"
            else:
                filtered_df["Impact"] = filtered_df["Impact"].astype(str).str.upper()
        else:
            filtered_df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "Impact", "DP", "AF"])
        # CACHE filtered DataFrame
        st.session_state.cached_filtered_df = filtered_df.copy()
        st.session_state["filtered_df"] = filtered_df.copy()
        st.session_state["filter_done"] = True
        progress_bar.progress(100)
        progress_bar.empty()
        progress_text.empty()
    except Exception as e:
        st.error(f"❌ Variant Filtering failed: {e}")
# Always display filtering results from cache
if st.session_state.get("filter_done"):
    st.markdown("<h3 style='color:white;'>✅ Variant Filtering Completed</h3>", unsafe_allow_html=True)
    if st.session_state.cached_filtered_df is not None and not st.session_state.cached_filtered_df.empty:
        st.dataframe(st.session_state.cached_filtered_df, width=1600, height=450)
    else:
        st.warning("⚠️ No variants passed filtering")


# Backend Call Code for Genome Visualization
import plotly.graph_objects as go
import re
if st.session_state.get("filter_done"):
    st.markdown("<h3 style='color:white;'>✅ Genome View</h3>", unsafe_allow_html=True)
    temp_dir = st.session_state.get("temp_dir")
    fasta_path_list = glob.glob(os.path.join(temp_dir, "*.fa")) + \
                      glob.glob(os.path.join(temp_dir, "*.fasta"))
    sam_path_list = glob.glob(os.path.join(temp_dir, "*.sam"))
    fasta_path = fasta_path_list[0]
    sam_path = sam_path_list[0]
    # LOAD REFERENCE
    ref_seq = ""
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                ref_seq += line.strip().upper()
    genome_length = len(ref_seq)
    st.markdown("""
    <style>
    label[data-testid="stWidgetLabel"]{
        color: white !important;
        font-weight: bold;
    }
    </style>
    """, unsafe_allow_html=True)
    region_start, region_end = st.slider(
        "Select Genomic Region",
        1,
        genome_length,
        (1, min(500, genome_length))
    )
    region_size = region_end - region_start
    show_letters = region_size <= 10000000000000
    base_colors = {
        "A": "#2ecc71", "T": "#e74c3c", "U": "#e74c3c",
        "C": "#3498db", "G": "#f39c12", "N": "#bdc3c7"
    }
    ref_x, ref_y, ref_color, ref_text = [], [], [], []
    for pos in range(region_start, region_end):
        base = ref_seq[pos - 1]
        ref_x.append(pos + 0.5)
        ref_y.append(0)
        ref_color.append(base_colors.get(base, "#bdc3c7"))
        ref_text.append(base if show_letters else "")
    reads = []
    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            pos = int(fields[3])
            cigar = fields[5]
            seq = fields[9].upper()
            if cigar == "*":
                continue
            unmapped = (flag & 4) != 0
            reverse = (flag & 16) != 0
            if unmapped:
                continue
            ref_length = 0
            cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
            for length, op in cigar_tuples:
                length = int(length)
                if op in ["M", "D", "N", "=", "X"]:
                    ref_length += length
            ref_end = pos + ref_length
            if ref_end < region_start or pos > region_end:
                continue
            reads.append({"start": pos, "end": ref_end, "seq": seq, "cigar_tuples": cigar_tuples, "reverse": reverse})
    if not reads:
        st.warning("⚠ No reads found in this region")
    else:
        reads.sort(key=lambda x: x["start"])
        rows = []
        layout = []
        for read in reads:
            placed = False
            for i in range(len(rows)):
                if read["start"] > rows[i]:
                    rows[i] = read["end"]
                    layout.append((i + 1, read))
                    placed = True
                    break
            if not placed:
                rows.append(read["end"])
                layout.append((len(rows), read))
        read_x, read_y, read_color, read_text = [], [], [], []
        arrow_x, arrow_y, arrow_text = [], [], []
        for row_index, read in layout:
            ref_pointer = read["start"]
            seq_pointer = 0
            for length, op in read["cigar_tuples"]:
                length = int(length)
                if op in ["M", "=", "X"]:
                    for i in range(length):
                        genomic_pos = ref_pointer + i
                        if genomic_pos < region_start or genomic_pos >= region_end:
                            continue
                        base = read["seq"][seq_pointer + i]
                        ref_base = ref_seq[genomic_pos - 1]
                        color = base_colors.get(base, "#bdc3c7")
                        if base != ref_base:
                            color = "#4b44ad"
                        read_x.append(genomic_pos + 0.5)
                        read_y.append(row_index)
                        read_color.append(color)
                        read_text.append(base if show_letters else "")
                    ref_pointer += length
                    seq_pointer += length
                elif op == "I":
                    for i in range(length):
                        genomic_pos = ref_pointer
                        if genomic_pos < region_start or genomic_pos >= region_end:
                            continue
                        read_x.append(genomic_pos + 0.5)
                        read_y.append(row_index)
                        read_color.append("#b65991")
                        read_text.append(read["seq"][seq_pointer + i] if show_letters else "")
                    seq_pointer += length
                elif op in ["D", "N"]:
                    for i in range(length):
                        genomic_pos = ref_pointer + i
                        if genomic_pos < region_start or genomic_pos >= region_end:
                            continue
                        read_x.append(genomic_pos + 0.5)
                        read_y.append(row_index)
                        read_color.append("#000000")
                        read_text.append("")
                    ref_pointer += length
                elif op == "S":
                    seq_pointer += length
            arrow_x.append(read["start"] + 0.5)
            arrow_y.append(row_index)
            arrow_text.append("←" if read["reverse"] else "→")
        fig = go.Figure()
        square_size = 50
        row_spacing = 1
        ref_y = [0 for _ in ref_y]
        read_y = [row * row_spacing for row in read_y]
        arrow_y = [row * row_spacing for row in arrow_y]
        fig.add_trace(go.Scatter(
            x=ref_x, y=ref_y,
            mode="markers+text" if show_letters else "markers",
            marker=dict(symbol="square", size=square_size, color=ref_color, line=dict(width=0)),
            text=ref_text, textposition="middle center",
            textfont=dict(size=12, color="white"),
            hoverinfo="skip", name="Reference", showlegend=True
        ))
        fig.add_trace(go.Scatter(
            x=read_x, y=read_y,
            mode="markers+text" if show_letters else "markers",
            marker=dict(symbol="square", size=square_size, color=read_color, line=dict(width=0)),
            text=read_text, textposition="middle center",
            textfont=dict(size=11, color="white"),
            hoverinfo="skip", name="Aligned Reads", showlegend=True
        ))
        fig.add_trace(go.Scatter(
            x=arrow_x, y=arrow_y, mode="text",
            text=arrow_text, textposition="middle left",
            textfont=dict(size=14, color="black"),
            hoverinfo="skip", name="Strand Direction", showlegend=True
        ))
        fig.update_layout(
            height=350 + len(rows) * 30,
            plot_bgcolor="white", paper_bgcolor="white",
            margin=dict(l=10, r=10, t=40, b=20),
            xaxis=dict(range=[min(ref_x) - 0.5, max(ref_x) + 0.5], showgrid=False, zeroline=False, constrain="domain"),
            yaxis=dict(visible=False, range=[-0.5, len(rows) + 0.5], scaleanchor="x", scaleratio=0.5),
            dragmode="pan",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
        )
        st.plotly_chart(fig, use_container_width=True)


    # Backend Call Code for Report Generation
    st.markdown("---")
    st.markdown("<h3 style='color:white;'>📄 Report Generation</h3>", unsafe_allow_html=True)
    vcf_files = glob.glob(os.path.join(temp_dir, "*.vcf"))
    if not vcf_files:
        st.warning("⚠ No VCF file found for report generation")
    else:
        vcf_file = vcf_files[0]
        st.success(f"✅ VCF detected: {os.path.basename(vcf_file)}")
    if st.button("🚀 Generate Report", use_container_width=True):
        if not vcf_files:
            st.error("❌ Cannot generate report without VCF")
            st.stop()
        with st.spinner("Generating report..."):
            try:
                from report import run_report
                output_pdf = run_report(temp_dir, vcf_file)
                st.success("✅ Report Generated Successfully!")
                with open(output_pdf, "rb") as f:
                    st.download_button(
                        label="📥 Download Report",
                        data=f,
                        file_name="VariantXplorer_Report.pdf",
                        mime="application/pdf",
                        use_container_width=True
                    )
                pdf_path = run_report(temp_dir, vcf_file)
                st.success(f"Report saved to: {pdf_path}")            
            except Exception as e:
                st.error(f"❌ Report generation failed: {str(e)}")