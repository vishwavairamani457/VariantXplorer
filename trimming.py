# Trimming Backend Code
import gzip
import streamlit as st
# Parameters
def get_trimming_params(user_input):
    """
    Convert user input from frontend/tool into trimming parameters.
    """
    params = {
        "window_size": 5,
        "window_qual": 10,
        "leading_quality": 3,
        "trailing_quality": 3,
        "read_crop_length": 0,
        "read_head_crop_length": 0,
        "min_read_length": 10,
        "min_gc_content": 40,
        "max_gc_content": 60,
        "ns_max_p": 1,
        "ns_max_n": 3,
        "lc_threshold": 10
    }
    trim_params = user_input.get("Trimming", {})
    params["window_size"] = int(trim_params.get("Window Size", params["window_size"]))
    params["window_qual"] = int(trim_params.get("Window Qual", params["window_qual"]))
    params["leading_quality"] = int(trim_params.get("Leading Quality", params["leading_quality"]))
    params["trailing_quality"] = int(trim_params.get("Trailing Quality", params["trailing_quality"]))
    params["read_crop_length"] = int(trim_params.get("Read Crop Length", params["read_crop_length"]))
    params["read_head_crop_length"] = int(trim_params.get("Read Head Crop Length", params["read_head_crop_length"]))
    params["min_read_length"] = int(trim_params.get("Min Read Length", params["min_read_length"]))
    params["min_gc_content"] = int(trim_params.get("Min GC Content", params["min_gc_content"]))
    params["max_gc_content"] = int(trim_params.get("Max GC Content", params["max_gc_content"]))
    params["ns_max_p"] = float(trim_params.get("ns_max_p  (Percentage of Ambiguous Nucleotide)", params["ns_max_p"]))
    params["ns_max_n"] = int(trim_params.get("ns_max_n  (No.of Ambiguous Nucleotide)", params["ns_max_n"]))
    params["lc_threshold"] = int(trim_params.get("lc_threshold (Low Complexity Sequence)", params["lc_threshold"]))
    return params

# FASTQ Reader
def read_fastq(file_path):
    open_func = gzip.open if file_path.endswith(".gz") else open
    with open_func(file_path, "rt", newline=None) as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline().strip()
            qual = fh.readline().strip()
            if header and seq and plus and qual:
                yield header, seq, plus, qual

# Trimming Helpers
def trim_hits(seq, qual, hits):
    trimmed_flag = False
    count = 0
    for hit in hits:
        while hit in seq:
            idx = seq.find(hit)
            seq = seq[:idx] + seq[idx+len(hit):]
            qual = qual[:idx] + qual[idx+len(hit):]
            trimmed_flag = True
            count += 1
    return seq, qual, trimmed_flag, count

def trim_leading(seq, qual, threshold):
    cut = 0
    for q in qual:
        if ord(q)-33 < threshold:
            cut += 1
        else:
            break
    return seq[cut:], qual[cut:], cut > 0

def trim_trailing(seq, qual, threshold):
    cut = 0
    for q in reversed(qual):
        if ord(q)-33 < threshold:
            cut += 1
        else:
            break
    if cut > 0:
        return seq[:-cut], qual[:-cut], True
    return seq, qual, False

def trim_sliding_window(seq, qual, window_size, window_qual):
    cut_pos = len(seq)
    for i in range(len(seq)-window_size+1):
        window = [ord(q)-33 for q in qual[i:i+window_size]]
        if sum(window)/window_size < window_qual:
            cut_pos = i
            break
    if cut_pos < len(seq):
        return seq[:cut_pos], qual[:cut_pos], True
    return seq, qual, False

def check_gc(seq, min_gc, max_gc):
    if not seq:
        return False
    gc = (seq.count("G")+seq.count("C"))/len(seq)*100
    return min_gc <= gc <= max_gc

# Main Trimming Function (with progress callback)
def trim_reads(input_fastq, output_prefix, adapters=None, contaminants=None, user_params=None, progress_callback=None, step_name="Trimming"):
    adapters = adapters or []
    contaminants = contaminants or []
    params = {
        "window_size": 5,
        "window_qual": 10,
        "leading_quality": 3,
        "trailing_quality": 3,
        "read_crop_length": 0,
        "read_head_crop_length": 0,
        "min_read_length": 10,
        "min_gc_content": 40,
        "max_gc_content": 60,
        "ns_max_p": 1,
        "ns_max_n": 3,
        "lc_threshold": 10
    }
    if user_params:
        params.update(user_params)
    output_fastq = f"{output_prefix}_trimmed.fastq"
    stats = {
        "total_reads": 0,
        "trimmed_reads": 0,
        "dropped_reads": 0,
        "adapter_hits_trimmed": 0,
        "contaminant_hits_trimmed": 0,
        "total_bp_before": 0,
        "total_bp_after": 0,
    }
    # Count total reads first for % completion
    total_reads_est = sum(1 for _ in read_fastq(input_fastq))
    processed_reads = 0
    with open(output_fastq, "w") as out:
        for header, seq, plus, qual in read_fastq(input_fastq):
            processed_reads += 1
            stats["total_reads"] += 1
            stats["total_bp_before"] += len(seq)
            trimmed_flag = False
            # Adapters
            seq, qual, t1, adapter_count = trim_hits(seq, qual, adapters)
            stats["adapter_hits_trimmed"] += adapter_count
            if t1: trimmed_flag = True
            # Contaminants
            seq, qual, t2, contam_count = trim_hits(seq, qual, contaminants)
            stats["contaminant_hits_trimmed"] += contam_count
            if t2: trimmed_flag = True
            # Quality trimming
            seq, qual, t3 = trim_leading(seq, qual, params["leading_quality"]); trimmed_flag |= t3
            seq, qual, t4 = trim_trailing(seq, qual, params["trailing_quality"]); trimmed_flag |= t4
            seq, qual, t5 = trim_sliding_window(seq, qual, params["window_size"], params["window_qual"]); trimmed_flag |= t5
            if params["read_head_crop_length"] > 0:
                seq, qual = seq[params["read_head_crop_length"]:], qual[params["read_head_crop_length"]:]
                trimmed_flag = True
            if params["read_crop_length"] > 0 and len(seq) > params["read_crop_length"]:
                seq, qual = seq[:params["read_crop_length"]], qual[:params["read_crop_length"]]
                trimmed_flag = True
            # Filters
            if not check_gc(seq, params["min_gc_content"], params["max_gc_content"]):
                stats["dropped_reads"] += 1
                continue
            if len(seq) < params["min_read_length"]:
                stats["dropped_reads"] += 1
                continue
            if trimmed_flag:
                stats["trimmed_reads"] += 1
            stats["total_bp_after"] += len(seq)
            out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
            # Update progress if callback provided
            if progress_callback and total_reads_est > 0:
                percent = min(100, int((processed_reads / total_reads_est) * 100))  # clamp to 100
                progress_callback(percent, f"{step_name}... {percent}%")
    try:
        import streamlit as st
        st.session_state["output_fastq_1"] = output_fastq
        st.session_state["trimming_stats"] = stats   # 🔥 IMPORTANT
    except:
        pass    
    # Compute mean read length and GC content
    mean_read_length = round(
        stats["total_bp_after"] / (stats["total_reads"] - stats["dropped_reads"]), 2
    ) if (stats["total_reads"] - stats["dropped_reads"]) > 0 else 0
    gc_count = 0
    for _, seq, _, _ in read_fastq(output_fastq):
        gc_count += seq.count("G") + seq.count("C")
    mean_gc_content = round((gc_count / stats["total_bp_after"])*100, 2) if stats["total_bp_after"] > 0 else 0
    stats["mean_read_length"] = mean_read_length
    stats["mean_gc_content"] = mean_gc_content
    # Final 100% callback
    if progress_callback:
        progress_callback(100, f"{step_name} Done")
    import streamlit as st
    st.session_state["trimming_stats"] = stats
    return output_fastq, stats