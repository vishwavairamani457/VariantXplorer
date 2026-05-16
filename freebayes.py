# Freebayes Backend Code
import math, sys
import tempfile
import shutil
import threading
import struct, zlib
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import as_completed
from itertools import combinations_with_replacement
import numpy as np
import re, itertools
from collections import defaultdict

# Utilities
PHRED_CAP = 60
def phred_to_prob(q):
    return 10 ** (-q / 10.0)
def safe_char(seq, i, default="-"):
    if seq is None:
        return default
    if i < 0 or i >= len(seq):
        return default
    return seq[i]
def clamp_phred(q):
    return max(0, min(int(q), PHRED_CAP))
def phred_from_prob(p):
    p = max(p, 1e-6)
    q = -10.0 * math.log10(p)
    return int(min(999, round(q)))
def multinomial_log_likelihood(counts, probs, tiny=1e-300):
    ll = 0.0
    total = sum(counts)
    if total == 0: return 0.0
    for c, p in zip(counts, probs):
        p = max(p, tiny)
        if c > 0:
            ll += c * math.log(p)
    return ll
def phred_scale_from_probs(probs):
    tiny = 1e-300
    clean_probs = []
    for p in probs:
        if isinstance(p, str):
            p = p.replace(",", ".")  
            try:
                p = float(p)
            except:
                p = tiny
        else:
            p = float(p) if p is not None else tiny
        clean_probs.append(max(p, tiny)) 
    logs = [-10.0 * math.log10(p) for p in clean_probs]
    mn = min(logs) if logs else 0
    pls = [int(round(x - mn)) for x in logs]
    return pls
def genotype_order(n_alleles, ploidy):
    return list(combinations_with_replacement(range(n_alleles), ploidy))
def cigar_to_ops(cigar):
    if not cigar or cigar == "*":
        return []
    matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if not matches:
        return []
    return [(op, int(length)) for length, op in matches]
def parse_cigar(cigar):
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
def left_trim_ref_alt(ref, alt, pos, ref_sequence):
    r = ref
    a = alt
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        pos += 1
    changed = True
    while changed:
        changed = False
        if len(r) > len(a):  # deletion (ref longer)
            if pos > 1:
                prev_base = ref_sequence[pos - 2]  # pos is 1-based
                if r[-1] == prev_base:
                    r = prev_base + r[:-1]
                    a = prev_base + a
                    pos -= 1
                    changed = True
        elif len(a) > len(r):  # insertion (alt longer)
            if pos > 1:
                prev_base = ref_sequence[pos - 2]
                if a[-1] == prev_base:
                    r = prev_base + r
                    a = prev_base + a[:-1]
                    pos -= 1
                    changed = True
    return r, a, pos
# Smith-Waterman implementation (returns alignment and score)
def smith_waterman(ref, read, match=2, mismatch=-3, gap_open=-6, gap_extend=-1):
    if not ref or not read:
        return {
            "aligned_ref": "",
            "aligned_read": "",
            "score": 0,
            "ref_start": 1,
            "ref_end": 1
        }
    n = len(ref)
    m = len(read)
    H = np.zeros((n+1, m+1), dtype=int)
    E = np.zeros((n+1, m+1), dtype=int)
    F = np.zeros((n+1, m+1), dtype=int)
    TB = np.zeros((n+1, m+1), dtype=np.int8)
    best_score = 0
    best_i = best_j = 0
    for i in range(1, n+1):
        ri = ref[i-1]
        for j in range(1, m+1):
            s = match if ri == read[j-1] else mismatch
            diag = H[i-1, j-1] + s
            E[i, j] = max(H[i-1, j] + gap_open, E[i-1, j] + gap_extend)
            F[i, j] = max(H[i, j-1] + gap_open, F[i, j-1] + gap_extend)
            val = max(0, diag, E[i, j], F[i, j])
            H[i, j] = val
            if val == 0:
                TB[i, j] = 0
            elif val == diag:
                TB[i, j] = 1
            elif val == E[i, j]:
                TB[i, j] = 2
            else:
                TB[i, j] = 3
            if val > best_score:
                best_score = val
                best_i, best_j = i, j
    # backtrack
    i, j = best_i, best_j
    aligned_ref = []
    aligned_read = []
    while i > 0 and j > 0 and TB[i, j] != 0:
        t = TB[i, j]
        if t == 1:
            aligned_ref.append(ref[i-1]); aligned_read.append(read[j-1]); i -= 1; j -= 1
        elif t == 2:
            aligned_ref.append(ref[i-1]); aligned_read.append('-'); i -= 1
        else:
            aligned_ref.append('-'); aligned_read.append(read[j-1]); j -= 1
    aligned_ref = "".join(reversed(aligned_ref))
    aligned_read = "".join(reversed(aligned_read))
    # SAFETY: ensure same length
    min_len = min(len(aligned_ref), len(aligned_read))
    aligned_ref = aligned_ref[:min_len]
    aligned_read = aligned_read[:min_len]
    start_on_ref = i + 1
    end_on_ref = best_i
    return {"aligned_ref": aligned_ref, "aligned_read": aligned_read, "score": int(best_score), "ref_start": start_on_ref, "ref_end": end_on_ref}
def detect_complex_variant(read_seq, ref_seq):
    mismatches = []
    L = min(len(read_seq), len(ref_seq))
    for i in range(L):
        if read_seq[i] != ref_seq[i]:
            mismatches.append(i)
    if len(mismatches) < 2:
        return None
    start = mismatches[0]
    end = mismatches[-1] + 1
    ref = ref_seq[start:end]
    alt = read_seq[start:end]
    if len(ref) == 1 and len(alt) == 1:
        return None
    return ref, alt

class BGZFReader:
    def __init__(self, path):
        self.fh = open(path, "rb")
        self.buffer = b""
        self.offset = 0

    def _read_block(self):
        header = self.fh.read(18)
        if len(header) == 0:
            return b""
        _ = struct.unpack("<BBBBIBBH", header)
        xlen = struct.unpack("<BBBBIBBH", header)[-1]
        extra = self.fh.read(xlen)
        bsize = struct.unpack("<H", extra[4:6])[0] + 1
        block = header + extra + self.fh.read(bsize - (12 + xlen))
        data = block[:-8]  
        return zlib.decompress(data, wbits=-15)

    def read(self, n):
        while len(self.buffer) < n:
            block = self._read_block()
            if not block:
                break
            self.buffer += block
        out = self.buffer[:n]
        self.buffer = self.buffer[n:]
        return out

    def readline(self):
        while True:
            pos = self.buffer.find(b"\n")
            if pos >= 0:
                out = self.buffer[:pos+1]
                self.buffer = self.buffer[pos+1:]
                return out
            block = self._read_block()
            if not block:
                out = self.buffer
                self.buffer = b""
                return out
            self.buffer += block

    def close(self):
        self.fh.close()

def decode_cigar(raw, n_ops):
    ops = "MIDNSHP=XB"
    out = []
    for i in range(n_ops):
        val = raw[i]
        length = val >> 4
        op = ops[val & 0xF]
        out.append(f"{length}{op}")
    return "".join(out)

def decode_seq(nibbles):
    table = "=ACMGRSVTWYHKDBN"
    return "".join(table[(b >> 4)] + table[(b & 0xF)] for b in nibbles)[:len(nibbles)*2]

def decode_bam_record(data, refs):
    (
        ref_id, pos,
        read_name_bin,
        cigar_flag,
        l_seq,
        next_ref_id, next_pos, tlen
    ) = struct.unpack("<iiIIiiii", data[:32])
    # unpack bitfields
    l_read_name =  read_name_bin & 0xFF
    mapq        = (read_name_bin >> 8) & 0xFF
    bin_        = (read_name_bin >> 16) & 0xFFFF
    n_cigar     =  cigar_flag & 0xFFFF
    flag        = (cigar_flag >> 16) & 0xFFFF
    qname = data[32:32 + l_read_name - 1].decode()
    offset = 32 + l_read_name
    cigar_raw = struct.unpack("<" + "I" * n_cigar, data[offset:offset + 4 * n_cigar])
    offset += 4 * n_cigar
    cigar = decode_cigar(cigar_raw, n_cigar)
    seq_bytes = data[offset:offset + (l_seq + 1) // 2]
    offset += (l_seq + 1) // 2
    seq = decode_seq(seq_bytes)[:l_seq]
    qual_raw = data[offset:offset + l_seq]
    offset += l_seq
    qual = [q for q in qual_raw]
    tags = {}
    rname = refs[ref_id] if ref_id >= 0 else "*"
    next_rname = refs[next_ref_id] if next_ref_id >= 0 else "*"
    return {
        "qname": qname,
        "flag": flag,
        "rname": rname,
        "pos": pos + 1,
        "mapq": mapq,
        "cigar": cigar,
        "seq": seq,
        "qual": qual,
        "tags": tags,
        "next_ref_id": next_ref_id,
        "next_rname": next_rname,
        "next_pos": next_pos + 1 if next_pos >= 0 else 0,
        "tlen": tlen,
        "bin": bin_,
    }

def read_uncompressed_bam(bam_path):
    fh = open(bam_path, "rb")
    # Header magic
    fh.read(4)  # BAM\1
    # Read text header
    l_text = struct.unpack("<i", fh.read(4))[0]
    _ = fh.read(l_text)  # header string
    # Read reference sequences
    n_ref = struct.unpack("<i", fh.read(4))[0]
    refs = []
    ref_lengths = []
    for _ in range(n_ref):
        l_name = struct.unpack("<i", fh.read(4))[0]
        name = fh.read(l_name).rstrip(b"\0").decode()
        length = struct.unpack("<i", fh.read(4))[0]
        refs.append(name)
        ref_lengths.append(length)  # store length
    # Iterate alignments
    while True:
        raw = fh.read(4)
        if len(raw) < 4:
            break
        block_size = struct.unpack("<i", raw)[0]
        data = fh.read(block_size)
        if len(data) < block_size:
            break
        yield decode_bam_record(data, refs)
    fh.close()

def stream_alignments(sam_path=None, bam_path=None):
    if sam_path:
        with open(sam_path, "rt") as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 11:
                    continue
                qname = parts[0]; flag = int(parts[1]); rname = parts[2]
                pos = int(parts[3]) if parts[3].isdigit() else 0
                mapq = int(parts[4]) if parts[4].isdigit() else 0
                cigar = parts[5]
                seq = parts[9]
                qual_str = parts[10]
                qual = [ord(c) - 33 for c in qual_str] if qual_str != "*" else []
                tags = {}
                for t in parts[11:]:
                    if ":" in t:
                        k = t.split(":", 1)[0]
                        tags[k] = t
                yield {
                    "qname": qname, "flag": flag, "rname": rname,
                    "pos": pos, "mapq": mapq, "cigar": cigar,
                    "seq": seq, "qual": qual, "tags": tags, 
                    "next_ref_id": None, "next_pos": None,
                    "tlen": None, "bin": None
                }
        return
    if bam_path is None:
        raise ValueError("Require sam_path or bam_path")
    with open(bam_path, "rb") as raw:
        first4 = raw.read(4)
    if first4 == b"BAM\x01":
        for aln in read_uncompressed_bam(bam_path):
            yield aln
        return
    if first4 != b"\x1f\x8b\x08\x04":
        raise RuntimeError(
            f"Invalid BAM: expected BGZF or BAM\\1, got {first4!r}"
        )
    bgzf = BGZFReader(bam_path)
    bgzf.buffer = first4 + bgzf.buffer
    l_text = struct.unpack("<i", bgzf.read(4))[0]
    _ = bgzf.read(l_text).decode()
    # read reference sequences
    n_ref = struct.unpack("<i", bgzf.read(4))[0]
    refs = []
    ref_lengths = []
    for _ in range(n_ref):
        l_name = struct.unpack("<i", bgzf.read(4))[0]
        name = bgzf.read(l_name).rstrip(b"\x00").decode()
        length = struct.unpack("<i", bgzf.read(4))[0]
        refs.append(name)
        ref_lengths.append(length)
    # iterate alignment records
    while True:
        raw = bgzf.read(4)
        if len(raw) < 4:
            break   
        block_size = struct.unpack("<i", raw)[0]
        data = bgzf.read(block_size)
        if len(data) < block_size:
            break
        (
            ref_id, pos,
            read_name_bin,
            cigar_flag,
            l_seq,
            next_ref_id, next_pos, tlen
        ) = struct.unpack("<iiIIiiii", data[:32])
        # unpack bitfields
        l_read_name =  read_name_bin & 0xFF
        mapq        = (read_name_bin >> 8) & 0xFF
        bin_        = (read_name_bin >> 16) & 0xFFFF
        n_cigar     =  cigar_flag & 0xFFFF
        flag        = (cigar_flag >> 16) & 0xFFFF
        qname = data[32:32 + l_read_name - 1].decode()
        offset = 32 + l_read_name
        # CIGAR
        cigar_raw = struct.unpack("<" + "I" * n_cigar,data[offset:offset + 4 * n_cigar])
        offset += 4 * n_cigar
        cigar = decode_cigar(cigar_raw, n_cigar)
        # SEQ
        seq_bytes = data[offset:offset + ((l_seq + 1) // 2)]
        offset += (l_seq + 1) // 2
        seq = decode_seq(seq_bytes)[:l_seq]
        # QUAL
        qual_raw = data[offset:offset + l_seq]
        offset += l_seq
        qual = list(qual_raw)
        # TAGS
        tags = {}
        rest = data[offset:]
        i = 0
        while i < len(rest):
            tag = rest[i:i+2].decode()
            type_ = chr(rest[i+2])
            i += 3
            if type_ == "Z":
                end = rest.index(b"\x00", i)
                val = rest[i:end].decode()
                i = end + 1
            elif type_ == "i":
                val = struct.unpack("<i", rest[i:i+4])[0]
                i += 4
            else:
                # unsupported type — break to avoid errors
                break
            tags[tag] = f"{tag}:{type_}:{val}"
        rname = refs[ref_id] if ref_id >= 0 else "*"
        yield {
            "qname": qname,
            "flag": flag,
            "rname": rname,
            "pos": pos + 1,
            "mapq": mapq,
            "cigar": cigar,
            "seq": seq,
            "qual": qual,
            "tags": tags,
            "next_ref_id": next_ref_id,
            "next_pos": next_pos + 1,
            "tlen": tlen,
            "bin": bin_
        }
    bgzf.close()

class Read:
    def __init__(self, chrom, pos, base, base_qual, mapq):
        self.chrom = chrom
        self.pos = pos       
        self.base = base
        self.base_qual = base_qual
        self.mapq = mapq

def safe_div(a, b):
    return 0.0 if b == 0 else a / b

# Main FreeBayes-style caller class
class FreeBayesParallelCaller:

    def _normalize_dict(self, obj):
        if isinstance(obj, dict):
            return obj
        if hasattr(obj, "__dict__"):
            return dict(obj.__dict__)
        return {}

    def extract_variants_from_read(self, ref_seq, read_seq, cigar, ref_start):
        variants = []
        ref_pos = ref_start - 1
        read_pos = 0
        for length, op in parse_cigar(cigar):
            if op in ("M", "=", "X"):
                for _ in range(length):
                    if ref_pos >= len(ref_seq) or read_pos >= len(read_seq):
                        break
                    ref_base = ref_seq[ref_pos]
                    read_base = read_seq[read_pos]
                    ref = ref_base
                    alt = read_base
                    if ref != alt:
                        variants.append({
                            "pos": ref_pos + 1,
                            "ref": ref,
                            "alt": alt,
                            "type": "snp"
                        })

                    ref_pos += 1
                    read_pos += 1
            elif op == "I":
                ins = read_seq[read_pos:read_pos + length]
                ref = "-"
                alt = ins
                if ref != alt:
                    variants.append({
                        "pos": ref_pos,
                        "ref": ref,
                        "alt": alt,
                        "type": "ins"
                    })
                read_pos += length
            elif op == "D":
                dele = ref_seq[ref_pos:ref_pos + length]
                ref = dele
                alt = "-"
                if ref != alt:
                    variants.append({
                        "pos": ref_pos + 1,
                        "ref": ref,
                        "alt": alt,
                        "type": "del"
                    })
                ref_pos += length
            elif op in ("S", "H"):
                read_pos += length
            else:
                ref_pos += length
        return variants

    # Extract only Variant Calling block
    def _extract_variant_calling_params(self, params):
        if not isinstance(params, dict):
            return {}
        if "Variant Calling" in params and isinstance(params["Variant Calling"], dict):
            return params["Variant Calling"]
        return params

    # Load defaults + apply incoming overrides
    def get_reference_params(self, default_params=None):
        defaults = {
            "ploidy": 1,
            "min_mapping_quality": 0,
            "min_base_quality": 0,
            "reference-mapping-quality": 100,
            "reference-base-quality": 60,
            "min_coverage": 1,
            "min_alt_count": 1,
            "min_alt_fraction": 0.01,
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
        }
        # print("\n[DEBUG] Loading FreeBayes default parameters")
        params = defaults.copy()
        incoming = self._normalize_dict(default_params) or {}
        nested = incoming.get("Variant Calling", {}) if isinstance(incoming, dict) else {}
        # print("[DEBUG] Incoming override (flat):", incoming)
        # print("[DEBUG] Incoming override (nested):", nested)
        # Apply flat keys
        for k, v in incoming.items():
            if k == "Variant Calling":
                continue
            if k in params:
                params[k] = v
        # Apply nested keys LAST
        if isinstance(nested, dict):
            for k, v in nested.items():
                if k in params:
                    params[k] = v
        # Cast types according to defaults
        for key, default_val in defaults.items():
            if key in params:
                try:
                    params[key] = type(default_val)(params[key])
                except:
                    params[key] = default_val
        # print("[DEBUG] Final params after applying defaults + incoming:", params)
        self.params = params
        return params

    # Final user param override step
    def update_params(self, user_params=None):
        if not user_params:
            # print("[DEBUG] update_params(): No user overrides")
            return self.params
        user_params = self._normalize_dict(user_params)
        nested = user_params.get("Variant Calling", {})
        # print("\n[DEBUG] update_params() START")
        # print("[DEBUG] Incoming flat:", user_params)
        # print("[DEBUG] Incoming nested:", nested)
        # Apply flat
        for key, value in user_params.items():
            if key == "Variant Calling":
                continue
            if key in self.params:
                # print(f"[DEBUG] Overriding {key} -> {value}")
                self.params[key] = value
        # Apply nested
        if isinstance(nested, dict):
            for key, value in nested.items():
                if key in self.params:
                    # print(f"[DEBUG] Overriding nested {key} -> {value}")
                    self.params[key] = value
        # print("[DEBUG] update_params() FINAL:", self.params)
        return self.params
    # Constructor
    def __init__(self,
                 reference_fasta,
                 sam_path=None,
                 bam_path=None,
                 params=None,
                 threads=4,
                 window_size=1000000,
                 sample_name="SAMPLE",
                 keep_temp=False, progress_callback=None):
        # Store progress callback
        self.progress_callback = progress_callback 
        # Progress counters
        self.total_reads = 0      
        self.processed_reads = 0 
        # Ref load
        self.reference = Path(reference_fasta)
        if not self.reference.exists():
            raise FileNotFoundError("Reference FASTA not found: " + str(self.reference))
        self.ref_dict = self._load_reference()
        self.contigs = list(self.ref_dict.keys())
        # IO
        self.sam_path = sam_path
        self.bam_path = bam_path
        if not (sam_path or bam_path):
            raise ValueError("Provide --sam or --bam")
        self.threads = int(threads)
        self.window_size = int(window_size)
        self.sample_name = sample_name
        self.keep_temp = keep_temp
        self.tempdir = Path(tempfile.mkdtemp(prefix="pyfb_"))
        # Extract only the Variant Calling block
        normalized_defaults = self._extract_variant_calling_params(params)
        # print("\n[DEBUG] __init__ incoming params:", params)
        # print("[DEBUG] Extracted VC block:", normalized_defaults)
        # Load defaults + shallow incoming override
        self.get_reference_params(normalized_defaults)
        # Final override pass
        self.update_params(normalized_defaults)
        # Read index
        self.reads_by_contig = defaultdict(list)
        self._index_reads()
        # self.total_reads = sum(len(v) for v in self.reads_by_contig.values())
        # print(f"[DEBUG] Total reads indexed: {self.total_reads}")
        self.completed_haplotype_regions = 0
        self.total_haplotype_regions = 0
        self.lock = threading.Lock()

    def _count_reads(self):
        count = 0
        for _ in stream_alignments(bam_path=self.bam_path):
            count += 1
        return max(1, count)

    def _load_reference(self):
        ref = {}
        with open(self.reference, "r") as fh:
            name = None; parts = []
            for line in fh:
                line = line.rstrip("\n")
                if not line: continue
                if line.startswith(">"):
                    if name:
                        ref[name] = "".join(parts).upper()
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
            if name:
                ref[name] = "".join(parts).upper()
        return ref

    def _index_reads(self):
        # Reset counters
        self.reads_by_contig = defaultdict(list)
        self.total_reads = 0
        self.processed_reads = 0
        if self.params.get("debug"):
            print("[DEBUG] Counting total reads...")
        for aln in stream_alignments(sam_path=self.sam_path, bam_path=self.bam_path):
            rname = aln["rname"]
            pos = aln["pos"]
            if rname == "*" or pos == 0:
                continue
            self.total_reads += 1
        if self.params.get("debug"):
            print(f"[DEBUG] Total reads found: {self.total_reads}")
        if self.params.get("debug"):
            print("[DEBUG] Indexing reads into memory. This may be slow and memory-heavy for large BAMs.")
        # reset stream
        count_stream = stream_alignments(sam_path=self.sam_path, bam_path=self.bam_path)
        for aln in count_stream:
            rname = aln["rname"]
            pos = aln["pos"]
            if rname == "*" or pos == 0:
                continue
            # store fields (original logic untouched)
            self.reads_by_contig[rname].append({
                "qname": aln["qname"], "flag": aln["flag"], "pos": aln["pos"], "mapq": aln["mapq"],
                "cigar": aln["cigar"], "seq": aln["seq"], "qual": aln["qual"], "tags": aln["tags"]
            })
            # PROGRESS UPDATE
            self.processed_reads += 1
            if self.progress_callback and self.total_reads > 0:
                pct = int((self.processed_reads / self.total_reads) * 100)
                msg = f"Indexing reads... {self.processed_reads}/{self.total_reads}"
                self.progress_callback(pct, msg)
        # print per-contig debug summary
        if self.params.get("debug"):
            for c in self.contigs:
                print(f"[DEBUG] {c}: {len(self.reads_by_contig.get(c, []))} reads indexed")

    def make_windows(self):
        windows = []
        for chrom, seq in self.ref_dict.items():
            L = len(seq)
            for start in range(1, L+1, self.window_size):
                end = min(L, start + self.window_size - 1)
                windows.append((chrom, start, end))
        return windows

    def _reads_in_window(self, chrom, start1, end1):
        res = []
        for r in self.reads_by_contig.get(chrom, []):
            aln_start = r["pos"]
            end = aln_start
            for op, length in cigar_to_ops(r["cigar"]):
                if op in ("M", "D", "N", "=", "X"):
                    end += length
            if end - 1 < start1 or aln_start > end1:
                continue
            res.append(r)
        return res

    def build_pileup(self, reads, reference):
        pileup = defaultdict(list)
        ref_len = len(reference)
        for read in reads:
            seq = read.get("seq", "")
            qual = read.get("qual", [])
            cigar = read.get("cigar", "")
            pos = read.get("pos", 1)
            if not seq or not cigar:
                continue
            read_i = 0
            ref_i = pos - 1
            for length, op in parse_cigar(cigar):
                if op in ("M", "=", "X"):
                    for _ in range(length):
                        if read_i >= len(seq):
                            break
                        if 0 <= ref_i < ref_len:
                            base = seq[read_i]
                            pileup[ref_i + 1].append(base)
                        read_i += 1
                        ref_i += 1
                elif op == "I":
                    ins = seq[read_i:read_i + length]
                    anchor = max(ref_i - 1, 0)
                    pileup[anchor].append("+" + ins)
                    read_i += length
                elif op == "D":
                    del_seq = reference[ref_i:ref_i + length]
                    anchor = max(ref_i - 1, 0)
                    pileup[anchor].append("-" + del_seq)
                    ref_i += length
                elif op == "S":
                    read_i += length
                elif op in ("H", "N"):
                    ref_i += length
                if "I" in cigar or "D" in cigar:
                    print("INDEL READ:", read.pos, cigar)            
        return pileup
    
    def _project_read(self, read, reference):
        try:
            seq = read.get("seq", "") or ""
            qual = read.get("qual", []) or []
            mapq = read.get("mapq", 0)
            seq = str(seq)
            if seq == "*" or len(seq) == 0:
                return [], []
            # convert ASCII QUAL if needed
            if isinstance(qual, str):
                qual = [ord(q) - 33 for q in qual]
            # ensure QUAL length matches SEQ
            if len(qual) < len(seq):
                qual = qual + [0] * (len(seq) - len(qual))
            strand = "+" if not read.get("is_reverse", False) else "-"
            cigar = read.get("cigar")
            if not cigar:
                return [], []
            ops = cigar_to_ops(cigar)
            ref_pos = read.get("pos", 0)
            qptr = 0
            projected = []
            insertions = []
            for op, length in ops:
                if op in ("M", "=", "X"):
                    for _ in range(length):
                        if qptr >= len(seq):
                            break
                        b = seq[qptr]
                        # SAFE QUAL ACCESS
                        if qptr < len(qual):
                            bq = qual[qptr]
                        else:
                            bq = 0
                        projected.append((ref_pos, b, bq, mapq, strand))
                        ref_pos += 1
                        qptr += 1
                elif op == "I":
                    end = min(qptr + length, len(seq))
                    ins = seq[qptr:end]
                    if len(ins) == 0:
                        qptr = end
                        continue
                    ins_qs = qual[qptr:end] if qptr < len(qual) else []
                    if len(ins_qs) == 0:
                        ins_qs = [0] * len(ins)
                    anchor = max(1, ref_pos - 1)
                    avgq = int(sum(ins_qs) / len(ins_qs))
                    insertions.append((anchor, ins, avgq, mapq))
                    qptr = end
                elif op == "D":
                    del_seq = reference[ref_pos-1:ref_pos-1+length]
                    token = "-" + del_seq
                    projected.append((ref_pos, token, 0, mapq, strand))
                    ref_pos += length
                elif op == "N":
                    ref_pos += length
                elif op == "S":
                    qptr = min(qptr + length, len(seq))
                elif op == "H":
                    pass
                else:
                    qptr = min(qptr + length, len(seq))
                    ref_pos += length
            return projected, insertions
        except Exception as e:
            print("\n===== PROJECT READ ERROR =====")
            print("Read object:", read)
            try:
                print("SEQ length:", len(read.get("seq", "")))
                print("QUAL length:", len(read.get("qual", [])))
                print("CIGAR:", read.get("cigar"))
                print("POS:", read.get("pos"))
            except:
                pass
            raise e           

    def _discover_sites_in_window(self, chrom, start1, end1, reads):
        from collections import defaultdict, Counter
        counts = defaultdict(Counter)        
        base_records = defaultdict(list)
        indel_positions = set()
        min_mapq = int(self.params.get("min_mapping_quality"))
        min_bq = int(self.params.get("min_base_quality"))
        min_cov = int(self.params.get("min_coverage"))
        min_alt_count = int(self.params.get("min_alt_count"))
        ref_seq = self.ref_dict[chrom]
        for read in reads:
            if read["mapq"] < min_mapq:
                continue
            proj, insertions = self._project_read(read, ref_seq)
            strand = read.get("strand", "+")
            for rp, tok, bq, mq, strand in proj:
                if rp < start1 or rp >= end1:
                    continue
                if rp - 1 < 0 or rp - 1 >= len(ref_seq):
                    continue
                if isinstance(tok, str) and len(tok) == 1:
                    if bq < min_bq:
                        continue
                counts[rp][tok] += 1
                base_records[rp].append({
                    "tok": tok,
                    "bq": bq,
                    "mq": mq,
                    "strand": strand
                })
                if isinstance(tok, str) and tok.startswith("-"):
                    indel_positions.add(rp)
            # Insertions
            for anchor, ins, avgq, mq in insertions:
                if anchor < start1 or anchor >= end1:
                    continue
                token = "+" + ins
                counts[anchor][token] += 1
                base_records[anchor].append({
                    "tok": token,
                    "bq": avgq,
                    "mq": mq,
                    "strand": strand
                })
                indel_positions.add(anchor)
        # Variant filtering
        filtered_counts = {}
        filtered_records = {}
        for pos, counter in counts.items():
            if pos < start1 or pos >= end1:
                continue
            total = sum(counter.values())
            if total < min_cov:
                continue
            if pos <= 0 or pos > len(ref_seq):
                continue
            ref_base = ref_seq[pos - 1]
            alts = defaultdict(int)
            for tok, c in counter.items():
                # INSERTION
                if isinstance(tok, str) and tok.startswith("+"):
                    ins = tok[1:]
                    ref_base = ref_seq[pos - 1]
                    ref = ref_base
                    alt = ref_base + ins
                    alts[(ref, alt)] += c   
                    continue
                # DELETION
                if isinstance(tok, str) and tok.startswith("-"):
                    del_seq = tok[1:]
                    ref_base = ref_seq[pos - 1]
                    ref = ref_base + del_seq
                    alt = ref_base
                    alts[(ref, alt)] += c   
                    continue
                # SNP
                if isinstance(tok, str) and len(tok) == 1:
                    ref = ref_seq[pos - 1]
                    alt = tok
                    if alt == ref:
                        continue
                    alts[(ref, alt)] += c
                    continue
                # complex
                if isinstance(tok, str) and len(tok) > 1 and not tok.startswith(("+", "-")):
                    ref = ref_seq[pos - 1]
                    alt = tok
                    alts[(ref, alt)] += c
            if not alts:
                continue
            valid_alt = False
            for (ref, alt), c in alts.items():
                is_indel = len(ref) != 1 or len(alt) != 1
                if is_indel:
                    if c >= 1:
                        valid_alt = True
                        break
                else:
                    alt_fraction = c / total
                    if c >= min_alt_count and alt_fraction >= self.params.get("min_alt_fraction"):
                        valid_alt = True
                        break
            if not valid_alt:
                continue
            filtered_counts[pos] = counter
            filtered_records[pos] = base_records[pos]
        return filtered_counts, filtered_records, indel_positions
    
    def _allele_prior(self, depth):
        theta = float(self.params.get("theta"))
        pvar = float(self.params.get("pvar"))
        rdf  = float(self.params.get("read_dependence_factor"))
        # Effective mutation rate
        mu = theta * pvar
        # Read dependence scaling of evidence
        scale = max(0.1, min(1.0, rdf ** max(1, depth)))
        return mu * scale

    def _assemble_haplotypes(
            self, 
            chrom, 
            reads, 
            core_positions, 
            window_start=1, 
            window_end=None, 
            region_start=None, 
            region_end=None):
        max_hap_len = min(200, int(self.params.get("max_haplotype_length")))
        flanking = min(80, max_hap_len // 2)
        if window_end is None:
            window_end = len(self.ref_dict[chrom])
        L = region_start if region_start is not None else max(window_start, min(core_positions) - flanking)
        R = region_end   if region_end   is not None else min(window_end, max(core_positions) + flanking)
        region_ref = self.ref_dict[chrom][L-1:R]
        hap_counter = Counter()
        hap_info = {}
        min_sup = int(self.params.get("min_haplotype_support"))
        for read in reads:
            proj, insertions = self._project_read(read, region_ref)
            seq_parts = []
            mq_sum = 0
            bq_sum = 0
            support = 0
            has_alt = False
            obs_map = {}
            for rp, base, bq, mq, strand in proj:
                if L <= rp <= R:
                    if rp not in obs_map or bq > obs_map[rp][1]:
                        obs_map[rp] = (base, bq, mq)
            ins_map = {}
            for anchor, ins_seq, avgq, mq in insertions:
                if L <= anchor <= R:
                    ins_map[anchor] = (ins_seq, avgq, mq)
            hap_seq_builder = []
            pos = L
            while pos <= R:
                if pos in obs_map:
                    base, bq, mq = obs_map[pos]
                    hap_seq_builder.append(base)
                    seq_parts.append(base)
                    bq_sum += bq
                    mq_sum += mq
                    support += 1
                    idx = pos - L
                    if idx < 0 or idx >= len(region_ref):
                        ref_base = "N"
                    else:
                        ref_base = safe_char(region_ref, idx, "N")
                    if base != ref_base:
                        has_alt = True
                else:
                    idx = pos - L
                    if 0 <= idx < len(region_ref):
                        hap_seq_builder.append(region_ref[idx])
                    else:
                        hap_seq_builder.append("N")                    
                if pos in ins_map:
                    ins_seq, avgq, mq = ins_map[pos]
                    # append actual inserted sequence
                    hap_seq_builder.append(ins_seq)
                    seq_parts.append(ins_seq)
                    bq_sum += avgq
                    mq_sum += mq
                    support += 1
                    has_alt = True
                pos += 1
            hap_string = "".join(hap_seq_builder)
            # Normalize haplotype length to region length
            target_len = len(region_ref)
            if len(hap_string) > target_len:
                hap_string = hap_string[:target_len]
            elif len(hap_string) < target_len:
                hap_string = hap_string + "-" * (target_len - len(hap_string))
            is_alt = has_alt
            if support == 0 and is_alt and min_sup == 1:
                support = 1
            if support < min_sup:
                continue
            key = tuple(seq_parts)
            if key not in hap_info:
                hap_info[key] = {
                    "seq": hap_string,
                    "support": 0,
                    "bq": 0,
                    "mq": 0,
                    "seq_parts": []
                }
            hap_info[key]["support"] += 1
            hap_info[key]["bq"] += bq_sum
            hap_info[key]["mq"] += mq_sum
            hap_info[key]["seq_parts"].extend(seq_parts)
        max_haps = int(self.params.get("max_haplotypes"))
        if len(hap_info) > max_haps * 4:
            hap_info = dict(sorted(hap_info.items(), key=lambda x: -x[1]["support"])[:max_haps*4])
        hap_items = sorted(hap_info.values(), key=lambda x: -x["support"])
        haplotypes = []
        for info in hap_items[:max_haps]:
            seq = info.get("seq", "")
            if not seq:
                continue
            haplotypes.append({
                "seq": info["seq"],
                "start": L,
                "support": info["support"],
                "bq_sum": info["bq"],
                "mq_sum": info["mq"],
                "seq_parts": info["seq_parts"]
            })
        ref_seq = region_ref
        ref_support = sum(1 for h in haplotypes if h["seq"] == ref_seq)
        ref_hap = {
            "seq": ref_seq,
            "start": L,
            "support": len(reads),
            "bq_sum": 0,
            "mq_sum": 0,
            "seq_parts": list(ref_seq)
        }
        if not any(h["seq"] == ref_seq for h in haplotypes):
            haplotypes.insert(0, ref_hap)
        max_haplotypes = 128
        if len(haplotypes) > max_haplotypes:
            haplotypes = haplotypes[:max_haplotypes]
        clean_haps = []
        for h in haplotypes:
            seq = h.get("seq", "")
            if not seq:
                continue
            if set(seq) == {"-"}:
                continue
            gap_fraction = seq.count("-") / len(seq)
            if gap_fraction > 0.8:
                continue
            if len(seq) > len(region_ref) + 20:
                continue
            clean_haps.append(h)
        haplotypes = clean_haps
        if not haplotypes:
            haplotypes = [ref_hap]
        max_len = max(len(h["seq"]) for h in haplotypes)
        for h in haplotypes:
            if len(h["seq"]) < max_len:
                pad = "-" * (max_len - len(h["seq"]))
                h["seq"] = h["seq"] + pad      
        # print("Haplotypes generated:", haplotypes)
        return haplotypes, L, R
  
    def _haplotype_to_alleles(self, ref_hap_seq, hap_seq, ref_start, ref_sequence, ref_name):
        alleles = []
        i = 0
        while i < len(hap_seq):
            ref_base = ref_hap_seq[i] if i < len(ref_hap_seq) else "-"
            alt_base = hap_seq[i]
            pos = ref_start + i
            # SNP
            if alt_base != ref_base and alt_base != "-" and ref_base != "-":
                alleles.append({
                    "pos": pos,
                    "ref": ref_base,
                    "alt": alt_base
                })
                i += 1
                continue
            # DELETION
            if alt_base == "-" and ref_base != "-":
                j = i
                while j < len(hap_seq) and hap_seq[j] == "-":
                    j += 1
                del_len = j - i
                ref = ref_sequence[pos-1 : pos-1 + del_len + 1]
                alt = ref[0]
                alleles.append({
                    "pos": pos,
                    "ref": ref,
                    "alt": alt
                })
                i = j
                continue
            # INSERTION
            if ref_base == "-" and alt_base != "-" and i > 0:
                j = i
                ins = ""
                while j < len(hap_seq) and ref_hap_seq[j] == "-":
                    ins += hap_seq[j]
                    j += 1
                anchor = pos - 1
                ref = ref_sequence[anchor]
                alt = ref + ins
                alleles.append({
                    "pos": anchor + 1,
                    "ref": ref,
                    "alt": alt
                })
                i = j
                continue
            i += 1
            print("REF HAP:", ref_hap_seq)
            print("ALT HAP:", hap_seq)
        return {
            "alleles": alleles,
            "ref_alignment": ref_hap_seq,
            "alt_alignment": hap_seq,
            "ref_name": ref_name
        }

    def ref_base_at_position(self, pos1_based, contig=None):
        if contig is None:
            raise ValueError("Contig must be specified.")
        if contig not in self.ref_dict:
            return None  
        seq = self.ref_dict[contig]
        pos0 = pos1_based - 1  
        if pos0 < 0 or pos0 >= len(seq):
            return None  
        return seq[pos0]
    
    def _compute_genotype_likelihoods_from_haplotypes(
        self, alleles_refalt, reads, chrom, region_start, region_end
    ):
        region_ref = self.ref_dict[chrom][region_start-1:region_end]
        n_alleles = 1 + len(alleles_refalt) 
        # Build haplotype sequences for each allele
        allele_hap_seqs = [region_ref]       
        for pos, ref, alt in alleles_refalt:
            rel = pos - region_start           
            prefix = region_ref[:rel]            
            # For insertions and deletions, need to handle length changes
            ref_len = len(ref)
            alt_len = len(alt)      
            if alt_len > ref_len:
                suffix_start = rel + ref_len
                suffix = region_ref[suffix_start:] if suffix_start < len(region_ref) else ""
                newseq = prefix + alt + suffix
            elif ref_len > alt_len:
                suffix_start = rel + ref_len
                suffix = region_ref[suffix_start:] if suffix_start < len(region_ref) else ""
                newseq = prefix + alt + suffix
            else:
                suffix_start = rel + ref_len
                suffix = region_ref[suffix_start:] if suffix_start < len(region_ref) else ""
                newseq = prefix + alt + suffix         
            allele_hap_seqs.append(newseq)        
        # Normalize haplotype lengths (pad shorter ones with reference bases)
        max_len = max(len(h) for h in allele_hap_seqs)
        normalized = []
        for i, h in enumerate(allele_hap_seqs):
            if len(h) < max_len:
                # Pad with the actual reference sequence, not dashes
                pad_len = max_len - len(h)
                if i == 0:
                    # Reference haplotype - extend with reference
                    pad_start = len(h)
                    if region_start - 1 + pad_start < len(self.ref_dict[chrom]):
                        pad = self.ref_dict[chrom][region_start - 1 + pad_start: region_start - 1 + pad_start + pad_len]
                        if len(pad) < pad_len:
                            pad = pad + "N" * (pad_len - len(pad))
                    else:
                        pad = "N" * pad_len
                else:
                    pad = "N" * pad_len
                h = h + pad
            normalized.append(h)      
        allele_hap_seqs = normalized        
        min_bq = int(self.params.get("min_base_quality"))
        read_probs = []      
        # Compute read-to-haplotype alignment probabilities
        for read in reads:
            proj, insertions = self._project_read(read, region_ref)          
            read_seq = []
            read_quals = []
            read_mapq = read.get("mapq", 0) or 0            
            # Build read sequence from projection
            for rp, tok, bq, mq, strand in proj:
                if rp < region_start or rp > region_end:
                    continue               
                if isinstance(tok, str) and len(tok) == 1:
                    if bq >= min_bq:
                        read_seq.append(tok)
                        read_quals.append(bq)
                elif isinstance(tok, str) and tok.startswith("+"):
                    # Insertion in read
                    ins_bases = tok[1:]
                    for c in ins_bases:
                        if bq >= min_bq:
                            read_seq.append(c)
                            read_quals.append(bq)
                elif isinstance(tok, str) and tok.startswith("-"):
                    # Deletion in read - mark with gap indicator
                    del_len = len(tok) - 1
                    for _ in range(del_len):
                        read_seq.append("-")  
                        read_quals.append(bq)
            # Add separate insertions
            for anchor, ins_seq, avgq, ins_mq in insertions:
                if region_start <= anchor <= region_end:
                    if avgq >= min_bq:
                        read_seq.extend(list(ins_seq))
                        read_quals.extend([avgq] * len(ins_seq))
                    read_mapq = max(read_mapq, ins_mq)
            read_seq_str = "".join(read_seq)           
            # Handle empty reads
            if not read_seq_str or read_seq_str == "*":
                read_probs.append([1.0 / n_alleles] * n_alleles)
                continue            
            # Compute likelihood against each haplotype
            log_probs = []            
            for aseq in allele_hap_seqs:
                logL = 0.0                
                # Use Smith-Waterman for better alignment if sequences differ significantly
                if len(read_seq_str) != len(aseq) or read_seq_str != aseq:
                    # Align read to haplotype
                    sw_result = smith_waterman(aseq, read_seq_str)
                    aligned_ref = sw_result.get("aligned_ref", aseq)
                    aligned_read = sw_result.get("aligned_read", read_seq_str)                    
                    # Compute likelihood from alignment
                    L = min(len(aligned_ref), len(aligned_read), len(read_quals))
                    for i in range(L):
                        a = aligned_ref[i] if i < len(aligned_ref) else "-"
                        b = aligned_read[i] if i < len(aligned_read) else "-"
                        q = read_quals[i] if i < len(read_quals) else 20                      
                        error = 10 ** (-q / 10.0)                        
                        if a == b and a != "-":
                            p = 1.0 - error
                        elif a == "-" or b == "-":
                            p = error  # Gap penalty
                        else:
                            p = error / 3.0  # Mismatch                        
                        logL += math.log(max(p, 1e-12))
                else:
                    # Direct comparison for identical length sequences
                    L = min(len(read_seq_str), len(read_quals), len(aseq))
                    for i in range(L):
                        a = aseq[i]
                        b = read_seq_str[i]
                        q = read_quals[i] if i < len(read_quals) else 20                       
                        error = 10 ** (-q / 10.0)                        
                        if a == b:
                            p = 1.0 - error
                        else:
                            p = error / 3.0                        
                        logL += math.log(max(p, 1e-12))                
                log_probs.append(logL)           
            # Normalize to probabilities
            max_log = max(log_probs) if log_probs else 0
            probs = [math.exp(lp - max_log) for lp in log_probs]
            total = sum(probs)           
            if total > 0:
                probs = [p / total for p in probs]
            else:
                probs = [1.0 / n_alleles] * n_alleles           
            read_probs.append(probs)      
        # Compute genotype likelihoods
        ploidy = int(self.params.get("ploidy", 2))
        gts = genotype_order(n_alleles, ploidy)        
        gl_map = {}
        for gt in gts:
            log_total = 0.0
            for probs in read_probs:
                # Average likelihood across alleles in genotype
                p_obs = sum(probs[a] for a in gt) / len(gt)
                p_obs = max(p_obs, 1e-300)
                log_total += math.log(p_obs)
            gl_map[gt] = log_total        
        # Normalize genotype likelihoods
        max_gl = max(gl_map.values()) if gl_map else 0
        gl_norm = {}
        total = 0.0        
        for gt, gl in gl_map.items():
            val = math.exp(gl - max_gl)
            gl_norm[gt] = val
            total += val        
        if total > 0:
            for gt in gl_norm:
                gl_norm[gt] /= total        
        gl_map = gl_norm        
        # Compute allele-level MQ/BQ statistics
        allele_mq_sum = [0] * n_alleles
        allele_bq_sum = [0] * n_alleles        
        for read, probs in zip(reads, read_probs):
            best = int(np.argmax(probs))
            quals = read.get("qual", [])
            mean_bq = int(np.mean([q for q in quals if q >= min_bq])) if quals else 0
            allele_mq_sum[best] += read.get("mapq", 0)            
            proj, insertions = self._project_read(read, region_ref)
            for _, _, _, ins_mq in insertions:
                allele_mq_sum[best] += ins_mq            
            allele_bq_sum[best] += mean_bq        
        return gl_map, allele_mq_sum, allele_bq_sum

    def compute_allele_depths(self, ref_base, alt_alleles, pileup_entries):
        ref_count = 0
        alt_counts = {alt: 0 for alt in alt_alleles}
        alt_forward = {alt: 0 for alt in alt_alleles}
        alt_reverse = {alt: 0 for alt in alt_alleles}
        min_bq = self.params.get("min_base_quality")
        min_mq = self.params.get("min_mapping_quality")
        for rec in pileup_entries:
            tok = rec["tok"]
            bq = rec["bq"]
            mq = rec["mq"]
            strand = rec.get("strand", "+")
            # Convert quality if needed
            if isinstance(bq, str):
                bq = ord(bq) - 33
            if isinstance(mq, str):
                mq = ord(mq) - 33
            if bq < min_bq or mq < min_mq:
                continue
            # Case 1: Simple base (SNP or reference)
            if len(tok) == 1:
                if tok == ref_base:
                    ref_count += 1
                elif tok in alt_counts:
                    alt_counts[tok] += 1
                    if strand == "+":
                        alt_forward[tok] += 1
                    else:
                        alt_reverse[tok] += 1
                continue
            # Case 2: Insertion token (e.g., "+AG" means AG inserted after current position)
            if tok.startswith("+"):
                inserted_seq = tok[1:]  # The inserted bases
                # For insertion allele, alt should be ref_base + inserted_seq
                expected_alt = ref_base + inserted_seq
                if expected_alt in alt_counts:
                    alt_counts[expected_alt] += 1
                    if strand == "+":
                        alt_forward[expected_alt] += 1
                    else:
                        alt_reverse[expected_alt] += 1
                continue
            # Case 3: Deletion token (e.g., "-AG" means AG deleted starting at next position)
            if tok.startswith("-"):
                deleted_seq = tok[1:]  # The deleted bases
                # For deletion allele, ref should be ref_base + deleted_seq, alt = ref_base
                expected_ref = ref_base + deleted_seq
                for alt_allele in alt_alleles:
                    # Check if this is a matching deletion
                    if len(alt_allele) < len(expected_ref):
                        if expected_ref.startswith(alt_allele):
                            alt_counts[alt_allele] += 1
                            if strand == "+":
                                alt_forward[alt_allele] += 1
                            else:
                                alt_reverse[alt_allele] += 1
                continue
            # Case 4: Complex token - try direct match
            if tok in alt_counts:
                alt_counts[tok] += 1
                if strand == "+":
                    alt_forward[tok] += 1
                else:
                    alt_reverse[tok] += 1
        return ref_count, alt_counts, alt_forward, alt_reverse
    
    def compute_indel_depths(self, pos, alt_alleles, insertions_by_pos, deletions_by_pos):
        alt_counts = {alt: 0 for alt in alt_alleles}
        # Insertions
        if pos in insertions_by_pos:
            for ins, ins_q, mapq in insertions_by_pos[pos]:
                tag = f"+{ins}"
                if tag in alt_counts:
                    if ins_q >= self.params["min_base_quality"] and mapq >= self.params["min_mapping_quality"]:
                        alt_counts[tag] += 1
        # Deletions
        if pos in deletions_by_pos:
            for del_seq, del_q, mapq in deletions_by_pos[pos]:
                tag = f"-{del_seq}"
                if tag in alt_counts:
                    if del_q >= self.params["min_base_quality"] and mapq >= self.params["min_mapping_quality"]:
                        alt_counts[tag] += 1
        return alt_counts

    def _phred_to_error_prob(self, q):
        return 10 ** (-q / 10.0)

    def read_likelihood(self, read_base, allele, base_quality, mapq):
        # Base error probability
        e_base = 10 ** (-base_quality / 10.0)
        # Mapping error probability
        e_map = 10 ** (-mapq / 10.0)
        # Combine errors
        e = min(0.5, e_base + e_map - (e_base * e_map))
        if read_base == allele:
            return max(1.0 - e, 1e-6)
        else:
            return max(e / 3.0, 1e-6)

    def assign_read(self, read, haplotypes):
        read_base = read["base"]
        bq = int(read["bq"])
        mq = int(read["mq"])
        best_allele = None
        best_lk = -1.0
        for allele, hap_base in haplotypes.items():
            lk = self.read_likelihood(
                read_base,
                hap_base,
                bq,
                mq
            )
            if lk > best_lk:
                best_lk = lk
                best_allele = allele
        return best_allele, best_lk
        
    def genotype_likelihood_from_reads(self, reads, ref, alt, gt):
        ll = 0.0
        for base, bq, mq in reads:
            if isinstance(bq, str):
                bq = ord(bq) - 33
            if isinstance(mq, str):
                mq = ord(mq) - 33
            bq = max(1, int(bq))
            mq = max(1, int(mq))
            if gt == "RR":
                p = self.read_likelihood(
                    base,
                    ref,
                    bq,
                    mq
                )
            elif gt == "RA":
                p = (
                    0.5 * self.read_likelihood(base, ref, bq, mq)
                    + 0.5 * self.read_likelihood(base, alt, bq, mq)
                )
            else:  # AA
                p = self.read_likelihood(
                    base,
                    alt,
                    bq,
                    mq
                )
            if p <= 0:
                continue
            ll += math.log(p)
        return ll

    def genotype_prior(self, gt, het_rate=0.001):
        if gt == "RR":
            return math.log(1.0 - het_rate)
        if gt == "RA":
            return math.log(het_rate)
        return math.log(het_rate / 2.0)

    def call_variants_in_window(self, chrom, start1, end1, output_vcf=None, count_only=False):
        if not hasattr(self, "global_counts"):
            self.global_counts = defaultdict(Counter)
        if not hasattr(self, "global_records"):
            self.global_records = defaultdict(list)
        if not hasattr(self, "variants_seen"):
            with self.lock:
                if not hasattr(self, "variants_seen"):
                    self.variants_seen = set()
        import re
        pattern = re.compile(r'\[\+[^\]]+\]|\[\-[^\]]+\]|.')  # ✅ regex cache
        debug = bool(self.params.get("debug"))
        reads = self._reads_in_window(chrom, start1, end1)
        if not reads:
            return []
        counts, base_records, indel_positions = self._discover_sites_in_window(
            chrom, start1, end1, reads
        )
        br_get = base_records.get  # ✅ faster lookup
        candidate_positions = sorted(counts.keys())
        exclusion = int(self.params.get("indel_exclusion_window"))
        if exclusion < 0:
            exclusion = 0
        excluded = set()
        for ip in indel_positions:
            for x in range(max(1, ip - exclusion), ip + exclusion + 1):
                excluded.add(x)
        min_bq = int(self.params.get("min_base_quality"))
        min_mapq = int(self.params.get("min_mapping_quality"))
        min_coverage = int(self.params.get("min_coverage"))
        min_bq_local = min_bq
        min_mapq_local = min_mapq
        vcf_records = []
        processed_positions = set()
        i = 0
        ref_seq_full = self.ref_dict.get(chrom, "")
        if not ref_seq_full:
            return []
        # ✅ SW cache
        sw_cache = {}
        # ✅ project_read cache
        projected_cache = {}
        while i < len(candidate_positions):
            pos = candidate_positions[i]
            if pos in excluded or pos in processed_positions:
                i += 1
                continue
            group = [pos]
            j = i + 1
            while j < len(candidate_positions):
                nxt = candidate_positions[j]
                if nxt in excluded or nxt in processed_positions:
                    j += 1
                    continue
                if nxt - group[-1] <= 10:
                    group.append(nxt)
                    j += 1
                else:
                    break
            for g in group:
                processed_positions.add(g)
            i = j
            if len(group) > 8:
                group = group[:8]
            try:
                haplotypes, L, R = self._assemble_haplotypes(
                    chrom=chrom,
                    reads=reads,
                    core_positions=group,
                    window_start=start1,
                    window_end=end1,
                    region_start=None,
                    region_end=None
                )
                with self.lock:
                    self.processed_reads += 1
            except Exception as exc:
                raise RuntimeError(f"Error assembling haplotypes for {chrom}:{group} -> {exc}") from exc
            L_clamped = max(1, min(L, len(ref_seq_full)))
            R_clamped = max(L_clamped, min(R, len(ref_seq_full)))
            ref_region = ref_seq_full[L_clamped - 1:R_clamped]
            # ✅ expand haplotypes
            hapseqs_expanded = []
            for h in haplotypes:
                seq = h.get("seq", "")
                parts = pattern.findall(seq)
                expanded = []
                for t in parts:
                    if t.startswith("[+"):
                        expanded.append(t[2:-1])
                    elif t.startswith("[-"):
                        pass
                    else:
                        expanded.append(t)
                hapseqs_expanded.append("".join(expanded))
            allele_sets = []
            allele_mq_sums = []
            allele_bq_sums = []
            allele_supports = []
            for idx, h in enumerate(haplotypes):
                if idx == 0:
                    continue
                hapseq = hapseqs_expanded[idx]
                # ✅ SW cache
                if hapseq in sw_cache:
                    sw = sw_cache[hapseq]
                else:
                    sw = smith_waterman(ref_region, hapseq)
                    sw_cache[hapseq] = sw
                complex_var = detect_complex_variant(hapseq, ref_region)
                if complex_var:
                    ref_v, alt_v = complex_var
                if not isinstance(sw, dict) or "aligned_ref" not in sw:
                    continue
                aref = sw["aligned_ref"]
                aread = sw["aligned_read"]
                if not aref or not aread:
                    continue
                p_ref = int(sw.get("ref_start", 1))
                ref_pos_global = L_clamped + p_ref - 1
                k = 0
                alleles_here = []
                while k < len(aref) and k < len(aread):
                    rchar = aref[k]
                    dchar = aread[k]
                    if rchar == dchar:
                        ref_pos_global += 1
                        k += 1
                        continue
                    if rchar != "-" and dchar != "-":
                        alleles_here.append((ref_pos_global, rchar, dchar))
                        ref_pos_global += 1
                        k += 1
                        continue
                    if rchar != "-" and dchar == "-":
                        kk = k
                        while kk < len(aref) and aread[kk] == "-":
                            kk += 1
                        del_len = kk - k
                        refseq = ref_seq_full[ref_pos_global - 1: ref_pos_global + del_len]
                        if refseq:
                            altseq = refseq[0]
                            nref, nalt, npos = left_trim_ref_alt(refseq, altseq, ref_pos_global, ref_seq_full)
                            alleles_here.append((npos, nref, nalt))
                        ref_pos_global += del_len + 1
                        k = kk
                        continue
                    if rchar == "-" and dchar != "-":
                        kk = k
                        ins_seq = []
                        while kk < len(aref) and aref[kk] == "-":
                            ins_seq.append(aread[kk])
                            kk += 1
                        ins_seq = "".join(ins_seq)
                        if 1 <= ref_pos_global <= len(ref_seq_full):
                            ref_v = ref_seq_full[ref_pos_global - 1]
                            alt_v = ref_v + ins_seq
                            nref, nalt, npos = left_trim_ref_alt(ref_v, alt_v, ref_pos_global, ref_seq_full)
                            alleles_here.append((npos, nref, nalt))
                        k = kk
                        continue
                    k += 1
                uniq = {(p, r, a): (p, r, a) for (p, r, a) in alleles_here}
                allele_list = sorted(uniq.values(), key=lambda x: (x[0], x[1], x[2]))
                if allele_list:
                    allele_sets.extend(allele_list)
                    allele_mq_sums.append(int(h.get("mq_sum", 0)))
                    allele_bq_sums.append(int(h.get("bq_sum", 0)))
                    allele_supports.append(int(h.get("support", 0)))
            if count_only:
                with self.lock:
                    self.total_reads += len(reads)
                return []
            if not allele_sets:
                continue
            allele_dict = {}
            hap_idx = 0
            for idx, h in enumerate(haplotypes):
                if idx == 0:
                    continue  # skip reference haplotype
                mq = int(h.get("mq_sum", 0))
                bq = int(h.get("bq_sum", 0))
                supp = int(h.get("support", 0))
                hapseq = hapseqs_expanded[idx]
                sw = smith_waterman(ref_region, hapseq)
                if not sw or "aligned_ref" not in sw:
                    continue
                aref = sw["aligned_ref"]
                aread = sw["aligned_read"]
                if not aref or not aread:
                    continue
                p_ref = int(sw.get("ref_start", 1))
                ref_pos_global = L_clamped + p_ref - 1
                k = 0
                while k < len(aref) and k < len(aread):
                    rchar = aref[k]
                    dchar = aread[k]
                    if rchar == dchar:
                        ref_pos_global += 1
                        k += 1
                        continue
                    # SNP
                    if rchar != "-" and dchar != "-":
                        key = (ref_pos_global, rchar, dchar)
                        if key not in allele_dict:
                            allele_dict[key] = {"mq": 0, "bq": 0, "support": 0}
                        allele_dict[key]["mq"] += mq
                        allele_dict[key]["bq"] += bq
                        allele_dict[key]["support"] += supp
                        ref_pos_global += 1
                        k += 1
                        continue
                    # deletion
                    if rchar != "-" and dchar == "-":
                        start_pos = ref_pos_global
                        kk = k
                        while kk < len(aref) and aread[kk] == "-":
                            kk += 1
                        del_len = kk - k
                        ref_seq_vcf = ref_seq_full[start_pos - 1: start_pos - 1 + del_len + 1]
                        if ref_seq_vcf:
                            alt_seq_vcf = ref_seq_vcf[0]
                            nref, nalt, npos = left_trim_ref_alt(
                                ref_seq_vcf, alt_seq_vcf, start_pos, ref_seq_full
                            )
                            key = (npos, nref, nalt)
                            if key not in allele_dict:
                                allele_dict[key] = {"mq": 0, "bq": 0, "support": 0}
                            allele_dict[key]["mq"] += mq
                            allele_dict[key]["bq"] += bq
                            allele_dict[key]["support"] += supp
                        ref_pos_global += del_len + 1
                        k = kk
                        continue
                    # insertion
                    if rchar == "-" and dchar != "-":
                        kk = k
                        ins_seq = []
                        while kk < len(aref) and aref[kk] == "-" and aread[kk] != "-":
                            ins_seq.append(aread[kk])
                            kk += 1
                        ins_seq = "".join(ins_seq)
                        pos_anchor = ref_pos_global
                        if 1 <= pos_anchor <= len(ref_seq_full):
                            ref_v = ref_seq_full[pos_anchor - 1]
                            alt_v = ref_v + ins_seq
                            nref, nalt, npos = left_trim_ref_alt(
                                ref_v, alt_v, pos_anchor, ref_seq_full
                            )
                            ins_part_check = nalt[len(nref):].upper()
                            has_pileup_support = any(
                                rec["tok"].startswith("+")
                                and rec["tok"][1:].upper() == ins_part_check
                                for offset in range(-2, 3)
                                for rec in base_records.get(npos + offset, [])
                            )
                            if has_pileup_support:
                                key = (npos, nref, nalt)
                                if key not in allele_dict:
                                    allele_dict[key] = {"mq": 0, "bq": 0, "support": 0}
                                allele_dict[key]["mq"] += mq
                                allele_dict[key]["bq"] += bq
                                allele_dict[key]["support"] += supp
                        k = kk
                        continue
                    k += 1
            alleles_final = [(p, r, a, v["mq"], v["bq"], v["support"]) for (p, r, a), v in allele_dict.items()]
            alleles_final = sorted(alleles_final, key=lambda x: (x[0], x[1], x[2]))
            # DEBUG HERE
            # if debug:
            #     print("\n===== FINAL ALLELES BEFORE VCF =====")
            #     for p, r, a, mq, bq, supp in alleles_final:
            #         print(f"POS={p} REF={r} ALT={a} MQ_SUM={mq} BQ_SUM={bq} SUPPORT={supp}")
            # print("\n[DEBUG] alleles_final BEFORE VCF:")
            # for p, r, a, mq, bq, supp in alleles_final:
            #     vtype = "INS" if len(a) > len(r) else "DEL" if len(r) > len(a) else "SNP"
            #     print(f"  {vtype} pos={p} REF={r} ALT={a} mq={mq} bq={bq} supp={supp}")
            alleles_refalt = [(p, r, a) for p, r, a, mq, bq, s in alleles_final]
            # print("[DEBUG] alleles_refalt:", alleles_refalt)
            # print("[DEBUG] alleles_final:", alleles_final)
            allele_AD  = []
            allele_REF = []
            for p, r, a in alleles_refalt:
                is_ins = len(a) > len(r)
                is_del = len(r) > len(a)
                # For insertions: anchor may be stored at p-1 or p-2 due to
                # left-normalisation, so search a small window.
                if is_ins:
                    pile = []
                    for offset in range(-2, 3):
                        pile.extend(base_records.get(p + offset, []))
                else:
                    pile = base_records.get(p, [])
                ref_count   = 0
                alt_count   = 0
                alt_forward = 0
                alt_reverse = 0
                for rec in pile:
                    tok    = rec["tok"]
                    bq     = rec["bq"]
                    mq     = rec["mq"]
                    strand = rec.get("strand", "+")
                    if isinstance(bq, str):
                        bq = ord(bq) - 33
                    if isinstance(mq, str):
                        mq = ord(mq) - 33
                    if bq < min_bq or mq < min_mapq:
                        continue
                    # SNV
                    if len(tok) == 1:
                        if tok.upper() == r.upper():
                            ref_count += 1
                        elif not is_ins and tok.upper() == a.upper():
                            alt_count += 1
                            if strand == "+":
                                alt_forward += 1
                            else:
                                alt_reverse += 1
                    # INSERTION
                    # Token format from _discover_sites_in_window: "+seq"
                    # inserted bases = a[len(r):]
                    elif tok.startswith("+") and is_ins:
                        inserted = tok[1:].upper()
                        ins_part = a[len(r):].upper()
                        if inserted == ins_part:
                            alt_count += 1
                            if strand == "+":
                                alt_forward += 1
                            else:
                                alt_reverse += 1
                    # DELETION 
                    # Token format: "-seq"  deleted bases = r[len(a):]
                    elif tok.startswith("-") and is_del:
                        del_tok  = tok[1:].upper()
                        del_part = r[len(a):].upper()
                        if del_tok == del_part:
                            alt_count += 1
                            if strand == "+":
                                alt_forward += 1
                            else:
                                alt_reverse += 1           
                if is_ins and alt_count == 0:
                    ins_part = a[len(r):].upper()   # inserted bases only
                    anchor_base = r[0].upper()
                    #  once per allele, NOT inside the read loop              
                    # print(f"  [PRE-FALLBACK] pos={p} REF={r} ALT={a} "
                    #     f"alt_count={alt_count} triggering_fallback=True")
                    for read in reads:
                        if read.get("mapq", 0) < min_mapq:
                            continue
                        strand = read.get("strand", "+")
                        proj, insertions = self._project_read(read, ref_seq_full)
                        #  Log nearby insertions found in this read 
                        if insertions:
                            for anchor, ins_seq, avgq, mq_val in insertions:
                                if abs(anchor - p) <= 10:
                                    print(f"  [FALLBACK DEBUG] pos={p} ins_part={ins_part} "
                                        f"anchor={anchor} ins_seq={ins_seq}")
                        # Check insertion list for this read
                        for anchor, ins_seq, avgq, mq_val in insertions:
                            if abs(anchor - p) > 2:
                                continue
                            if isinstance(avgq, str):
                                avgq = ord(avgq) - 33
                            if isinstance(mq_val, str):
                                mq_val = ord(mq_val) - 33
                            if avgq < min_bq or mq_val < min_mapq:
                                continue
                            if ins_seq.upper() == ins_part:
                                alt_count += 1
                                if strand == "+":
                                    alt_forward += 1
                                else:
                                    alt_reverse += 1
                                break  # one insertion hit per read max                           
                allele_AD.append(alt_count)
                allele_REF.append(ref_count)
                vtype = "INS" if is_ins else "DEL" if is_del else "SNP"
                # print(
                #     f"[DEBUG] {vtype} pos={p} REF={r} ALT={a} "
                #     f"pile={len(pile)} "
                #     f"ref={ref_count} alt={alt_count} "
                #     f"fwd={alt_forward} rev={alt_reverse}"
                # )
            # total_dp = sum(allele_REF) + sum(allele_AD)
            total_dp = len(base_records.get(p, []))
            if total_dp < min_coverage:
                continue
            # Filter alleles
            filtered_alleles = []
            filtered_AD = []
            filtered_REF = []
            for idx, (p, r, a, mq, bq, supp) in enumerate(alleles_final):
                ad = allele_AD[idx]
                refc = allele_REF[idx]
                # Allow SNPs and INDELs
                if r != a and ad >= 1:
                    filtered_alleles.append((p, r, a))
                    filtered_AD.append(ad)
                    filtered_REF.append(refc)
            # Genotyping (same as yours)
            alleles_refalt_filtered = filtered_alleles
            gl_map, allele_mq_sum, allele_bq_sum = self._compute_genotype_likelihoods_from_haplotypes(
                alleles_refalt_filtered, reads, chrom, L_clamped, R_clamped
            )
            ploidy = int(self.params.get("ploidy"))
            n_alleles = 1 + len(alleles_refalt_filtered)
            gorder = genotype_order(n_alleles, ploidy)
            mx = max(gl_map.values())
            probs = {gt: math.exp(gl_map.get(gt, -1e300) - mx) for gt in gorder}
            total = sum(probs.values())
            probs = {gt: p/total for gt,p in probs.items()}
            raw_pls = phred_scale_from_probs([probs.get(gt, 1e-300) for gt in gorder])
            pls = []
            for x in raw_pls:
                try:
                    if isinstance(x, (int, float)):
                        pls.append(int(x))
                    elif isinstance(x, str):
                        pls.append(int(x.split(",")[0]))
                    else:
                        pls.append(0)
                except:
                    pls.append(0)
            best_idx = int(np.argmin(pls))
            best_gt = gorder[best_idx]
            gt_str = "/".join(map(str, best_gt))
            total_dp_out = int(total_dp) if total_dp > 0 else 0
            allref = tuple([0] * ploidy)
            if allref in gorder:
                ref_idx = gorder.index(allref)
                qual_val = max(0, pls[ref_idx] - pls[best_idx])
            else:
                qual_val = 30
            alleles_by_pos = defaultdict(list)
            for idx, (p, r, a) in enumerate(filtered_alleles):
                alleles_by_pos[(p, r)].append((a, filtered_AD[idx]))
            for (p, r), alt_list in alleles_by_pos.items():
                alts = [a for a,_ in alt_list]
                acs  = [ac for _,ac in alt_list]
                alts = [str(a).replace(",", "") for a in alts]
                alt_str = ",".join(alts)
                afs = [
                    ac/total_dp_out if total_dp_out>0 else 0
                    for ac in acs
                ]
                key = (chrom, p, r, alt_str)
                with self.lock:
                    if key in self.variants_seen:
                        continue
                    self.variants_seen.add(key)
                safe_acs = [str(int(x)) for x in acs]
                safe_afs = [f"{float(x):.6f}" for x in afs]
                types = set([get_variant_type(r, a) for a in alts])
                if len(types) == 1:
                    vtype = list(types)[0]
                else:
                    vtype = "mixed"
                info_str = ";".join([
                    f"DP={int(total_dp_out)}",
                    f"AC={','.join(safe_acs)}",
                    f"AF={','.join(safe_afs)}",
                    f"SAF={int(alt_forward)}",
                    f"SAR={int(alt_reverse)}",
                    f"TYPE={vtype}"
                ])
                pos_vcf = p  # already 1-based in your code
                sample_pl = ",".join(str(int(p)) for p in pls)
                record = "\t".join([
                    chrom,
                    str(pos_vcf),
                    ".",
                    r,
                    alt_str,
                    str(int(qual_val)),
                    "PASS",
                    info_str,
                    "GT:PL",
                    f"{gt_str}:{sample_pl}"
                ])
                vcf_records.append(record)
                with self.lock:
                    self.completed_haplotype_regions += 1
                if debug:
                    print(f"\\===== EXTRACTED ALLELES FOR REGION {chrom}:{L}-{R} =====")
                    for p, r, a in alleles_refalt:
                        vtype = "SNP"
                        if len(r) > 1 or len(a) > 1:
                            if len(r) > len(a):
                                vtype = f"DEL({len(r)-len(a)}bp)"
                            elif len(a) > len(r):
                                vtype = f"INS({len(a)-len(r)}bp)"
                            else:
                                vtype = f"MNP({len(r)}bp)"
                        print(f"  POS={p} TYPE={vtype} REF={r} ALT={a}")
        return vcf_records

    # VCF header writer
    def _vcf_header(self):
        hdr = []
        hdr.append("##fileformat=VCFv4.2")
        hdr.append(f"##source=pyfreebayes_trial4_haplobuild")
        hdr.append(f"##reference={self.reference}")
        hdr.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
        hdr.append('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">')
        hdr.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        hdr.append('##INFO=<ID=MQ_SUM,Number=A,Type=Integer,Description="Per-allele sum of mapping qualities">')
        hdr.append('##INFO=<ID=BQ_SUM,Number=A,Type=Integer,Description="Per-allele sum of base qualities">')
        hdr.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        hdr.append('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">')
        hdr.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + self.sample_name)
        return "\n".join(hdr) + "\n"

    def _update_progress(self, step_name="Processing"):
        if not self.progress_callback or self.total_haplotype_regions == 0:
            return
        percent = int((self.completed_haplotype_regions / self.total_haplotype_regions) * 100)
        percent = max(0, min(100, percent))
        if getattr(self, "_last_percent", None) == percent:
            return  # skip redundant updates
        self._last_percent = percent
        self.progress_callback(percent, step_name)

    def run(self, output_vcf=None):
        from concurrent.futures import ThreadPoolExecutor
        import time
        import traceback
        windows = self.make_windows()
        vcf_files = []
        # Progress bookkeeping
        self.completed_haplotype_regions = 0
        self.total_haplotype_regions = 0
        self.total_reads = 0
        self.processed_reads = 0
        total_windows = len(windows)   
        for idx, (chrom, s, e) in enumerate(windows, start=1):
            self.call_variants_in_window(chrom, s, e, count_only=True)
            percent = int((idx / total_windows) * 100)
            if self.progress_callback:
                self.progress_callback(
                    percent,
                    f"Initializing variant calling ({idx}/{total_windows})"
                )
            time.sleep(0.02)
        if self.progress_callback:
            self.progress_callback(0, "Initializing FreeBayes")

        def worker(chrom, s, e, outpath):
            try:
                recs = self.call_variants_in_window(chrom, s, e)
                outdir = Path(outpath).parent
                outdir.mkdir(parents=True, exist_ok=True)
                with open(outpath, "w") as fo:
                    fo.write(self._vcf_header())
                    if recs:
                        for r in recs:
                            fo.write(r + "\n")
                return str(outpath)
            except Exception as exc:
                err_msg = f"Exception in worker for window {chrom}:{s}-{e} -> {exc!r}"
                print("[ERROR] " + err_msg)
                traceback.print_exc()
                raise RuntimeError(err_msg) from exc
        # ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=self.threads) as ex:
            futures = [
                ex.submit(
                    worker,
                    chrom, s, e,
                    str(self.tempdir / f"{chrom}_{s}_{e}.vcf")
                )
                for chrom, s, e in windows
                    ]
            # MAIN THREAD progress updater
            while any(not fut.done() for fut in futures):
                with self.lock:
                    done = self.processed_reads
                    total = max(1, self.total_reads)
                percent = int((done / total) * 100)
                if self.progress_callback:
                    self.progress_callback(
                        percent,
                        f"Processing haplotype region: {done} haplotypes /{total} reads"
                    )
                time.sleep(0.1)
            # Collect results AFTER completion
            for fut in futures:
                try:
                    res = fut.result()
                    if res:
                        vcf_files.append(res)
                except Exception as exc:
                    print(f"[ERROR] worker failed -> {exc}")
                    traceback.print_exc()
        # Merge VCFs
        merged_path = Path(self.tempdir) / f"{self.sample_name}_merged.vcf"
        lines = []
        for p in vcf_files:
            try:
                with open(p, "r") as fh:
                    for l in fh:
                        if not l.startswith("#"):
                            lines.append(l.rstrip("\n"))
            except Exception as exc:
                print(f"[WARN] could not read {p}: {exc!r}")

        def keyfunc(line):
            parts = line.split("\t")
            chrom = parts[0]
            try:
                pos_raw = parts[1]
                if "," in pos_raw:
                    pos_raw = pos_raw.split(",")[0]
                pos = int(pos_raw)
            except Exception:
                pos = 0
            return (chrom, pos)
        lines.sort(key=keyfunc)
        with open(merged_path, "w") as fo:
            fo.write(self._vcf_header())
            for l in lines:
                fo.write(l + "\n")
        out_final = merged_path
        if output_vcf:
            try:
                shutil.move(str(merged_path), output_vcf)
                out_final = Path(output_vcf)
            except Exception as exc:
                print(f"[WARN] move failed: {exc!r}")
        # Final UI update (main thread)
        if self.progress_callback:
            self.progress_callback(100, "Variant Calling Completed")
        return str(out_final)

def get_variant_type(ref, alt):
    """Helper to determine variant type for VCF INFO field."""
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    elif len(ref) == len(alt) and len(ref) > 1:
        return "mnp"
    elif len(alt) > len(ref):
        return "ins"
    elif len(ref) > len(alt):
        return "del"
    else:
        return "complex"
