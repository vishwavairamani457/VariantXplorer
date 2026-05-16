# Alignment Backend Code
import os
import struct
import re, shutil
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import bamnostic.bgzf as bgzf
from bamnostic.bgzf import BgzfWriter, BgzfReader

def sais(text):
    if not text.endswith("$"):
        text += "$"
    n = len(text)
    sa = list(range(n))
    rank = [ord(c) for c in text]
    tmp = [0] * n
    k = 1
    while k < n:
        sa.sort(key=lambda i: (rank[i], rank[i + k] if i + k < n else -1))
        tmp[sa[0]] = 0
        for i in range(1, n):
            prev = sa[i - 1]
            curr = sa[i]
            tmp[curr] = tmp[prev] + (
                (rank[prev], rank[prev + k] if prev + k < n else -1) <
                (rank[curr], rank[curr + k] if curr + k < n else -1)
            )
        rank = tmp[:]
        k <<= 1
    return sa

# COMPRESSED FM-INDEX
class FMIndex:
    def __init__(self, reference, checkpoint=256):
        self.reference = reference
        self.checkpoint = checkpoint
        self.sa = sais(reference)
        self.bwt = self._build_bwt()
        self.C = self._build_C()
        self.occ = self._build_occ()

    def _build_bwt(self):
        text = self.reference + "$"
        bwt = []
        for i in self.sa:
            bwt.append("$" if i == 0 else text[i - 1])
        return "".join(bwt)

    def _build_C(self):
        counts = {}
        for c in self.bwt:
            counts[c] = counts.get(c, 0) + 1

        total = 0
        C = {}
        for c in sorted(counts):
            C[c] = total
            total += counts[c]
        return C

    def _build_occ(self):
        self.checkpoint = 256
        alphabet = sorted(set(self.bwt))
        occ = {c: [] for c in alphabet}
        running = {c: 0 for c in alphabet}
        for i, char in enumerate(self.bwt):
            if i % self.checkpoint == 0:
                for c in alphabet:
                    occ[c].append(running[c])
            running[char] += 1
        return occ

    def _occ_count(self, char, pos):
        if char not in self.occ:
            return 0
        chk = min(pos // self.checkpoint, len(self.occ[char]) - 1)
        count = self.occ[char][chk]
        start = chk * self.checkpoint
        for i in range(start, pos):
            if self.bwt[i] == char:
                count += 1
        return count

    def search(self, pattern):
        l = 0
        r = len(self.bwt)
        for char in reversed(pattern):
            if char not in self.C:
                return []
            l = self.C[char] + self._occ_count(char, l)
            r = self.C[char] + self._occ_count(char, r)
            if l >= r:
                return []
        return self.sa[l:r]

# MINIMIZER SEEDING
from collections import deque
def compute_minimizers(seq, k=21, w=20):
    n = len(seq)
    if n < k:
        return []
    kmers = [seq[i:i+k] for i in range(n - k + 1)]
    minimizers = []
    dq = deque()
    for i, kmer in enumerate(kmers):
        while dq and dq[-1][0] >= kmer:
            dq.pop()
        dq.append((kmer, i))
        # Remove out-of-window
        while dq and dq[0][1] <= i - w:
            dq.popleft()
        if i >= w - 1:
            minimizers.append(dq[0])
    # Remove duplicates
    seen = set()
    result = []
    for kmer, pos in minimizers:
        if pos not in seen:
            seen.add(pos)
            result.append((kmer, pos))
    return result


# FAST K-MER HASH INDEX
class KmerIndex:
    """
    Stores all k-mer positions from the reference in a dict.
    For large genomes (> 5 MB) this replaces both the FM-index and
    the O(n) str.find() fallback, giving O(1) seed lookups.
    """
    def __init__(self, reference: str, k: int = 21, max_occ: int = 500):
        self.k = k
        self.max_occ = max_occ
        self._index: dict[str, list[int]] = defaultdict(list)
        self._build(reference)

    def _build(self, reference: str):
        k = self.k
        n = len(reference)
        for i in range(n - k + 1):
            kmer = reference[i:i + k]
            lst = self._index[kmer]
            # Cap repetitive k-mers early to avoid memory explosion
            if len(lst) < self.max_occ:
                lst.append(i)

    def search(self, kmer: str) -> list[int]:
        return self._index.get(kmer, [])

class SimpleBWAligner:
    # BAM nucleotide encoding (as per BAM spec 4-bit codes)
    nt_to_code = {
        "=": 0, "A": 1, "C": 2, "M": 3, "G": 4, "R": 5, "S": 6, "V": 7,
        "T": 8, "W": 9, "Y": 10, "H": 11, "K": 12, "D": 13, "B": 14, "N": 15
    }
    # inverse for decoding
    code_to_nt = {v: k for k, v in nt_to_code.items()}
    # BAM CIGAR op -> int
    op_to_int = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    int_to_op = {v: k for k, v in op_to_int.items()}

    def _normalize_dict(self, obj):
        if isinstance(obj, dict):
            return obj
        if hasattr(obj, "__dict__"):
            return dict(obj.__dict__)
        return {}

    def get_reference_params(self, default_params=None):
        defaults = {
            "Mismatch Penalty": 4,
            "Threads": 8,
            "Seed length": 20,
            "Band Width": 100,
            "Z-dropoff": 100,
            "Seed Split ratio": 1.5,
            "Max Occurence": 10000,
            "Match Score": 1,
            "Gap Open Penalty": 6,
            "Gap Extension Penalty": 1,
            "Clipping Penalty": 5,
            "Min Map Score": 1,
            "LowComplexityThreshold": 0.8,
        }
        params = defaults.copy()
        incoming = self._normalize_dict(default_params) or {}
        nested = {}
        if isinstance(incoming, dict) and "Reference Mapping" in incoming and isinstance(incoming["Reference Mapping"], dict):
            nested = incoming["Reference Mapping"]
        if isinstance(incoming, dict):
            for k, v in incoming.items():
                if k == "Reference Mapping":
                    continue
                if k in params:
                    params[k] = v
        for k, v in nested.items():
            if k in params:
                params[k] = v

        def _coerce(key, cast, fallback):
            val = params.get(key, fallback)
            try:
                return cast(val)
            except Exception:
                return fallback

        params["Mismatch Penalty"] = _coerce("Mismatch Penalty", int, defaults["Mismatch Penalty"])
        params["Threads"] = _coerce("Threads", int, defaults["Threads"])
        params["Seed length"] = _coerce("Seed length", int, defaults["Seed length"])
        params["Band Width"] = _coerce("Band Width", int, defaults["Band Width"])
        params["Z-dropoff"] = _coerce("Z-dropoff", int, defaults["Z-dropoff"])
        params["Seed Split ratio"] = _coerce("Seed Split ratio", float, defaults["Seed Split ratio"])
        params["Max Occurence"] = _coerce("Max Occurence", int, defaults["Max Occurence"])
        params["Match Score"] = _coerce("Match Score", float, defaults["Match Score"])
        params["Gap Open Penalty"] = _coerce("Gap Open Penalty", int, defaults["Gap Open Penalty"])
        params["Gap Extension Penalty"] = _coerce("Gap Extension Penalty", int, defaults["Gap Extension Penalty"])
        params["Clipping Penalty"] = _coerce("Clipping Penalty", int, defaults["Clipping Penalty"])
        params["Min Map Score"] = _coerce("Min Map Score", float, defaults["Min Map Score"])
        params["LowComplexityThreshold"] = _coerce("LowComplexityThreshold", float, defaults["LowComplexityThreshold"])
        self.params = params
        return params

    def __init__(self, reference_fasta, params=None):
        self.params = params.copy() if isinstance(params, dict) else {}
        self.ref_dict = self._load_reference_dict(reference_fasta)
        if not self.ref_dict:
            raise ValueError("Reference FASTA empty or not found.")
        self.ref_names = list(self.ref_dict.keys())
        self.reference = self.ref_dict[self.ref_names[0]]
        normalized_defaults = self._normalize_dict(self.params)
        self.get_reference_params(normalized_defaults)
        if isinstance(params, dict) and "Reference Mapping" in params:
            self.update_params(params=params)
        if isinstance(params, dict):
            for k, v in params.items():
                if k in self.params:
                    try:
                        self.params[k] = type(self.params[k])(v)
                    except:
                        self.params[k] = v
        self.params.setdefault("Band Width", 100)
        if len(self.reference) > 5_000_000:
            bw = int(self.params.get("Band Width", 100))
            self.params["Band Width"] = min(bw, 50)
        self._comp_table = str.maketrans("ATCGUNatcgun", "TAGCANtagcan")
        self.revcomp = self._revcomp_fast

        # INDEX SELECTION
        # Small genome  → FM-index (exact, memory-efficient)
        # Large genome  → KmerIndex (O(1) lookup, built once, fast for WGS)
        seed_k = int(self.params.get("Seed length", 21))
        max_occ = int(self.params.get("Max Occurence", 500))
        if len(self.reference) < 5_000_000:
            print("Small genome -> building FM-index")
            self.fm = FMIndex(self.reference)
            self.kmer_index = None
            self.use_fm = True
        else:
            print("🧬 Large genome detected → building k-mer hash index (fast seeding)")
            self.fm = None
            # Use seed_k for the kmer index so seeds match params
            self.kmer_index = KmerIndex(self.reference, k=seed_k, max_occ=max_occ)
            self.use_fm = False

    def _revcomp_fast(self, seq: str) -> str:
        return seq.translate(self._comp_table)[::-1]

    def update_params(self, default_params=None, params=None, debug=False):
        if default_params:
            default_params = self._normalize_dict(default_params)
        else:
            default_params = {}
        if params and isinstance(params, dict):
            for k, v in params.items():
                if k in self.params:
                    try:
                        self.params[k] = type(self.params[k])(v)
                    except Exception:
                        self.params[k] = v
        if "Reference Mapping" in default_params and isinstance(default_params["Reference Mapping"], dict):
            ref_params = default_params["Reference Mapping"]
        else:
            ref_params = {}
        for key in list(self.params.keys()):
            if key in ref_params:
                try:
                    self.params[key] = type(self.params[key])(ref_params[key])
                except Exception:
                    pass
        if debug:
            print("Updated params:", self.params)
        return self.params

    def _load_reference_dict(self, fasta_file):
        ref_dict = {}
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"{fasta_file} not found")
        name = None
        seqs = []
        with open(fasta_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if name:
                        ref_dict[name] = "".join(seqs)
                    name = line[1:].strip().split()[0]
                    seqs = []
                else:
                    seqs.append(line.strip())
            if name:
                ref_dict[name] = "".join(seqs)
        return ref_dict

    def is_low_complexity(self, seq):
        if not seq:
            return True
        from collections import Counter
        freq = Counter(seq)
        maxf = max(freq.values()) / len(seq)
        return maxf >= float(self.params.get("LowComplexityThreshold"))

    def align_and_cigar(self, read, ref_subseq):
        params = self.params
        match = float(params["Match Score"])
        mismatch = -float(params["Mismatch Penalty"])
        gap_open = -float(params["Gap Open Penalty"])
        gap_extend = -float(params["Gap Extension Penalty"])
        z_drop = float(params["Z-dropoff"])
        band = min(int(params["Band Width"]), max(20, int(len(read) * 0.05)))      
        m = len(read)
        n = len(ref_subseq)
        if m == 0 or n == 0:
            return 0, "*", 0, 0
        NEG_INF = float("-inf")
        H_prev = [0] * (n + 1)
        H_curr = [0] * (n + 1)
        E = [NEG_INF] * (n + 1)
        F = [NEG_INF] * (n + 1)
        traceback = {}
        best_score = 0
        best_pos = (0, 0)
        early_stop_score = m * match * 0.9
        for i in range(1, m + 1):
            j_start = max(1, i - band)
            j_end = min(n, i + band)
            row_best = 0
            F[j_start - 1] = NEG_INF
            for j in range(j_start, j_end + 1):
                score_sub = match if read[i-1] == ref_subseq[j-1] else mismatch
                E[j] = max(H_prev[j] + gap_open, E[j] + gap_extend)
                F[j] = max(H_curr[j-1] + gap_open, F[j-1] + gap_extend)
                diag = H_prev[j-1] + score_sub
                h_val = max(0, diag, E[j], F[j])
                if best_score > 0 and (best_score - h_val) > z_drop:
                    H_curr[j] = 0
                    continue
                H_curr[j] = h_val
                row_best = max(row_best, h_val)
                if h_val > 0:
                    if h_val == diag:
                        traceback[(i, j)] = 1
                    elif h_val == E[j]:
                        traceback[(i, j)] = 2
                    elif h_val == F[j]:
                        traceback[(i, j)] = 3
                if h_val > best_score:
                    best_score = h_val
                    best_pos = (i, j)
            if best_score >= early_stop_score:
                break
            if row_best == 0 and best_score > 0:
                break
            H_prev, H_curr = H_curr, [0] * (n + 1)
        if best_score <= 0:
            return 0, "*", 0, 0
        i, j = best_pos
        cigar = []
        ref_consumed = 0
        while i > 0 and j > 0:
            key = (i, j)
            if key not in traceback:
                break
            tb = traceback[key]
            if tb == 1:
                cigar.append("M")
                i -= 1
                j -= 1
                ref_consumed += 1
            elif tb == 2:
                cigar.append("I")
                i -= 1
            elif tb == 3:
                cigar.append("D")
                j -= 1
                ref_consumed += 1
            else:
                break
        ref_start_offset = j
        cigar.reverse()
        cigar_ops = []
        prev = None
        count = 0
        for c in cigar:
            if c == prev:
                count += 1
            else:
                if prev:
                    cigar_ops.append(f"{count}{prev}")
                prev = c
                count = 1
        if prev:
            cigar_ops.append(f"{count}{prev}")
        cigar_str = "".join(cigar_ops) if cigar_ops else "*"
        return best_score, cigar_str, ref_start_offset, ref_consumed

    def _kmer_search(self, kmer: str) -> list[int]:
        """Unified k-mer lookup: uses KmerIndex for large genomes, FM-index for small."""
        if self.use_fm:
            return self.fm.search(kmer)
        else:
            return self.kmer_index.search(kmer)

    def seed_positions(self, read):
        seed_len = int(self.params.get("Seed length"))
        max_occ = int(self.params.get("Max Occurence"))
        split_ratio = float(self.params.get("Seed Split ratio"))
        if len(read) < seed_len:
            return []
        positions = []
        step = seed_len  # non-overlapping seeds — fastest
        # Primary seeding
        for i in range(0, len(read) - seed_len + 1, step):
            seed = read[i:i + seed_len]
            hits = self._kmer_search(seed)
            if hits:
                for h in hits[:max_occ]:
                    start = h - i
                    if start >= 0:
                        positions.append((start, '+'))
            rc = self.revcomp(seed)
            hits = self._kmer_search(rc)
            if hits:
                for h in hits[:max_occ]:
                    start = h - (len(read) - i - seed_len)
                    if start >= 0:
                        positions.append((start, '-'))
            if len(positions) > max_occ:
                break
        # Fallback seeds (shorter k-mers)
        if not positions:
            for new_seed_len in [
                max(12, int(seed_len / split_ratio)),
                max(10, int(seed_len / (split_ratio * 2))),
                8, 6
            ]:
                if new_seed_len >= seed_len:
                    continue
                step2 = max(1, new_seed_len // 2)
                for i in range(0, len(read) - new_seed_len + 1, step2):
                    seed = read[i:i + new_seed_len]
                    if self.use_fm:
                        hits = self.fm.search(seed)
                    else:
                        # Short fallback: use str.find once (still O(n) but rare path)
                        pos = self.reference.find(seed)
                        hits = [pos] if pos != -1 else []
                    if hits:
                        for h in hits[:max_occ]:
                            start = h - i
                            if start >= 0:
                                positions.append((start, '+'))
                    rc = self.revcomp(seed)
                    if self.use_fm:
                        hits = self.fm.search(rc)
                    else:
                        pos = self.reference.find(rc)
                        hits = [pos] if pos != -1 else []
                    if hits:
                        for h in hits[:max_occ]:
                            start = h - (len(read) - i - new_seed_len)
                            if start >= 0:
                                positions.append((start, '-'))
                    if len(positions) > max_occ:
                        break
                if positions:
                    break
        # Remove duplicates + cluster
        uniq = {}
        for p, s in positions:
            if p not in uniq:
                uniq[p] = s
        result = sorted([(p, uniq[p]) for p in uniq.keys()])
        clustered = []
        CLUSTER_DIST = 20
        for pos, strand in result:
            if not clustered or abs(pos - clustered[-1][0]) > CLUSTER_DIST:
                clustered.append((pos, strand))
        MAX_CANDIDATES = 30
        return clustered[:min(MAX_CANDIDATES, max_occ)]
  
    def best_alignment(self, read):
        quick_pos = self.reference.find(read[:20])
        if quick_pos != -1:
            return quick_pos, 60, '+', f"{len(read)}M"        
        if self.is_low_complexity(read):
            return -1, 0, '+', "*"
        read_len = len(read)
        if read_len > 80000:
            return -1, 0, '+', "*"
        seed_hits = self.seed_positions(read)
        if not seed_hits:
            return -1, 0, '+', "*"
        scored = []
        for pos, strand in seed_hits:
            if pos < 0 or pos + read_len >= len(self.reference):
                continue
            read_seq = self.revcomp(read) if strand == '-' else read
            ref_seg = self.reference[pos:pos + read_len]
            matches = sum(1 for a, b in zip(read_seq, ref_seg) if a == b)
            scored.append((matches, pos, strand))
        scored.sort(reverse=True)
        top_candidates = [(p, s) for _, p, s in scored[:5]]
        if not top_candidates:
            return -1, 0, '+', "*"
        best_score = -1
        best_pos = -1
        best_strand = '+'
        best_cigar = "*"
        band = min(50, max(20, read_len // 20))
        for pos, strand in top_candidates:
            start = max(0, pos - band)
            end = min(len(self.reference), pos + read_len + band)
            ref_subseq = self.reference[start:end]
            read_seq = self.revcomp(read) if strand == '-' else read
            if read_seq in ref_subseq:
                abs_pos = start + ref_subseq.index(read_seq)
                return abs_pos, 60, strand, f"{read_len}M"
            score, cigar, ref_offset, _ = self.align_and_cigar(read_seq, ref_subseq)
            if score > best_score:
                best_score = score
                best_pos = start + ref_offset
                best_strand = strand
                best_cigar = cigar
            if score > read_len * 0.9:
                break
        if best_score <= 0:
            return -1, 0, '+', "*"
        return best_pos, 60, best_strand, best_cigar

    def read_fastq_chunks(self, fq_file, chunk_size=10000):
        chunk = []
        with open(fq_file, "r") as fq:
            while True:
                h = fq.readline()
                if not h:
                    if chunk:
                        yield chunk
                    break
                header = h.strip()
                seq = fq.readline().strip()
                fq.readline()  # plus line
                qual = fq.readline().strip()
                chunk.append((header[1:] if header.startswith("@") else header, seq, qual))
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
    @staticmethod
    def _estimate_read_count(fq_file: str) -> int:
        file_size = os.path.getsize(fq_file)
        total_bytes = 0
        records = 0
        with open(fq_file, "rb") as f:
            for _ in range(400):
                line = f.readline()
                if not line:
                    break
                total_bytes += len(line)
                records += 1
        if records < 4:
            return max(1, file_size // 200)   # very rough fallback
        avg_record_bytes = total_bytes / (records / 4)  # 4 lines per FASTQ record
        return max(1, int(file_size / avg_record_bytes))

    def align_reads_single(self, fq_file, output_sam, progress_callback=None):
        total_est_reads = self._estimate_read_count(fq_file)
        total_processed = 0
        threads = int(self.params.get("Threads"))
        with open(output_sam, "w") as sam:
            sam.write("@HD\tVN:1.0\tSO:unsorted\n")
            sam.write(f"@SQ\tSN:{self.ref_names[0]}\tLN:{len(self.reference)}\n")
            with ThreadPoolExecutor(max_workers=threads) as ex:
                for chunk in self.read_fastq_chunks(fq_file, 10000):
                    results = [None] * len(chunk)
                    def process(idx, read_tuple):
                        qname, seq, qual = read_tuple
                        pos, mapq, strand, cigar = self.best_alignment(seq)
                        return idx, (qname, pos, mapq, strand, seq, qual, cigar)
                    futures = {
                        ex.submit(process, i, r): i
                        for i, r in enumerate(chunk)
                    }
                    for fut in as_completed(futures):
                        i, result = fut.result()
                        results[i] = result
                        total_processed += 1
                        if progress_callback and total_processed % 100 == 0:
                            percent = int((total_processed / total_est_reads) * 100)
                            percent = max(0, min(99, percent))
                            progress_callback(
                                percent,
                                f"Aligning {total_processed}/{total_est_reads} reads..."
                            )
                    for qname, pos, mapq, strand, seq, qual, cigar in results:
                        if pos < 0:
                            flag, pos_out, cigar_out, mapq_out = 4, 0, "*", 0
                        else:
                            flag = 0 if strand == '+' else 16
                            pos_out = pos + 1
                            cigar_out = cigar if cigar != "*" else f"{len(seq)}M"
                            mapq_out = mapq
                        sam.write(
                            f"{qname}\t{flag}\t{self.ref_names[0]}"
                            f"\t{pos_out}\t{mapq_out}\t{cigar_out}"
                            f"\t*\t0\t0\t{seq}\t{qual}\n"
                        )
                    if progress_callback:
                        percent = int((total_processed / total_est_reads) * 100)
                        progress_callback(
                            percent,
                            f"Aligning {total_processed}/{total_est_reads} reads..."
                        )
        if progress_callback:
            progress_callback(100, f"Finalizing {total_processed} reads...")
        self.convert_sam_to_bam(output_sam)
        return True
    def align_reads_paired(self, fq1, fq2, output_sam, progress_callback=None):
        total_est_reads = min(
            self._estimate_read_count(fq1),
            self._estimate_read_count(fq2)
        )
        total_processed = 0
        last_printed = 0
        it1 = self.read_fastq_chunks(fq1, 10000)
        it2 = self.read_fastq_chunks(fq2, 10000)
        with open(output_sam, "w") as sam:
            sam.write("@HD\tVN:1.0\tSO:unsorted\n")
            sam.write(f"@SQ\tSN:{self.ref_names[0]}\tLN:{len(self.reference)}\n")
            while True:
                try:
                    chunk1 = next(it1)
                    chunk2 = next(it2)
                except StopIteration:
                    break
                n = min(len(chunk1), len(chunk2))
                chunk1 = chunk1[:n]
                chunk2 = chunk2[:n]
                results = [None] * n
                def process_pair(idx, pair):
                    (qname1, seq1, qual1), (qname2, seq2, qual2) = pair
                    pos1, mapq1, strand1, cigar1 = self.best_alignment(seq1)
                    pos2, mapq2, strand2, cigar2 = self.best_alignment(seq2)
                    return idx, (
                        qname1, seq1, qual1, pos1, mapq1, strand1, cigar1,
                        qname2, seq2, qual2, pos2, mapq2, strand2, cigar2
                    )
                with ThreadPoolExecutor(max_workers=int(self.params.get("Threads", 8))) as ex:
                    futures = {ex.submit(process_pair, i, (chunk1[i], chunk2[i])): i for i in range(n)}
                    for fut in as_completed(futures):
                        idx, res = fut.result()
                        results[idx] = res
                        total_processed += 1
                        if progress_callback and total_processed - last_printed >= 100:
                            last_printed = total_processed
                            percent = int((total_processed / total_est_reads) * 100)
                            percent = max(0, min(99, percent))
                            progress_callback(percent,
                                f"Aligning {total_processed}/{total_est_reads} paired reads...")
                for (q1, s1, ql1, p1, mq1, st1, cg1,
                    q2, s2, ql2, p2, mq2, st2, cg2) in results:
                    flag1 = 0x1 | 0x40
                    flag2 = 0x1 | 0x80
                    if p1 < 0: flag1 |= 0x4
                    if p2 < 0: flag2 |= 0x4
                    if st1 == "-": flag1 |= 0x10
                    if st2 == "-": flag2 |= 0x10
                    p1_out = p1 + 1 if p1 >= 0 else 0
                    p2_out = p2 + 1 if p2 >= 0 else 0
                    cigar1 = cg1 if p1 >= 0 else "*"
                    cigar2 = cg2 if p2 >= 0 else "*"
                    mapq1 = mq1 if p1 >= 0 else 0
                    mapq2 = mq2 if p2 >= 0 else 0
                    rnext1 = "=" if p2 >= 0 else "*"
                    rnext2 = "=" if p1 >= 0 else "*"
                    sam.write(f"{q1}\t{flag1}\t{self.ref_names[0]}\t{p1_out}\t{mapq1}\t{cigar1}\t"
                              f"{rnext1}\t{p2_out}\t0\t{s1}\t{ql1}\n")
                    sam.write(f"{q2}\t{flag2}\t{self.ref_names[0]}\t{p2_out}\t{mapq2}\t{cigar2}\t"
                              f"{rnext2}\t{p1_out}\t0\t{s2}\t{ql2}\n")
                if progress_callback:
                    percent = int((total_processed / total_est_reads) * 100)
                    progress_callback(percent,
                        f"Aligning {total_processed}/{total_est_reads} paired reads...")
        if progress_callback:
            progress_callback(100,
                f"Finalizing {total_processed} paired reads...")
        self.convert_sam_to_bam(output_sam)
        return True

    def encode_seq(self, seq: str) -> bytes:
        out = bytearray()
        seq = seq.upper()
        i = 0
        while i < len(seq):
            b1 = self.nt_to_code.get(seq[i], 15)
            if i + 1 < len(seq):
                b2 = self.nt_to_code.get(seq[i + 1], 15)
            else:
                b2 = 0
            packed = (b1 << 4) | (b2 & 0xF)
            out.append(packed)
            i += 2
        return bytes(out)

    def encode_cigar(self, cigar_str: str) -> bytes:
        if cigar_str == "*" or cigar_str == "":
            return b""
        ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar_str)
        out = bytearray()
        for length, op in ops:
            length = int(length)
            op_int = self.op_to_int.get(op)
            if op_int is None:
                raise ValueError(f"Invalid CIGAR op: {op}")
            encoded = (length << 4) | (op_int & 0xF)
            out.extend(struct.pack("<I", encoded))
        return bytes(out)

    def reg2bin(self, beg: int, end: int) -> int:
        end -= 1
        if beg >> 14 == end >> 14:
            return ((1 << 15) - 1) // 7 + (beg >> 14)
        if beg >> 17 == end >> 17:
            return ((1 << 12) - 1) // 7 + (beg >> 17)
        if beg >> 20 == end >> 20:
            return ((1 << 9) - 1) // 7 + (beg >> 20)
        if beg >> 23 == end >> 23:
            return ((1 << 6) - 1) // 7 + (beg >> 23)
        if beg >> 26 == end >> 26:
            return ((1 << 3) - 1) // 7 + (beg >> 26)
        return 0

    def compute_nm_md(self, read_seq, ref_subseq, cigar):
        md_parts = []
        nm = 0
        run_match = 0
        read_i = 0
        ref_i = 0
        if not cigar or cigar == "*":
            for a, b in zip(read_seq, ref_subseq):
                if a == b:
                    run_match += 1
                else:
                    if run_match:
                        md_parts.append(str(run_match))
                    md_parts.append(b)
                    run_match = 0
                    nm += 1
            if run_match:
                md_parts.append(str(run_match))
            md = "".join(md_parts)
            return nm, md
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        for length_s, op in ops:
            length = int(length_s)
            if op in ("M", "=", "X"):
                for k in range(length):
                    rbase = read_seq[read_i] if read_i < len(read_seq) else 'N'
                    refbase = ref_subseq[ref_i] if ref_i < len(ref_subseq) else 'N'
                    if rbase == refbase:
                        run_match += 1
                    else:
                        if run_match:
                            md_parts.append(str(run_match))
                        md_parts.append(refbase)
                        run_match = 0
                        nm += 1
                    read_i += 1
                    ref_i += 1
            elif op == "I":
                nm += length
                read_i += length
            elif op == "D":
                if run_match:
                    md_parts.append(str(run_match))
                    run_match = 0
                del_seq = ref_subseq[ref_i:ref_i + length]
                md_parts.append("^" + del_seq)
                nm += length
                ref_i += length
            elif op in ("S", "H"):
                if op == "S":
                    read_i += length
            elif op == "N":
                ref_i += length
        if run_match:
            md_parts.append(str(run_match))
        md = "".join(md_parts)
        return nm, md

    def convert_sam_to_bam(self, sam_file, reference_fasta=None):
        ref_dict = self.ref_dict
        if reference_fasta:
            ref_dict = self._load_reference_dict(reference_fasta)
        bam_file = sam_file.replace(".sam", ".bam")
        with open(sam_file, "r") as sam, open(bam_file, "wb") as bam:
            bam.write(b"BAM\x01")
            header_lines = []
            ref_names = []
            ref_lengths = []
            alignments = []
            for line in sam:
                if line.startswith("@"):
                    header_lines.append(line)
                    if line.startswith("@SQ"):
                        parts = dict(f.split(":", 1) for f in line.strip().split("\t")[1:])
                        ref_names.append(parts.get("SN"))
                        ref_lengths.append(int(parts.get("LN", 0)))
                else:
                    if line.strip():
                        alignments.append(line.rstrip("\n"))

            header_text = "".join(header_lines)
            header_bytes = header_text.encode()
            bam.write(struct.pack("<i", len(header_bytes)))
            bam.write(header_bytes)
            bam.write(struct.pack("<i", len(ref_names)))
            for name, length in zip(ref_names, ref_lengths):
                name_b = name.encode() + b"\x00"
                bam.write(struct.pack("<i", len(name_b)))
                bam.write(name_b)
                bam.write(struct.pack("<i", length))
            for line in alignments:
                fields = line.split("\t")
                if len(fields) < 11:
                    continue
                qname, flag_s, rname, pos_s, mapq_s, cigar, rnext, pnext_s, tlen_s, seq, qual = fields[:11]
                try:
                    flag = int(flag_s)
                except:
                    flag = 0
                pos = int(pos_s) - 1 if pos_s.lstrip("-").isdigit() else 0
                mapq = int(mapq_s) if mapq_s.isdigit() else 255
                tid = ref_names.index(rname) if rname in ref_names else -1
                if rnext == "*":
                    rnext_tid = -1
                elif rnext == "=":
                    rnext_tid = tid
                else:
                    rnext_tid = ref_names.index(rnext) if rnext in ref_names else -1
                pnext = int(pnext_s) if pnext_s.lstrip("-").isdigit() else 0
                tlen = int(tlen_s) if tlen_s.lstrip("-").isdigit() else 0
                ref_seq = ref_dict.get(rname, "")
                ref_consumed = 0
                for length_s, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
                    l = int(length_s)
                    if op in ("M", "=", "X", "D", "N"):
                        ref_consumed += l
                ref_seq_sub = ""
                if ref_seq and pos >= 0:
                    ref_seq_sub = ref_seq[pos:pos + ref_consumed]
                cigar_bytes = self.encode_cigar(cigar)
                n_cigar_op = len(cigar_bytes) // 4
                seq_bytes = self.encode_seq(seq)
                if qual == "*" or qual == "":
                    qual_bytes = b"\xff" * len(seq)
                else:
                    qual_bytes = bytes([max(0, ord(c) - 33) for c in qual])
                l_read_name = len(qname) + 1
                bin_field = self.reg2bin(max(0, pos), max(0, pos + max(1, ref_consumed)))
                nm, md = self.compute_nm_md(seq, ref_seq_sub, cigar)
                md_val = md.encode()
                optional_bytes = b""
                optional_bytes += b"NM" + b"i" + struct.pack("<i", nm)
                optional_bytes += b"MD" + b"Z" + md_val
                block_size = 32 + l_read_name + len(cigar_bytes) + len(seq_bytes) + len(qual_bytes) + len(optional_bytes)
                bam.write(struct.pack("<i", block_size))
                bam.write(struct.pack("<i", tid))
                bam.write(struct.pack("<i", pos))
                bam.write(struct.pack("<I", (bin_field << 16) | (mapq << 8) | (l_read_name & 0xff)))
                bam.write(struct.pack("<I", (flag << 16) | (n_cigar_op & 0xffff)))
                bam.write(struct.pack("<i", len(seq)))
                bam.write(struct.pack("<i", rnext_tid))
                bam.write(struct.pack("<i", pnext))
                bam.write(struct.pack("<i", tlen))
                bam.write(qname.encode() + b"\x00")
                bam.write(cigar_bytes)
                bam.write(seq_bytes)
                bam.write(qual_bytes)
                bam.write(optional_bytes)
        return bam_file
    