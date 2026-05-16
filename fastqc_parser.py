
# import os
# import io
# import zipfile
# from typing import Dict, Any, Tuple, List

# def _open_fastqc_file(result_path: str) -> Tuple[str, io.TextIOBase]:
#     """
#     Return tuple (sample_name, handle) for fastqc_data.txt inside a FastQC result.
#     `result_path` may be either a directory produced by FastQC or the .zip file.
#     """
#     if os.path.isdir(result_path):
#         # Expect structure: result_path/fastqc_data.txt
#         data_fp = os.path.join(result_path, "fastqc_data.txt")
#         if not os.path.exists(data_fp):
#             # Sometimes folder name ends with "_fastqc", ensure we are inside it
#             alt = os.path.join(result_path, os.path.basename(result_path) + ".txt")
#             if os.path.exists(alt):
#                 data_fp = alt
#         if not os.path.exists(data_fp):
#             raise FileNotFoundError(f"fastqc_data.txt not found in {result_path}")
#         sample_name = os.path.basename(result_path).replace("_fastqc", "")
#         return sample_name, open(data_fp, "r", encoding="utf-8", errors="ignore")
#     else:
#         # Assume ZIP
#         sample_name = os.path.basename(result_path).replace("_fastqc.zip", "")
#         zf = zipfile.ZipFile(result_path, "r")
#         # Inside zip: <name>_fastqc/fastqc_data.txt
#         inner = None
#         for name in zf.namelist():
#             if name.endswith("fastqc_data.txt"):
#                 inner = name
#                 break
#         if inner is None:
#             raise FileNotFoundError("fastqc_data.txt not found inside zip")
#         return sample_name, io.TextIOWrapper(zf.open(inner, "r"), encoding="utf-8", errors="ignore")

# def parse_fastqc_data(result_path: str) -> Dict[str, Any]:
#     """
#     Parse selected sections from fastqc_data.txt.
#     Returns a dict with keys:
#       - 'basic_statistics', 'per_base_quality', 'per_sequence_gc_content',
#         'sequence_length_distribution', 'overrepresented_sequences', 'module_status'
#     """
#     sample, handle = _open_fastqc_file(result_path)
#     with handle:
#         lines = [ln.rstrip("\n") for ln in handle]

#     sections = {}
#     module_status = {}
#     i = 0
#     n = len(lines)

#     def read_table(start_idx: int) -> List[str]:
#         table = []
#         j = start_idx + 1  # skip header line
#         while j < n and not lines[j].startswith(">>END_MODULE"):
#             if not lines[j].startswith("#"):
#                 table.append(lines[j])
#             j += 1
#         return table

#     while i < n:
#         line = lines[i]
#         if line.startswith(">>"):
#             # e.g., >>Per base sequence quality\tpass
#             header = line[2:]
#             if "\t" in header:
#                 name, status = header.split("\t", 1)
#             else:
#                 name, status = header, "UNKNOWN"
#             module_status[name] = status.lower()

#             # Capture content until >>END_MODULE
#             table = read_table(i + 1)

#             if header.startswith("Basic Statistics"):
#                 stats = {}
#                 for row in table:
#                     if "\t" in row:
#                         k, v = row.split("\t", 1)
#                         stats[k] = v
#                 sections["basic_statistics"] = stats

#             elif header.startswith("Per base sequence quality"):
#                 # columns: Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile
#                 data = []
#                 for row in table:
#                     parts = row.split("\t")
#                     if len(parts) >= 2:
#                         data.append({
#                             "Base": parts[0],
#                             "Mean": float(parts[1])
#                         })
#                 sections["per_base_quality"] = data

#             elif header.startswith("Per sequence GC content"):
#                 # columns: GC Content\tCount
#                 data = []
#                 for row in table:
#                     parts = row.split("\t")
#                     if len(parts) >= 2:
#                         data.append({
#                             "GC": float(parts[0]),
#                             "Count": float(parts[1])
#                         })
#                 sections["per_sequence_gc_content"] = data

#             elif header.startswith("Sequence Length Distribution"):
#                 data = []
#                 for row in table:
#                     parts = row.split("\t")
#                     if len(parts) >= 2:
#                         data.append({
#                             "Length": parts[0],
#                             "Count": float(parts[1])
#                         })
#                 sections["sequence_length_distribution"] = data

#             elif header.startswith("Overrepresented sequences"):
#                 # columns: Sequence\tCount\tPercentage\tPossible Source
#                 data = []
#                 for row in table:
#                     parts = row.split("\t")
#                     if len(parts) >= 4:
#                         data.append({
#                             "Sequence": parts[0],
#                             "Count": float(parts[1]),
#                             "Percentage": parts[2],
#                             "Source": parts[3],
#                         })
#                 sections["overrepresented_sequences"] = data

#             # fast-forward to end
#             while i < n and not lines[i].startswith(">>END_MODULE"):
#                 i += 1
#         i += 1

#     sections["module_status"] = module_status
#     sections["sample"] = sample
#     return sections



