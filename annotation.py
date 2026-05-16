# Annotation Backend Code
import pandas as pd
from collections import defaultdict
# VEP-like Annotation Engine
class VEPAnnotation:
    CODON_TABLE = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
    }
    AA3 = {
        "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys",
        "Q":"Gln","E":"Glu","G":"Gly","H":"His","I":"Ile",
        "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro",
        "S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val",
        "*":"Ter"
    }
    @staticmethod
    def codon_to_aa(codon):
        return VEPAnnotation.CODON_TABLE.get(codon.upper(), "?")
    @staticmethod
    def parse_gff(gff_file):
        genes = defaultdict(list)
        with open(gff_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip().split("\t")
                if len(cols) < 9:
                    continue
                chrom, _, feature, start, end, _, strand, _, attrs = cols
                if feature not in ("gene", "CDS"):
                    continue
                attr = {}
                for a in attrs.split(";"):
                    if "=" in a:
                        k, v = a.split("=", 1)
                        attr[k] = v
                gene = (
                    attr.get("gene")
                    or attr.get("gene_name")
                    or attr.get("locus_tag")
                    or "-"
                )
                genes[chrom].append({
                    "type": feature,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene": gene
                })
        return genes
    @staticmethod
    def annotate_variants(variants, gff_file, fasta):
        genes = VEPAnnotation.parse_gff(gff_file)
        with open(fasta) as f:
            ref_seq = "".join(l.strip().upper() for l in f if not l.startswith(">"))
        results = []
        if isinstance(variants, pd.DataFrame):
            variants = variants.values.tolist()
        for v in variants:
            if not isinstance(v, (list, tuple)) or len(v) != 8:
                continue
            chrom, pos, end, ref, alt, dp, refc, altc = v
            gene = "-"
            region = "Intergenic"
            consequence = "intergenic_variant"
            impact = "LOW"
            hgvsc = f"c.{pos}{ref}>{alt}"
            hgvsp = "-"
            aa_change = "-"
            if len(ref) == len(alt) == 1:
                variant_name = f"{chrom}:g.{pos}{ref}>{alt}"
            else:
                variant_name = f"{chrom}:g.{pos}_{end}delins{ref}>{alt}"
            if chrom in genes:
                for g in genes[chrom]:
                    if g["start"] <= pos <= g["end"]:
                        gene = g["gene"]
                        if g["type"] == "CDS":
                            region = "CDS"
                            if len(ref) != len(alt):
                                consequence = "frameshift_variant"
                                impact = "HIGH"
                                codon_start = pos - ((pos - 1) % 3)
                                ref_codon = ref_seq[codon_start-1:codon_start+2]
                                ref_aa = VEPAnnotation.codon_to_aa(ref_codon)
                                prot_pos = ((pos - 1) // 3) + 1
                                hgvsp = f"p.{VEPAnnotation.AA3.get(ref_aa, '?')}{prot_pos}fs"
                                aa_change = f"{ref_aa}/-" if ref_aa != "*" else "-/-"
                            else:
                                codon = ref_seq[pos-1:pos+2]
                                if len(codon) == 3:
                                    alt_codon = alt + codon[1:]
                                    ref_aa = VEPAnnotation.codon_to_aa(codon)
                                    alt_aa = VEPAnnotation.codon_to_aa(alt_codon)
                                    aa_change = f"{ref_aa}>{alt_aa}"
                                    hgvsp = f"p.{VEPAnnotation.AA3.get(ref_aa)}{(pos//3)+1}{VEPAnnotation.AA3.get(alt_aa)}"
                                    if ref_aa == alt_aa:
                                        consequence = "synonymous_variant"
                                        impact = "LOW"
                                    elif alt_aa == "*":
                                        consequence = "nonsense_variant"
                                        impact = "HIGH"
                                    else:
                                        consequence = "missense_variant"
                                        impact = "MODERATE"
                        break
            results.append([
                chrom, pos, end,
                ref, alt,
                dp, refc, altc,
                consequence, impact,
                gene, region,
                hgvsc, hgvsp, aa_change,
                variant_name 
            ])
        return results
def parse_vcf(vcf_path):
    variants = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4].split(",")[0]
            info = cols[7]
            end = pos + len(ref) - 1
            dp = refc = altc = 0
            for field in info.split(";"):
                if field.startswith("DP="):
                    dp = int(field.split("=")[1])
                elif field.startswith("AC="):
                    ac_val = field.split("=")[1]                   
                    # handle multi-allelic case like "2,1"
                    if "," in ac_val:
                        altc = int(ac_val.split(",")[0])   # take first ALT only
                    else:
                        altc = int(ac_val)
            if dp > 0:
                refc = dp - altc
            variants.append((chrom, pos, end, ref, alt, dp, refc, altc))
    return variants
# WRITE ANNOTATED VCF
def write_annotated_vcf(input_vcf, annotated_df, output_vcf):
    ann_map = {
        (r["Chromosome"], int(r["Start"]), r["Ref"], r["Alt"]): r
        for _, r in annotated_df.iterrows()
    }
    with open(input_vcf) as fin, open(output_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            cols = line.strip().split("\t")
            key = (cols[0], int(cols[1]), cols[3], cols[4].split(",")[0])
            if key in ann_map:
                r = ann_map[key]
                ann = (
                    f"{r['Consequence']}|{r['Impact']}|{r['Gene']}|"
                    f"{r['Region']}|{r['HGVS.c']}|{r['HGVS.p']}|{r['VariantName']}"
                )
                cols[7] += f";ANN={ann}"
            fout.write("\t".join(cols) + "\n")
# FRONTEND SAFE API
def annotate(vcf_path, gtf=None, gff=None, fasta=None, output_vcf="annotated.vcf"):
    annotation_file = gtf if gtf else gff
    if annotation_file is None:
        raise ValueError("GTF or GFF required")
    if fasta is None:
        raise ValueError("FASTA required")
    variants = parse_vcf(vcf_path)
    rows = VEPAnnotation.annotate_variants(variants, annotation_file, fasta)
    df = pd.DataFrame(
        rows,
        columns=[
            "Chromosome","Start","End",
            "Ref","Alt",
            "Depth","RefCount","AltCount",
            "Consequence","Impact",
            "Gene","Region",
            "HGVS.c","HGVS.p","AA_Change",
            "VariantName"  
        ]
    )
    # ✅ STORE ANNOTATED VCF
    write_annotated_vcf(vcf_path, df, output_vcf)
    import streamlit as st
    st.session_state["annotation_df"] = df  
    return df
