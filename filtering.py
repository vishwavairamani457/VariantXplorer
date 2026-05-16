# Filtering Backend Code
import pandas as pd
import streamlit as st
class VariantFilter:
    def __init__(self, annotated_vcf_path, filter_params=None, progress_callback=None):
        self.vcf_path = annotated_vcf_path
        self.progress_callback = progress_callback
        # Default parameters
        self.filter_params = self.default_params()
        # Override with frontend values
        if filter_params and isinstance(filter_params, dict):
            for key, value in filter_params.items():
                if key in self.filter_params:
                    self.filter_params[key] = value
        # Normalize Allowed Impacts
        allowed_impacts = self.filter_params.get("Allowed Impacts", [])
        if isinstance(allowed_impacts, str):
            try:
                allowed_impacts = eval(allowed_impacts)
            except:
                allowed_impacts = []
        self.filter_params["Allowed Impacts"] = [
            str(x).strip().upper() for x in allowed_impacts
        ]
        # Normalize Variant Types
        allowed_types = self.filter_params.get("Allowed Variant Types", [])
        if isinstance(allowed_types, str):
            try:
                allowed_types = eval(allowed_types)
            except:
                allowed_types = []
        self.filter_params["Allowed Variant Types"] = [
            str(x).strip().lower() for x in allowed_types
        ]
        self.variants = []
    # Default Parameters
    @staticmethod
    def default_params():
        return {
            "Min Depth": 0,
            "Max Allele Frequency": 1.0,
            "Allowed Impacts": ["HIGH", "MODERATE", "LOW", "MODIFIER"],
            "Allowed Variant Types": [
                "missense_variant",
                "synonymous_variant",
                "stop_gained",
                "frameshift_variant",
                "nonsense_variant", 
                "intergenic_variant",
                "intron_variant",
                "splice_donor_variant",
                "splice_acceptor_variant",
                "splice_region_variant"
            ]
        }
    # Load VCF
    def load_vcf(self):
        df = pd.read_csv(
            self.vcf_path,
            sep="\t",
            engine="python",
            comment="#",
            header=None
        )
        df = df.iloc[:, :8]
        df.columns = [
            "CHROM", "POS", "ID", "REF", "ALT",
            "QUAL", "FILTER", "INFO"
        ]
        # Parse INFO field
        def parse_info(info):
            info_dict = {}
            for entry in str(info).split(";"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    info_dict[k] = v
            return info_dict
        info_parsed = df["INFO"].apply(parse_info)
        # Depth
        df["DP"] = info_parsed.apply(lambda x: int(x.get("DP", 0)))
        # Allele Frequency
        def parse_af(x):
            af_val = x.get("AF", "0.0")
            if isinstance(af_val, str) and "," in af_val:
                return float(af_val.split(",")[0])
            try:
                return float(af_val)
            except:
                return 0.0
        df["AF"] = info_parsed.apply(parse_af)
        # FIXED: Extract Annotation + Impact
        def extract_annotation(info):
            ann = info.get("ANN")
            if not ann:
                return "unknown", "LOW"
            try:
                first_ann = ann.split(",")[0]
                fields = first_ann.split("|")
                # Detect format automatically
                if fields[0] in ["A","T","G","C"]:  
                    # VEP format
                    annotation = fields[1].lower()
                    impact = fields[2].upper()
                else:
                    # Your custom format
                    annotation = fields[0].lower()
                    impact = fields[1].upper()
                return annotation, impact
            except:
                return "unknown", "LOW"
        df["Annotation"], df["Impact"] = zip(*info_parsed.apply(extract_annotation))

        # Region classification (optional but useful)
        def classify_region(annotation):
            if "splice_donor" in annotation or "splice_acceptor" in annotation:
                return "Splice_Site"
            elif "splice_region" in annotation:
                return "Splice_Region"
            elif "intron" in annotation:
                return "Intronic"
            elif "missense" in annotation or "synonymous" in annotation:
                return "Exonic"
            else:
                return "Other"
        df["Region"] = df["Annotation"].apply(classify_region)
        df["FILTER"] = "PASS"
        self.variants = df.to_dict(orient="records")
        unique_types = set([v["Annotation"] for v in self.variants])
        # print("Available Variant Types:", unique_types)
    # Apply Filters
    def pass_filters(self):
        passed = []
        total = len(self.variants)
        p = self.filter_params
        for i, var in enumerate(self.variants):
            dp = int(var.get("DP", 0))
            af = float(var.get("AF", 0.0))
            impact = str(var.get("Impact", "LOW")).strip().upper()
            variant_type = str(var.get("Annotation", "")).lower()
            # DP filter
            if dp < p["Min Depth"]:
                continue
            # AF filter (FIXED)
            if af > p["Max Allele Frequency"]:
                continue
            # Impact filter
            if impact not in p["Allowed Impacts"]:
                continue
            # Variant type filter (FIXED - partial match)
            allowed_types = p["Allowed Variant Types"]
            if not any(t in variant_type for t in allowed_types):
                continue
            passed.append(var)
            if self.progress_callback:
                self.progress_callback(i + 1, total)
            # print(f"DP={dp}, AF={af}, Impact={impact}, Type={variant_type}")
        return pd.DataFrame(passed)

    # Export
    @staticmethod
    def export_to_tsv(variants, output_path="filtered_variants.tsv"):
        rows = []
        for v in variants:
            rows.append({
                "CHROM": v.get("CHROM"),
                "POS": v.get("POS"),
                "REF": v.get("REF"),
                "ALT": v.get("ALT"),
                "FILTER": v.get("FILTER", "PASS"),
                "Annotation": v.get("Annotation"),
                "Region": v.get("Region"),
                "Impact": v.get("Impact"),
                "DP": v.get("DP"),
                "AF": v.get("AF")
            })
        df = pd.DataFrame(rows)
        df.to_csv(output_path, sep="\t", index=False)
        st.write(f"✅ Exported filtered variants to {output_path}")
        st.session_state["filtering_df"] = df
        return output_path