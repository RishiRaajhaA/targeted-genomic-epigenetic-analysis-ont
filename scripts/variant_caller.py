#!/usr/bin/env python3
"""
scripts/variant_caller.py
High-Accuracy Long-Read Variant Calling & Dual-Strand Filtering for nCATS Nanopore Data.
"""

import os
import sys
import math
import pysam
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

BAM_PATH = "results/alignment/ncats.sorted.bam"
FASTA_PATH = "results/alignment/ncats_targets.fa"
BED_PATH = "results/alignment/targets.bed"
OUT_DIR = "results/variants"

os.makedirs(OUT_DIR, exist_ok=True)

# Target annotation metadata (from targets.bed and GRCh38)
TARGET_INFO = {
    "NC_000015.10:27983281-27993166": {"id": "T1", "chr": "chr15", "start": 27983281, "end": 27993166, "marker": "rs1800407", "gene": "OCA2"},
    "NC_000015.10:28112702-28130250": {"id": "T2", "chr": "chr15", "start": 28112702, "end": 28130250, "marker": "rs12913832", "gene": "HERC2/OCA2"},
    "NC_000014.9:92303403-92323757":   {"id": "T3", "chr": "chr14", "start": 92303403, "end": 92323757, "marker": "rs12896399", "gene": "SLC24A4"},
    "NC_000005.10:33944711-33959555":  {"id": "T4", "chr": "chr5",  "start": 33944711, "end": 33959555, "marker": "rs16891982", "gene": "SLC45A2"},
    "NC_000011.10:89259000-89295942":  {"id": "T5", "chr": "chr11", "start": 89259000, "end": 89295942, "marker": "rs1393350",  "gene": "TYR"},
    "NC_000006.12:392229-401463":      {"id": "T6", "chr": "chr6",  "start": 392229,   "end": 401463,   "marker": "rs12203592", "gene": "IRF4"},
    "NC_000002.12:1480364-1494141":    {"id": "T7", "chr": "chr2",  "start": 1480364,  "end": 1494141,  "marker": "TPOX",       "gene": "TPOX"},
    "NC_000021.9:43627563-43644088":   {"id": "T8", "chr": "chr21", "start": 43627563, "end": 43644088, "marker": "PentaD",     "gene": "PentaD"},
    "NC_000001.11:159199781-159212236": {"id": "T9", "chr": "chr1", "start": 159199781, "end": 159212236, "marker": "rs2814778", "gene": "ACKR1/CADM3"},
    "NC_000004.12:99314773-99323024":  {"id": "T10", "chr": "chr4", "start": 99314773, "end": 99323024, "marker": "rs1229984",  "gene": "ADH1B"}
}

def call_variants():
    print("[1/4] Loading reference FASTA and indexed BAM...")
    fasta = pysam.FastaFile(FASTA_PATH)
    bam = pysam.AlignmentFile(BAM_PATH, "rb")

    raw_variants = []
    filtered_variants = []

    print("[2/4] Iterating pileups across all 10 target loci...")

    for ref_name in bam.references:
        t_meta = TARGET_INFO.get(ref_name, {"id": ref_name, "chr": ref_name, "start": 1, "end": 1, "marker": "", "gene": ""})
        ref_seq = fasta.fetch(ref_name).upper()
        ref_len = len(ref_seq)

        print(f"  -> Processing {t_meta['id']} ({t_meta['gene']}) | Length: {ref_len:,} bp...")

        for pileupcolumn in bam.pileup(ref_name, 0, ref_len, truncate=True, min_base_quality=10, min_mapping_quality=10):
            pos = pileupcolumn.pos  # 0-indexed
            if pos >= ref_len:
                continue
            ref_base = ref_seq[pos]
            if ref_base not in ["A", "C", "G", "T"]:
                continue

            # Forward and reverse base counts
            counts_f = {"A": 0, "C": 0, "G": 0, "T": 0, "DEL": 0}
            counts_r = {"A": 0, "C": 0, "G": 0, "T": 0, "DEL": 0}
            base_qualities = []

            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del:
                    if pileupread.alignment.is_reverse:
                        counts_r["DEL"] += 1
                    else:
                        counts_f["DEL"] += 1
                elif not pileupread.is_refskip:
                    b = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    q = pileupread.alignment.query_qualities[pileupread.query_position]
                    base_qualities.append(q)
                    if b in ["A", "C", "G", "T"]:
                        if pileupread.alignment.is_reverse:
                            counts_r[b] += 1
                        else:
                            counts_f[b] += 1

            total_dp_f = sum(counts_f.values())
            total_dp_r = sum(counts_r.values())
            dp = total_dp_f + total_dp_r

            if dp < 8:
                continue

            ref_f = counts_f.get(ref_base, 0)
            ref_r = counts_r.get(ref_base, 0)
            ref_dp = ref_f + ref_r

            # Check each possible alt base
            for alt_base in ["A", "C", "G", "T"]:
                if alt_base == ref_base:
                    continue

                alt_f = counts_f.get(alt_base, 0)
                alt_r = counts_r.get(alt_base, 0)
                alt_dp = alt_f + alt_r

                if alt_dp < 3:
                    continue

                vaf = alt_dp / dp
                if vaf < 0.15:
                    continue

                # Genotype classification
                if vaf >= 0.75:
                    gt = "1/1"
                elif vaf >= 0.20:
                    gt = "0/1"
                else:
                    gt = "0/1"

                # Dual strand check: require alt on both forward and reverse strands
                has_dual_strand = (alt_f >= 2 and alt_r >= 2)
                # Strand bias test (Fisher exact test on [[ref_f, ref_r], [alt_f, alt_r]])
                table = [[ref_f, ref_r], [alt_f, alt_r]]
                try:
                    _, pval = fisher_exact(table)
                except Exception:
                    pval = 1.0

                qual = round(min(99, -10 * math.log10(max(1e-10, pval if pval > 0 else 1e-10)) * (vaf * 10)))
                qual = max(20, qual)

                # Genomic coordinate (1-based GRCh38)
                genomic_pos = t_meta["start"] + pos
                filter_status = "PASS" if (has_dual_strand and vaf >= 0.18 and alt_dp >= 4) else "StrandBias"

                var_record = {
                    "Target_ID": t_meta["id"],
                    "Target_Name": t_meta["gene"],
                    "Chr": t_meta["chr"],
                    "Genomic_Pos": genomic_pos,
                    "Target_Pos": pos + 1,
                    "Ref": ref_base,
                    "Alt": alt_base,
                    "Type": "SNV",
                    "DP": dp,
                    "AD": alt_dp,
                    "AD_Forward": alt_f,
                    "AD_Reverse": alt_r,
                    "Ref_Forward": ref_f,
                    "Ref_Reverse": ref_r,
                    "VAF": round(vaf, 4),
                    "GT": gt,
                    "QUAL": qual,
                    "Fisher_P": round(pval, 6),
                    "Dual_Strand": has_dual_strand,
                    "Filter": filter_status
                }

                raw_variants.append(var_record)
                if filter_status == "PASS":
                    filtered_variants.append(var_record)

    df_raw = pd.DataFrame(raw_variants)
    df_pass = pd.DataFrame(filtered_variants)

    print(f"[3/4] Identified {len(df_raw)} raw candidate variants, {len(df_pass)} passed dual-strand concordance filter.")

    # Save TSV summaries
    df_raw.to_csv(f"{OUT_DIR}/variant_calls_raw_all.tsv", sep="\t", index=False)
    df_pass.to_csv(f"{OUT_DIR}/variant_calls_summary.tsv", sep="\t", index=False)

    # Write standard VCF files
    write_vcf(f"{OUT_DIR}/ncats_raw_variants.vcf", df_raw)
    write_vcf(f"{OUT_DIR}/ncats_filtered_variants.vcf", df_pass)

    # Annotated target markers check
    annotated = []
    for ref_name, t_meta in TARGET_INFO.items():
        marker = t_meta["marker"]
        gene = t_meta["gene"]
        # Find variants in this target
        t_vars = df_pass[df_pass["Target_ID"] == t_meta["id"]]
        annotated.append({
            "Target_ID": t_meta["id"],
            "Target_Locus": gene,
            "Associated_Marker": marker,
            "Chr_Region": f"{t_meta['chr']}:{t_meta['start']}-{t_meta['end']}",
            "High_Conf_SNVs_Found": len(t_vars),
            "Mean_VAF": round(t_vars["VAF"].mean(), 3) if len(t_vars) > 0 else 0.0,
            "Mean_Depth": round(t_vars["DP"].mean(), 1) if len(t_vars) > 0 else 0.0,
            "Heterozygous_SNVs": len(t_vars[t_vars["GT"] == "0/1"]),
            "Homozygous_SNVs": len(t_vars[t_vars["GT"] == "1/1"])
        })
    df_annot = pd.DataFrame(annotated)
    df_annot.to_csv(f"{OUT_DIR}/annotated_target_snps.tsv", sep="\t", index=False)

    print("[4/4] Variant calling complete. Files generated in results/variants/")
    print(df_annot.to_string())


def write_vcf(vcf_path, df):
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=nCATS_Nanopore_VariantCaller\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n")
        f.write("##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Alternate Allele Depth\">\n")
        f.write("##INFO=<ID=ADF,Number=1,Type=Integer,Description=\"Alternate Allele Depth on Forward Strand\">\n")
        f.write("##INFO=<ID=ADR,Number=1,Type=Integer,Description=\"Alternate Allele Depth on Reverse Strand\">\n")
        f.write("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant Allele Frequency\">\n")
        f.write("##INFO=<ID=TID,Number=1,Type=String,Description=\"Target Locus Identifier\">\n")
        f.write("##FILTER=<ID=StrandBias,Description=\"Variant fails dual-strand forward/reverse concordance filter\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        f.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allelic Depths\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n")

        for _, row in df.iterrows():
            chrom = row["Chr"]
            pos = row["Genomic_Pos"]
            var_id = f"{row['Target_ID']}_{row['Target_Pos']}"
            ref = row["Ref"]
            alt = row["Alt"]
            qual = row["QUAL"]
            flt = row["Filter"]
            info = f"DP={row['DP']};AD={row['AD']};ADF={row['AD_Forward']};ADR={row['AD_Reverse']};VAF={row['VAF']};TID={row['Target_ID']}"
            fmt = "GT:DP:AD"
            sample = f"{row['GT']}:{row['DP']}:{row['AD']}"
            f.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\t{fmt}\t{sample}\n")


if __name__ == "__main__":
    call_variants()
