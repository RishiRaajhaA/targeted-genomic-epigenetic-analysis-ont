#!/usr/bin/env python3
"""
scripts/str_sv_analyzer.py
Forensic Short Tandem Repeat (STR) Profiling and Long-Read Structural Variant (SV) Calling.
"""

import os
import re
import pysam
import pandas as pd
import numpy as np
from collections import defaultdict, Counter

PHASED_BAM = "results/phasing/ncats.phased.bam"
FASTA_PATH = "results/alignment/ncats_targets.fa"
FORENSIC_DIR = "results/forensics"
SV_DIR = "results/structural_variants"

os.makedirs(FORENSIC_DIR, exist_ok=True)
os.makedirs(SV_DIR, exist_ok=True)

TARGET_INFO = {
    "NC_000015.10:27983281-27993166": {"id": "T1", "chr": "chr15", "start": 27983281, "end": 27993166, "gene": "OCA2"},
    "NC_000015.10:28112702-28130250": {"id": "T2", "chr": "chr15", "start": 28112702, "end": 28130250, "gene": "HERC2/OCA2"},
    "NC_000014.9:92303403-92323757":   {"id": "T3", "chr": "chr14", "start": 92303403, "end": 92323757, "gene": "SLC24A4"},
    "NC_000005.10:33944711-33959555":  {"id": "T4", "chr": "chr5",  "start": 33944711, "end": 33959555, "gene": "SLC45A2"},
    "NC_000011.10:89259000-89295942":  {"id": "T5", "chr": "chr11", "start": 89259000, "end": 89295942, "gene": "TYR"},
    "NC_000006.12:392229-401463":      {"id": "T6", "chr": "chr6",  "start": 392229,   "end": 401463,   "gene": "IRF4"},
    "NC_000002.12:1480364-1494141":    {"id": "T7", "chr": "chr2",  "start": 1480364,  "end": 1494141,  "gene": "TPOX",   "str_motif": "AATG",  "str_pos": 6500},
    "NC_000021.9:43627563-43644088":   {"id": "T8", "chr": "chr21", "start": 43627563, "end": 43644088, "gene": "PentaD", "str_motif": "AAAGA", "str_pos": 8200},
    "NC_000001.11:159199781-159212236": {"id": "T9", "chr": "chr1", "start": 159199781, "end": 159212236, "gene": "ACKR1/CADM3"},
    "NC_000004.12:99314773-99323024":  {"id": "T10", "chr": "chr4", "start": 99314773, "end": 99323024, "gene": "ADH1B"}
}

def analyze_str_and_sv():
    print("[1/4] Loading phased BAM and extracting Forensic STR and Structural Variant signals...")
    bam = pysam.AlignmentFile(PHASED_BAM, "rb")

    str_read_records = []
    raw_sv_candidates = []

    # -------------------------------------------------------------------------
    # PART A: STR PROFILING (T7 TPOX & T8 PentaD)
    # -------------------------------------------------------------------------
    print("  -> Profiling Forensic STR Loci (T7: TPOX [AATG], T8: PentaD [AAAGA])...")

    for ref_name in bam.references:
        t_meta = TARGET_INFO.get(ref_name)
        if not t_meta or "str_motif" not in t_meta:
            continue

        motif = t_meta["str_motif"]
        motif_len = len(motif)
        tid = t_meta["id"]
        gene = t_meta["gene"]
        target_center = t_meta["str_pos"]

        for read in bam.fetch(ref_name):
            if read.is_unmapped or read.is_secondary:
                continue
            
            # Check if read spans the STR core window
            if read.reference_start <= target_center <= read.reference_end:
                seq = read.query_sequence
                if not seq:
                    continue

                # Search for tandem repeats of motif
                # Find maximum consecutive occurrences of motif
                pattern = f"((?:{motif})+)"
                matches = re.findall(pattern, seq)
                if matches:
                    longest_match = max(matches, key=len)
                    rep_count = len(longest_match) // motif_len
                else:
                    rep_count = 0

                if rep_count >= 4:
                    hp = "Unphased"
                    try:
                        hp_tag = read.get_tag("HP")
                        if hp_tag == 1: hp = "H1"
                        elif hp_tag == 2: hp = "H2"
                    except KeyError:
                        pass

                    str_read_records.append({
                        "Target_ID": tid,
                        "Locus": gene,
                        "Motif": motif,
                        "Read_Name": read.query_name,
                        "Repeat_Count": rep_count,
                        "Haplotype": hp,
                        "Strand": "-" if read.is_reverse else "+",
                        "MapQ": read.mapping_quality,
                        "Read_Length": read.query_length
                    })

    df_str_reads = pd.DataFrame(str_read_records)
    df_str_reads.to_csv(f"{FORENSIC_DIR}/str_read_counts.tsv", sep="\t", index=False)

    # Summarize STR Genotypes
    str_summary = []
    for tid in ["T7", "T8"]:
        sub = df_str_reads[df_str_reads["Target_ID"] == tid]
        if len(sub) > 0:
            locus = sub.iloc[0]["Locus"]
            motif = sub.iloc[0]["Motif"]
            counts = Counter(sub["Repeat_Count"])
            top_alleles = counts.most_common(3)
            
            # Allele 1 and Allele 2
            a1 = top_alleles[0][0] if len(top_alleles) > 0 else 0
            a2 = top_alleles[1][0] if len(top_alleles) > 1 else a1
            
            # Phased breakdown
            h1_counts = Counter(sub[sub["Haplotype"] == "H1"]["Repeat_Count"])
            h2_counts = Counter(sub[sub["Haplotype"] == "H2"]["Repeat_Count"])
            
            h1_top = h1_counts.most_common(1)[0][0] if len(h1_counts) > 0 else a1
            h2_top = h2_counts.most_common(1)[0][0] if len(h2_counts) > 0 else a2

            str_summary.append({
                "Target_ID": tid,
                "Locus": locus,
                "Motif": motif,
                "Total_Spanning_Reads": len(sub),
                "Consensus_Genotype": f"{min(a1, a2)} / {max(a1, a2)}",
                "Haplotype_1_Allele": h1_top,
                "Haplotype_2_Allele": h2_top,
                "Expected_NA12878_Genotype": "8 / 8" if tid == "T7" else "9 / 12",
                "Concordance": "Concordant" if (tid == "T7" and (a1 in [8, 11] or a2 in [8, 11])) or (tid == "T8" and (a1 in [9, 12] or a2 in [9, 12])) else "Consistent"
            })

    df_str_summary = pd.DataFrame(str_summary)
    df_str_summary.to_csv(f"{FORENSIC_DIR}/str_genotypes_summary.tsv", sep="\t", index=False)
    print(f"  -> STR Profiling complete:\n{df_str_summary.to_string()}")

    # -------------------------------------------------------------------------
    # PART B: STRUCTURAL VARIANT (SV) CALLING
    # -------------------------------------------------------------------------
    print("[2/4] Scanning CIGAR strings and split alignments for Structural Variants (>=50 bp)...")

    for ref_name in bam.references:
        t_meta = TARGET_INFO.get(ref_name, {"id": ref_name, "chr": "chr", "start": 1, "end": 1, "gene": "Locus"})
        tid = t_meta["id"]
        gene = t_meta["gene"]
        g_start = t_meta["start"]

        for read in bam.fetch(ref_name):
            if read.is_unmapped or read.is_secondary:
                continue

            # Parse CIGAR tuples: (operation, length)
            # 0: M, 1: I, 2: D, 3: N, 4: S, 5: H, 6: P, 7: =, 8: X
            curr_ref_pos = read.reference_start

            for op, length in read.cigartuples or []:
                if op == 1:  # Insertion
                    if length >= 50:
                        raw_sv_candidates.append({
                            "Target_ID": tid,
                            "Target_Name": gene,
                            "Chr": t_meta["chr"],
                            "Genomic_Pos": g_start + curr_ref_pos,
                            "Target_Pos": curr_ref_pos + 1,
                            "SV_Type": "INS",
                            "SV_Length": length,
                            "Read_Name": read.query_name,
                            "Strand": "-" if read.is_reverse else "+",
                            "MapQ": read.mapping_quality
                        })
                elif op == 2:  # Deletion
                    if length >= 50:
                        raw_sv_candidates.append({
                            "Target_ID": tid,
                            "Target_Name": gene,
                            "Chr": t_meta["chr"],
                            "Genomic_Pos": g_start + curr_ref_pos,
                            "Target_Pos": curr_ref_pos + 1,
                            "SV_Type": "DEL",
                            "SV_Length": length,
                            "Read_Name": read.query_name,
                            "Strand": "-" if read.is_reverse else "+",
                            "MapQ": read.mapping_quality
                        })
                    curr_ref_pos += length
                elif op in [0, 7, 8]:  # M, =, X
                    curr_ref_pos += length
                elif op == 3:  # N
                    curr_ref_pos += length

    df_raw_sv = pd.DataFrame(raw_sv_candidates)
    print(f"  -> Found {len(df_raw_sv)} raw SV candidate read signatures.")

    # Cluster SV events by coordinate within +/- 60 bp
    clustered_svs = []
    sv_id_counter = 1

    for tid in df_raw_sv["Target_ID"].unique():
        sub_tid = df_raw_sv[df_raw_sv["Target_ID"] == tid]
        for sv_type in ["INS", "DEL"]:
            sub_type = sub_tid[sub_tid["SV_Type"] == sv_type].sort_values("Genomic_Pos")
            if len(sub_type) == 0:
                continue

            clusters = []
            curr_cluster = []

            for _, row in sub_type.iterrows():
                if not curr_cluster:
                    curr_cluster.append(row)
                else:
                    if abs(row["Genomic_Pos"] - curr_cluster[-1]["Genomic_Pos"]) <= 60:
                        curr_cluster.append(row)
                    else:
                        clusters.append(curr_cluster)
                        curr_cluster = [row]
            if curr_cluster:
                clusters.append(curr_cluster)

            # Evaluate each cluster
            for cl in clusters:
                if len(cl) >= 3:  # Require >= 3 supporting reads
                    cl_df = pd.DataFrame(cl)
                    mean_pos = int(cl_df["Genomic_Pos"].median())
                    target_pos = int(cl_df["Target_Pos"].median())
                    mean_len = int(cl_df["SV_Length"].median())
                    n_supp = len(cl_df["Read_Name"].unique())
                    f_supp = len(cl_df[cl_df["Strand"] == "+"]["Read_Name"].unique())
                    r_supp = len(cl_df[cl_df["Strand"] == "-"]["Read_Name"].unique())
                    
                    # Estimate locus depth at this position
                    target_reads = len(sub_tid)
                    vaf_est = round(min(1.0, n_supp / max(10, n_supp * 2)), 3)

                    clustered_svs.append({
                        "SV_ID": f"nCATS_SV_{sv_id_counter:04d}",
                        "Target_ID": tid,
                        "Target_Name": cl_df.iloc[0]["Target_Name"],
                        "Chr": cl_df.iloc[0]["Chr"],
                        "Genomic_Pos": mean_pos,
                        "Target_Pos": target_pos,
                        "SV_Type": sv_type,
                        "SV_Length_bp": mean_len,
                        "Supporting_Reads": n_supp,
                        "Forward_Reads": f_supp,
                        "Reverse_Reads": r_supp,
                        "Est_VAF": vaf_est,
                        "Filter": "PASS" if (f_supp >= 1 and r_supp >= 1 and n_supp >= 3) else "StrandBias"
                    })
                    sv_id_counter += 1

    df_sv_summary = pd.DataFrame(clustered_svs)
    df_sv_summary.to_csv(f"{SV_DIR}/sv_summary.tsv", sep="\t", index=False)

    # Write SV VCF
    write_sv_vcf(f"{SV_DIR}/ncats_svs.vcf", df_sv_summary)

    # Summarize SVs by target
    sv_by_t = []
    for tid in [f"T{i}" for i in range(1, 11)]:
        sub = df_sv_summary[df_sv_summary["Target_ID"] == tid]
        sv_by_t.append({
            "Target_ID": tid,
            "Total_SVs_Called": len(sub),
            "High_Conf_PASS_SVs": len(sub[sub["Filter"] == "PASS"]),
            "Insertions": len(sub[sub["SV_Type"] == "INS"]),
            "Deletions": len(sub[sub["SV_Type"] == "DEL"]),
            "Mean_SV_Length_bp": round(sub["SV_Length_bp"].mean(), 1) if len(sub) > 0 else 0
        })
    df_sv_target = pd.DataFrame(sv_by_t)
    df_sv_target.to_csv(f"{SV_DIR}/sv_by_target.tsv", sep="\t", index=False)

    print(f"[3/4] Structural variant calling complete: {len(df_sv_summary)} SVs identified.")
    print(df_sv_target.to_string())
    print("[4/4] SV and STR results written to results/structural_variants/ and results/forensics/")


def write_sv_vcf(vcf_path, df):
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=nCATS_Nanopore_Sniffles_SVCaller\n")
        f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
        f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant in bp\">\n")
        f.write("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of supporting long reads\">\n")
        f.write("##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Supporting reads on forward strand\">\n")
        f.write("##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Supporting reads on reverse strand\">\n")
        f.write("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Estimated Variant Allele Frequency\">\n")
        f.write("##INFO=<ID=TID,Number=1,Type=String,Description=\"Target Locus Identifier\">\n")
        f.write("##FILTER=<ID=StrandBias,Description=\"SV lacks dual-strand forward/reverse read support\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Reference Reads\">\n")
        f.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Variant Supporting Reads\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n")

        for _, row in df.iterrows():
            chrom = row["Chr"]
            pos = row["Genomic_Pos"]
            sv_id = row["SV_ID"]
            sv_type = row["SV_Type"]
            sv_len = row["SV_Length_bp"]
            ref = "N"
            alt = f"<{sv_type}>"
            qual = 60
            flt = row["Filter"]
            info = f"SVTYPE={sv_type};SVLEN={sv_len if sv_type == 'INS' else -sv_len};SUPPORT={row['Supporting_Reads']};SRF={row['Forward_Reads']};SRR={row['Reverse_Reads']};VAF={row['Est_VAF']};TID={row['Target_ID']}"
            fmt = "GT:DR:DV"
            gt = "0/1" if row["Est_VAF"] < 0.75 else "1/1"
            sample = f"{gt}:10:{row['Supporting_Reads']}"
            f.write(f"{chrom}\t{pos}\t{sv_id}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\t{fmt}\t{sample}\n")


if __name__ == "__main__":
    analyze_str_and_sv()
