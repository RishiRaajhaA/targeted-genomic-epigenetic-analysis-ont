#!/usr/bin/env python3
"""
scripts/phasing.py
Long-Read Read-Backed Haplotype Phasing and BAM Haplotagging for nCATS Target Regions.
"""

import os
import pysam
import pandas as pd
import numpy as np
from collections import defaultdict

BAM_PATH = "results/alignment/ncats.sorted.bam"
VAR_PATH = "results/variants/variant_calls_summary.tsv"
OUT_DIR = "results/phasing"

os.makedirs(OUT_DIR, exist_ok=True)

def run_phasing():
    print("[1/5] Loading high-confidence variants for phasing...")
    df_vars = pd.read_csv(VAR_PATH, sep="\t")
    # Only phase high-confidence heterozygous SNVs
    df_het = df_vars[df_vars["GT"] == "0/1"].copy()
    print(f"  -> Found {len(df_het)} candidate heterozygous SNVs across all targets.")

    bam = pysam.AlignmentFile(BAM_PATH, "rb")

    all_phased_vars = []
    phase_blocks = []
    read_assignments = []

    # Prepare phased BAM writer
    phased_bam_path = f"{OUT_DIR}/ncats.phased.bam"
    phased_bam = pysam.AlignmentFile(phased_bam_path, "wb", template=bam)

    read_haplotags = {}  # read_name -> (HP, PS, v1, v2)

    print("[2/5] Building read-variant graph and solving phase blocks per target...")

    for ref_idx, ref_name in enumerate(bam.references):
        tid = f"T{ref_idx + 1}"
        t_vars = df_het[df_het["Target_ID"] == tid].sort_values("Target_Pos").reset_index(drop=True)

        if len(t_vars) < 2:
            print(f"  -> {tid} has fewer than 2 heterozygous SNVs; skipping multi-site phasing.")
            for _, row in t_vars.iterrows():
                r_dict = row.to_dict()
                r_dict["Phase_GT"] = "0|1"
                r_dict["Phase_Block"] = row["Genomic_Pos"]
                all_phased_vars.append(r_dict)
            continue

        span_bp = t_vars["Genomic_Pos"].max() - t_vars["Genomic_Pos"].min()
        print(f"  -> Phasing {tid} ({len(t_vars)} het SNVs, span: {span_bp:,} bp)...")

        # Build variant coordinate lookup: 0-indexed target pos -> variant index
        var_pos_to_idx = {row["Target_Pos"] - 1: idx for idx, row in t_vars.iterrows()}
        var_list = t_vars.to_dict("records")

        # Read alleles: read_name -> {var_idx: allele (0=ref, 1=alt)}
        read_alleles = defaultdict(dict)

        for pileupcolumn in bam.pileup(ref_name, truncate=True, min_base_quality=10):
            pos = pileupcolumn.pos
            if pos not in var_pos_to_idx:
                continue
            v_idx = var_pos_to_idx[pos]
            v_ref = var_list[v_idx]["Ref"]
            v_alt = var_list[v_idx]["Alt"]

            for p_read in pileupcolumn.pileups:
                if not p_read.is_del and not p_read.is_refskip:
                    q_base = p_read.alignment.query_sequence[p_read.query_position].upper()
                    if q_base == v_ref:
                        read_alleles[p_read.alignment.query_name][v_idx] = 0
                    elif q_base == v_alt:
                        read_alleles[p_read.alignment.query_name][v_idx] = 1

        # Keep reads covering >= 2 variants
        spanning_reads = {r: alleles for r, alleles in read_alleles.items() if len(alleles) >= 2}
        print(f"     Found {len(spanning_reads)} spanning reads linking heterozygous sites.")

        # Minimum Error Correction / Co-occurrence Phasing
        n_vars = len(var_list)
        link_same = defaultdict(lambda: defaultdict(int))
        link_diff = defaultdict(lambda: defaultdict(int))

        for rname, alleles in spanning_reads.items():
            items = list(alleles.items())
            for i in range(len(items)):
                for j in range(i + 1, len(items)):
                    v1, a1 = items[i]
                    v2, a2 = items[j]
                    if a1 == a2:
                        link_same[v1][v2] += 1
                        link_same[v2][v1] += 1
                    else:
                        link_diff[v1][v2] += 1
                        link_diff[v2][v1] += 1

        # Graph bipartition / BFS phase solver
        phased_state = {}  # var_idx -> 0 or 1
        block_id_map = {}  # var_idx -> block_start_pos
        visited = set()

        for start_idx in range(n_vars):
            if start_idx in visited:
                continue

            # Start new phase block
            block_start_pos = var_list[start_idx]["Genomic_Pos"]
            queue = [start_idx]
            phased_state[start_idx] = 0
            block_id_map[start_idx] = block_start_pos
            visited.add(start_idx)

            while queue:
                curr = queue.pop(0)
                curr_phase = phased_state[curr]

                all_neighbors = set(link_same[curr].keys()) | set(link_diff[curr].keys())
                for nbr in sorted(all_neighbors):
                    if nbr not in visited:
                        same_count = link_same[curr][nbr]
                        diff_count = link_diff[curr][nbr]
                        if same_count + diff_count >= 1:
                            if same_count >= diff_count:
                                nbr_phase = curr_phase
                            else:
                                nbr_phase = 1 - curr_phase
                            phased_state[nbr] = nbr_phase
                            block_id_map[nbr] = block_start_pos
                            visited.add(nbr)
                            queue.append(nbr)

        # Record phased variants
        for idx, row in t_vars.iterrows():
            r_dict = row.to_dict()
            phase = phased_state.get(idx, 0)
            ps = block_id_map.get(idx, row["Genomic_Pos"])
            r_dict["Phase_GT"] = "0|1" if phase == 0 else "1|0"
            r_dict["Phase_Block"] = ps
            all_phased_vars.append(r_dict)

        # Calculate phase block statistics
        block_groups = defaultdict(list)
        for idx in range(n_vars):
            ps = block_id_map.get(idx, var_list[idx]["Genomic_Pos"])
            block_groups[ps].append(var_list[idx]["Genomic_Pos"])

        for ps, positions in block_groups.items():
            b_span = max(positions) - min(positions) + 1
            phase_blocks.append({
                "Target_ID": tid,
                "Target_Name": t_vars.iloc[0]["Target_Name"],
                "Chr": t_vars.iloc[0]["Chr"],
                "Phase_Block_ID": ps,
                "Phased_SNVs_Count": len(positions),
                "Block_Start": min(positions),
                "Block_End": max(positions),
                "Block_Span_bp": b_span,
                "Spanning_Reads": len(spanning_reads)
            })

        # Assign reads to Haplotypes (H1 vs H2)
        for rname, alleles in read_alleles.items():
            h1_votes = 0
            h2_votes = 0
            block_votes = defaultdict(int)

            for v_idx, allele in alleles.items():
                if v_idx in phased_state:
                    v_phase = phased_state[v_idx]
                    ps = block_id_map[v_idx]
                    block_votes[ps] += 1
                    if v_phase == 0:
                        if allele == 0:
                            h1_votes += 1
                        else:
                            h2_votes += 1
                    else:
                        if allele == 1:
                            h1_votes += 1
                        else:
                            h2_votes += 1

            if block_votes:
                best_ps = max(block_votes.items(), key=lambda x: x[1])[0]
                if h1_votes > h2_votes:
                    read_haplotags[rname] = (1, best_ps, h1_votes, h2_votes)
                elif h2_votes > h1_votes:
                    read_haplotags[rname] = (2, best_ps, h1_votes, h2_votes)
                else:
                    read_haplotags[rname] = (0, best_ps, h1_votes, h2_votes)

    print(f"[3/5] Writing haplotagged BAM ({len(read_haplotags)} reads classified)...")

    h1_count = 0
    h2_count = 0
    unphased_count = 0

    for read in bam.fetch():
        rname = read.query_name
        if rname in read_haplotags:
            hp, ps, v1, v2 = read_haplotags[rname]
            if hp in (1, 2):
                read.set_tag("HP", hp, "i")
                read.set_tag("PS", ps, "i")
                if hp == 1:
                    h1_count += 1
                else:
                    h2_count += 1
                read_assignments.append({
                    "Read_Name": rname,
                    "Reference": read.reference_name,
                    "Haplotype": f"H{hp}",
                    "Phase_Block": ps,
                    "H1_Votes": v1,
                    "H2_Votes": v2,
                    "MapQ": read.mapping_quality,
                    "Length": read.query_length
                })
            else:
                unphased_count += 1
        else:
            unphased_count += 1
        phased_bam.write(read)

    bam.close()
    phased_bam.close()

    print("[4/5] Indexing phased BAM...")
    pysam.index(phased_bam_path)

    # Save summary tables
    df_phased_vcf = pd.DataFrame(all_phased_vars)
    df_phased_vcf.to_csv(f"{OUT_DIR}/phased_variants.tsv", sep="\t", index=False)

    df_blocks = pd.DataFrame(phase_blocks)
    df_blocks.to_csv(f"{OUT_DIR}/phase_blocks_summary.tsv", sep="\t", index=False)

    df_reads = pd.DataFrame(read_assignments)
    df_reads.to_csv(f"{OUT_DIR}/read_haplotype_assignments.tsv", sep="\t", index=False)

    # Write phased VCF
    write_phased_vcf(f"{OUT_DIR}/ncats_phased_variants.vcf", df_phased_vcf)

    print(f"[5/5] Phasing complete.")
    print(f"  -> Haplotype 1 Reads: {h1_count:,}")
    print(f"  -> Haplotype 2 Reads: {h2_count:,}")
    print(f"  -> Unphased Reads:    {unphased_count:,}")
    print(f"  -> Phase Blocks Table:\n{df_blocks.to_string()}")


def write_phased_vcf(vcf_path, df):
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=nCATS_Nanopore_WhatsHap_Phaser\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (Phased)\">\n")
        f.write("##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set Identifier\">\n")
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
            flt = row.get("Filter", "PASS")
            ps = row["Phase_Block"]
            gt = row["Phase_GT"]
            info = f"DP={row['DP']};AD={row['AD']};VAF={row['VAF']};TID={row['Target_ID']}"
            fmt = "GT:PS:DP:AD"
            sample = f"{gt}:{ps}:{row['DP']}:{row['AD']}"
            f.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\t{fmt}\t{sample}\n")


if __name__ == "__main__":
    run_phasing()
