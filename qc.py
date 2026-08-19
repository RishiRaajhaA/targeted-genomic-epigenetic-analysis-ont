import gzip
import statistics
import math
import matplotlib.pyplot as plt

FASTQ = "SRR12423814_1.fastq.gz"
OUT = "results/qc"

lengths = []
mean_qualities = []

total_bases = 0
read_count = 0


def phred_quality(quality_string):
    """
    Convert ASCII Phred+33 quality characters to numeric scores.
    """
    return [(ord(char) - 33) for char in quality_string]


with gzip.open(FASTQ, "rt") as f:

    while True:

        header = f.readline()

        if not header:
            break

        sequence = f.readline().strip()
        plus = f.readline()
        quality = f.readline().strip()

        length = len(sequence)

        qualities = phred_quality(quality)
        mean_q = statistics.mean(qualities)

        lengths.append(length)
        mean_qualities.append(mean_q)

        total_bases += length
        read_count += 1


# -------------------------
# Basic statistics
# -------------------------

mean_length = statistics.mean(lengths)
median_length = statistics.median(lengths)

# N50
sorted_lengths = sorted(lengths, reverse=True)

half_bases = total_bases / 2
running_bases = 0

for length in sorted_lengths:

    running_bases += length

    if running_bases >= half_bases:
        n50 = length
        break


print("========== NANOPore QC ==========")

print(f"Reads          : {read_count:,}")
print(f"Total bases    : {total_bases:,}")
print(f"Mean length    : {mean_length:.2f} bp")
print(f"Median length  : {median_length:.2f} bp")
print(f"Min length     : {min(lengths):,} bp")
print(f"Max length     : {max(lengths):,} bp")
print(f"N50            : {n50:,} bp")

print()
print("========== QUALITY ==========")

print(f"Mean read Q    : {statistics.mean(mean_qualities):.2f}")
print(f"Median read Q  : {statistics.median(mean_qualities):.2f}")
print(f"Min read Q     : {min(mean_qualities):.2f}")
print(f"Max read Q     : {max(mean_qualities):.2f}")


# -------------------------
# Read length histogram
# -------------------------

plt.figure(figsize=(10, 6))

plt.hist(lengths, bins=100)

plt.xlabel("Read length (bp)")
plt.ylabel("Number of reads")
plt.title("Nanopore Read Length Distribution")

plt.tight_layout()

plt.savefig(
    f"{OUT}/read_length_distribution.png",
    dpi=300
)

plt.close()


# -------------------------
# Quality histogram
# -------------------------

plt.figure(figsize=(10, 6))

plt.hist(mean_qualities, bins=50)

plt.xlabel("Mean Phred quality score")
plt.ylabel("Number of reads")
plt.title("Nanopore Read Quality Distribution")

plt.tight_layout()

plt.savefig(
    f"{OUT}/quality_distribution.png",
    dpi=300
)

plt.close()


print()
print("Plots saved to:")
print(f"  {OUT}/read_length_distribution.png")
print(f"  {OUT}/quality_distribution.png")