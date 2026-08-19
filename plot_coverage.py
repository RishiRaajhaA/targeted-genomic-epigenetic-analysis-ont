import pandas as pd
import matplotlib.pyplot as plt
import os

INPUT = "results/coverage/target_coverage.tsv"
OUTPUT_DIR = "results/coverage"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Read SAMtools coverage output
df = pd.read_csv(INPUT, sep="\t")

# Give short names to the targets
df["target"] = [f"T{i+1}" for i in range(len(df))]

# -----------------------------
# Plot 1: Mean sequencing depth
# -----------------------------

plt.figure(figsize=(10, 6))

plt.bar(df["target"], df["meandepth"])

plt.xlabel("Target region")
plt.ylabel("Mean sequencing depth (×)")
plt.title("Mean Sequencing Depth Across Target Regions")

plt.tight_layout()

plt.savefig(
    f"{OUTPUT_DIR}/mean_depth.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# -----------------------------
# Plot 2: Coverage breadth
# -----------------------------

plt.figure(figsize=(10, 6))

plt.bar(df["target"], df["coverage"])

plt.xlabel("Target region")
plt.ylabel("Covered bases (%)")
plt.title("Coverage Breadth Across Target Regions")

plt.ylim(99, 100.1)

plt.tight_layout()

plt.savefig(
    f"{OUTPUT_DIR}/coverage_breadth.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print("Plots generated:")
print(f"  {OUTPUT_DIR}/mean_depth.png")
print(f"  {OUTPUT_DIR}/coverage_breadth.png")