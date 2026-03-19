import pandas as pd
import ast

# Load sORF table
df = pd.read_csv("stage1_cleaned_sorfs.csv")

bed_lines = []

for idx, row in df.iterrows():
    orf_id = f"ORF_{idx}"
    strand = row["strand"]
    segments = ast.literal_eval(row["genomic_segments"])
    
    for chrom, start, end in segments:
        bed_lines.append([
            f"chr{chrom}" if not str(chrom).startswith("chr") else chrom,
            start,
            end,
            orf_id,
            0,
            strand
        ])

bed_df = pd.DataFrame(bed_lines)
bed_df.to_csv("sorfs_genomic.bed", sep="\t", header=False, index=False)

print("BED file created: sorfs_genomic.bed")
