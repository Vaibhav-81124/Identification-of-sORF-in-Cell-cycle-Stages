import pandas as pd

INPUT = "stage1_novel_sorfs.csv"
OUTPUT = "stage1_cleaned_sorfs.csv"

print("Loading data...")
df = pd.read_csv(INPUT)

print(f"Initial ORFs: {len(df):,}")

# ======================================================
# STEP 1 — Remove exact duplicate peptides
# ======================================================

df_unique_pep = df.drop_duplicates(subset=["aa_sequence"])
print(f"After removing identical peptides: {len(df_unique_pep):,}")

# ======================================================
# STEP 2 — Remove duplicate genomic loci
# ======================================================

df_unique_loci = df_unique_pep.drop_duplicates(
    subset=["chromosome", "strand", "transcript_start", "transcript_end"]
)
print(f"After removing duplicate loci: {len(df_unique_loci):,}")

# ======================================================
# STEP 3 — Remove nested ORFs (keep longest per region)
# ======================================================

df_unique_loci["region_key"] = (
    df_unique_loci["chromosome"].astype(str) + "_" +
    df_unique_loci["strand"].astype(str)
)

df_final = []

for key, group in df_unique_loci.groupby("region_key"):
    
    group = group.sort_values("aa_length", ascending=False)
    kept = []
    
    for _, row in group.iterrows():
        is_nested = False
        
        for kept_row in kept:
            if (
                row["transcript_start"] >= kept_row["transcript_start"] and
                row["transcript_end"] <= kept_row["transcript_end"]
            ):
                is_nested = True
                break
        
        if not is_nested:
            kept.append(row)
    
    df_final.extend(kept)

df_final = pd.DataFrame(df_final)
print(f"After removing nested ORFs: {len(df_final):,}")

# ======================================================
# STEP 4 — Raise minimum length to 15 aa
# ======================================================

df_final = df_final[df_final["aa_length"] >= 15]
print(f"After applying 15 aa minimum: {len(df_final):,}")

df_final.to_csv(OUTPUT, index=False)

print(f"\nCleaned dataset saved to: {OUTPUT}")
