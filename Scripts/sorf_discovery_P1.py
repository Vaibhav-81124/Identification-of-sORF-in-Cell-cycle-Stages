import os
import re
from collections import defaultdict
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# ==============================
# CONFIGURATION
# ==============================

FASTA_FILE = "Homo_sapiens.GRCh38.cdna.all.fa"
GTF_FILE = "Homo_sapiens.GRCh38.115.gtf"

MIN_AA = 10
MAX_AA = 100
STOP_CODONS = {"TAA", "TAG", "TGA"}

OUTPUT_FILE = "stage1_novel_sorfs.csv"

# ==============================
# DATA STRUCTURES
# ==============================

class Transcript:
    def __init__(self, chrom, strand):
        self.chrom = chrom
        self.strand = strand
        self.exons = []  # (start, end)
        self.cds = []    # (start, end)

# ==============================
# STEP 1: PARSE GTF
# ==============================

print("Parsing GTF...")

transcripts = {}

with open(GTF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        
        chrom, feature, start, end, strand = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6]
        
        match = re.search(r'transcript_id "([^"]+)"', fields[8])
        if not match:
            continue
        
        tid = match.group(1)
        
        if tid not in transcripts:
            transcripts[tid] = Transcript(chrom, strand)
        
        if feature == "exon":
            transcripts[tid].exons.append((start, end))
        
        if feature == "CDS":
            transcripts[tid].cds.append((start, end))

print(f"Loaded {len(transcripts)} transcripts")

# Sort exons properly
for t in transcripts.values():
    t.exons.sort(key=lambda x: x[0])

# ==============================
# STEP 2: LOAD FASTA
# ==============================

print("Loading transcript sequences...")

fasta_sequences = {}
for record in SeqIO.parse(FASTA_FILE, "fasta"):
    tid = record.id.split("|")[0]
    fasta_sequences[tid] = str(record.seq)

print(f"Loaded {len(fasta_sequences)} sequences")

# ==============================
# HELPER: MAP TRANSCRIPT → GENOME
# ==============================

def transcript_to_genome(tid: str, t_start: int, t_end: int) -> List[Tuple[str,int,int]]:
    """
    Convert transcript coordinates to genomic segments.
    Returns list of genomic segments (chrom, start, end).
    """
    t = transcripts[tid]
    exons = t.exons
    chrom = t.chrom
    strand = t.strand
    
    segments = []
    remaining_start = t_start
    remaining_end = t_end
    
    transcript_cursor = 0
    
    for exon_start, exon_end in exons:
        exon_length = exon_end - exon_start + 1
        
        exon_t_start = transcript_cursor
        exon_t_end = transcript_cursor + exon_length
        
        # Overlap with ORF
        if remaining_end <= exon_t_start:
            break
        
        if remaining_start >= exon_t_end:
            transcript_cursor += exon_length
            continue
        
        overlap_start = max(remaining_start, exon_t_start)
        overlap_end = min(remaining_end, exon_t_end)
        
        offset_start = overlap_start - exon_t_start
        offset_end = overlap_end - exon_t_start
        
        if strand == "+":
            g_start = exon_start + offset_start
            g_end = exon_start + offset_end - 1
        else:
            g_end = exon_end - offset_start
            g_start = exon_end - offset_end + 1
        
        segments.append((chrom, g_start, g_end))
        
        transcript_cursor += exon_length
    
    return segments

# ==============================
# STEP 3: FIND sORFs
# ==============================

print("Scanning for sORFs...")

results = []
total = len(fasta_sequences)

for idx, (tid, seq) in enumerate(fasta_sequences.items(), 1):
    
    # Progress indicator
    if idx % 5000 == 0:
        print(f"Processed {idx}/{total} transcripts | sORFs found: {len(results)}")
    
    if tid not in transcripts:
        continue
    
    # Skip very short transcripts
    if len(seq) < MIN_AA * 3:
        continue
    
    # Skip transcripts without ATG
    if "ATG" not in seq:
        continue
    
    seq = seq.upper()
    seq_len = len(seq)
    transcript_data = transcripts[tid]
    cds_regions = transcript_data.cds
    
    for frame in range(3):
        i = frame
        
        while i < seq_len - 2:
            codon = seq[i:i+3]
            
            if codon != "ATG":
                i += 3
                continue
            
            j = i + 3
            
            while j < seq_len - 2:
                stop = seq[j:j+3]
                
                if stop in STOP_CODONS:
                    nt_seq = seq[i:j+3]
                    aa_len = (j+3 - i) // 3 - 1
                    
                    if MIN_AA <= aa_len <= MAX_AA:
                        
                        genomic_segments = transcript_to_genome(tid, i, j+3)
                        
                        overlaps = False
                        for cds_start, cds_end in cds_regions:
                            for chrom, g_start, g_end in genomic_segments:
                                if g_start <= cds_end and g_end >= cds_start:
                                    overlaps = True
                                    break
                            if overlaps:
                                break
                        
                        if not overlaps:
                            results.append({
                                "transcript_id": tid,
                                "chromosome": transcript_data.chrom,
                                "strand": transcript_data.strand,
                                "transcript_start": i,
                                "transcript_end": j+3,
                                "aa_length": aa_len,
                                "genomic_segments": genomic_segments,
                                "nt_sequence": nt_seq,
                                "aa_sequence": str(Seq(nt_seq).translate(to_stop=True))
                            })
                    
                    i = j + 3
                    break
                
                j += 3
            else:
                i += 3
# ==============================
# STEP 4: SAVE RESULTS
# ==============================

df = pd.DataFrame(results)
df.to_csv(OUTPUT_FILE, index=False)

print(f"Saved to {OUTPUT_FILE}")
print("Stage 1 complete.")
