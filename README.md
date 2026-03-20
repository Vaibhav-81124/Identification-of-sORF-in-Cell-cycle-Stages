# Discovery and Characterization of Novel sORFs Using Ribosome Profiling and RNA-seq in HeLa Cells

**Undergraduate Thesis Project (B.E. Biotechnology)**  
Author: Vaibhav D  

---

## Overview
This project focuses on the **discovery and validation of novel small open reading frames (sORFs)** using an integrative pipeline combining **RNA-seq and Ribo-seq data** across different stages (Asynchronous,M, and S phase) of the cell cycle in HeLa cells.

The workflow identifies candidate sORFs from transcriptomic data and applies multiple layers of filtering, translation evidence, and functional annotation to obtain **high-confidence translated microproteins**.

Benchling link
https://benchling.com/s/etr-7D2igJSt59EQClohaoSW?m=slm-j8SAZyjmnDi34mRHEaJM
---

## Dataset
- Source: PRJNA316618  
- Data type:
  - RNA-seq (transcript abundance)
  - Ribo-seq (translation evidence)  
- Conditions: Multiple cell-cycle stages (with replicates)  
- Reference:
  - Genome: GRCh38  
  - Annotation: GTF (Ensembl)  

---

## Objectives
- Identify **novel sORFs** from transcript sequences  
- Validate translation using **ribosome profiling (Ribo-seq)**  
- Filter high-confidence **microproteins (≤100 aa)**  
- Analyze **reading frame periodicity**  
- Annotate genomic context and host genes  
- Perform downstream **functional and network analysis**

---

## Pipeline Overview

### Phase 1: sORF Discovery
- Input:
  - Transcript FASTA (cDNA)
  - GTF annotation  
- Steps:
  - Parse transcripts and exon structure  
  - Map transcript coordinates to genome  
  - Scan for ORFs:
    - Start codon: ATG  
    - Length: 10–100 aa (later refined ≥15 aa)  
  - Remove:
    - CDS-overlapping ORFs  
    - Short/invalid transcripts  


---

### Phase 1B: sORF Filtering & Cleaning
- Remove:
  - Duplicate peptides  
  - Duplicate genomic loci  
  - Nested ORFs (retain longest per region)  
- Apply stricter cutoff:
  - Minimum length ≥15 aa  


---

### Phase 2: RNA-seq Processing
- Tools:
  - FastQC  
  - Cutadapt  
  - STAR  

- Steps:
  - Quality control  
  - Adapter trimming  
  - Alignment to genome  
  - Generate sorted BAM files  

- Output:
  - Aligned RNA-seq BAM files  

---

### Phase 3: Ribo-seq Processing
- Tools:
  - Cutadapt  
  - SortMeRNA  
  - STAR  
  - Samtools  

- Steps:
  - Adapter trimming  
  - rRNA removal  
  - Alignment  
  - Filter ribosome-protected fragments (RPFs: 26–34 nt)  
  - P-site assignment  

- Output:
  - P-site BED files representing translation events  

---

### Phase 4: Translation Evidence & Quantification

#### Ribo-seq Counts
- Convert sORFs → BED format  
- Intersect with P-sites (BEDTools)  
- Aggregate counts per ORF  

#### Filtering Pipeline
- Remove:
  - CDS-overlapping ORFs  
  - Protein-coding gene overlaps  
- Apply translation threshold:
  - Minimum ribo reads ≥10  

#### Microprotein Selection
- Define:
  - Microproteins ≤100 aa  

- Outputs:
  - High-confidence translated ORFs  
  - High-confidence microproteins  

---

### Phase 5: Periodicity Analysis
- Evaluate 3-nt periodicity of ribosome footprints  
- Compute frame distribution (frame0, frame1, frame2)  
- Apply thresholds:
  - Minimum reads ≥8  
  - Frame0 fraction ≥0.55  

- Output:
  - `high_confidence_translated_orfs`

---

### Phase 6: RNA Expression Quantification
- Convert sORFs → GTF  
- Use `featureCounts` for RNA expression  
- Integrate RNA + Ribo signals  

---

### Phase 7: Gene Annotation & Functional Analysis
- Identify **host genes** via genomic overlap  
- Classify ORFs:
  - Intragenic  
  - Intergenic  
  - Pseudogene-associated  
  - lncRNA-associated  

- Downstream analysis:
  - Pathway enrichment (g:Profiler)  
  - Protein interaction networks (STRING)  

---

## Key Features of the Pipeline
- Transcript-to-genome coordinate mapping  
- Multi-step filtering for high-confidence ORFs  
- Integration of **RNA-seq + Ribo-seq**  
- Translation validation via:
  - Ribo counts  
  - 3-nt periodicity  
- Identification of **novel microproteins**  

---

## Outputs
- Novel sORF dataset  
- High-confidence translated ORFs  
- Microprotein candidates  
- Host gene annotations  
- Functional enrichment results  

---

## Conclusion
This pipeline enables systematic identification and validation of **novel translated sORFs**, revealing previously unannotated coding potential in the human transcriptome. The integration of ribosome profiling and transcriptomics provides strong evidence for **functional microproteins** and their potential biological roles.
