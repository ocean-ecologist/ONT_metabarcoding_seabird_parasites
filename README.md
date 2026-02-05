# Bioinformatic pipeline for ONT metabarcoding of seabird parasites

### Project description
This repository contains scripts, workflows, and datasets related to the molecular surveillance of apicomplexan parasites in seabirds, specifically the Endangered Cape cormorant and Cape gannet and the Critically Endangered African Penguin populations. Our goal is to develop and optimize bioinformatics tools for detecting and characterizing these parasites using molecular techniques. This pipeline has been adapted from the original source code created by Patrick Schacksen et al. (2024) on their GitHub repository https://github.com/MathiasEskildsen/ONT-AmpSeq/tree/main/workflow and follows on from the sequencing of amplicons of the 18S rRNA gene of apicomplexan blood parasites using Oxford Nanopore Technologies MinION Mk1C device. The general workflow, however, can be applied to any amplicons sequenced in the same manner and of varying lengths. In this project, sequences of around 378 bp and 1450 bp were processed. The pipeline starts off with raw sequences directly obtained from the MinION machine. The code embedded below is demonstrated for a single file. The full code as an iterative process, including loops for multiple files (barcodes), is available as a separate file. 

_________________________________________________________________________________________

 ## STEP 01-03: DOWNLOAD, BASECALL & QC REPORT OF RAW SEQUENCE FILES 
 
 ### 01 Download raw sequences from MinION 
 Sequences in the standard .fastq file output were downloaded from the sequencing machine (MinION, in this case) onto a laptop. Raw sequence data files are large and so a hard drive or local disk with a high storage capacity is required. Each .fastq file contained multiple sequences from a single barcode.

 ### 02 Basecall sequences
 All sequences were base-called using Dorado (https://github.com/nanoporetech/dorado ) version 0.9.6 with the DNA super accuracy configuration, dna_r10.4.1_e8.2_400bps_sup@v5.2.0. The code for this is detailed below. However, at this stage, other basecalling software such as Guppy may be used. 
```
 Code here
``` 

 ### 03 QC Report using pycoQC
 A pycoQC report (https://a-slide.github.io/pycoQC/) was used to check read quality (PHRED score), median read length, etc. This was used to inform parameters set during the subsequent QC processing of sequences. 
```
Code here
```

 ## STEP 04-07: SEQUENCE TRIMMING, QUALITY FILTERING, DEREPLICATION AND CHIMERA REMOVAL 
 
### 04 Sequence trimming
The tool Porechop (https://github.com/rrwick/Porechop) was used to trim adapters off the ends of the sequences in each .fastq file. 
```
porechop -i input_reads.fastq -o output_reads.fastq
```

### 05 Quality filtering
The tool NanoFilt (https://github.com/wdecoster/nanofilt) was used to filter out lower quality sequences with parameters chosen based off of the pycoQC report (in this case, quality threshold = 9; read length = 320-480bp for short reads and 1300-1700bp for long reads). The .fastq files were then converted to .fasta files for further processing.
```
NanoFilt -q 9 -l 1300 --maxlength 1700 input_reads.fastq > output_reads.fastq

seqtk seq -a input_reads.fastq > output_reads.fasta
```

### 06 Dereplication  
The quality filtered sequences underwent a dereplication step (grouping identical reads) using VSEARCH so that chimera detection could be performed. 
```
vsearch --fastx_uniques input_reads.fasta --fastaout output_reads.fasta
```

### 07 Chimera detection
Chimera detection refers to the process of identifying and removing artificial hybrid sequences generated during PCR. In this case, De novo chimera detection was performed using the UCHIME2 algorithm in VSEARCH. 
```
vsearch --uchime2_denovo input_reads.fasta --uchimeout output_reads.fasta
```

 ## STEP 08-10: SEQUENCE CLUSTERING, ALIGNMENT AND POLISHING 
The next steps took place on a HPC due to the processing power required, especially where numerous barcodes were involved. It could also be done on a laptop with >32GB RAM processing power. 
 
### 08 Sequence clustering
Draft consensus sequences were created from the quality reads using VSEARCH and a 90% identity threshold and a min abundance of 1 sequence per cluster. Specifically, a centroid-based approach was used per barcode. All the draft consensus sequences were then concatenated into a single file. 
```
vsearch --cluster_unoise input_reads.fasta -id 0.90 --minsize 1 --centroids centroids.fasta
cat *fasta > cat_centroids.fasta
```

### 09 Alignment 
To create alignments, the original filtered reads from STEP 06 were mapped to the draft consensus sequences from STEP 07 through minimap2 using the parameter for Nanopore Q20 genomic reads (v2.27+).
```
 minimap2 -ax lr:hq cat_centroids.fasta output_reads.fasta > aln.sam
```

### 10 Polishing
The alignments were then corrected using a polishing software called Racon to generate error corrected OTUs. This involved the use of the the draft consensus sequences from STEP 07, the alignment from STEP 08 and the original filtered reads from STEP 06. This was done using the fragment correction parameter and a quality threshold of 0.9. 
```
racon -f -q 9 cat_centroids.fasta aln.sam output_reads.fasta > polished_output_reads.fasta
```

 ## STEP 11-13: RELABELLING, OTU GENERATION AND TAXONOMIC CLASSIFICATION 

### 11 Relabelling
Polished sequences were sorted and given new headers indicating barcode and sequence number using VSEARCH. This was done so the final OTUs could be traced back to the sample that they came from. The .fasta files containing the relabelled sequences were then concatenated. 
```
vsearch --sortbysize "$input_file" --relabel "$relabel_prefix" --output "$output_file"
cat *fasta > cat_polished_output_reads.fasta
```

### 12 OTU Generation
A final clustering step was performed with VSEARCH to create the final biologically significant core OTUs. This time, a 97% identity threshold was utilised, the sequences were relabelled with OTU numbers and an OTU table was created.
```
vsearch --cluster_fast cat_polished_output_reads.fasta -id 0.97 --relabel OTU_ --sizeout --otutabout otutable.tsv --centroids centroids.fasta output_OTU_reads.fasta 
```

### 13 Taxonomic classification 
VSEARCH was used to taxonomically classify the polished OTUs using a SINTAX algorithm along with a custom SINTAX-formatted fasta database (see custom building of this database in this repository). 
```
vsearch --sintax output_OTU_reads.fasta --db custom_db.fasta --tabbedout taxonomic_classifications.tsv
```





 
