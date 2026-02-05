#########################################################################################

 o-o  o   o o-O-o       O              o-o           
o   o |\  |   |        / \            |              
|   | | \ |   |   o-o o---oo-O-o o-o   o-o  o-o  o-o 
o   o |  \|   |       |   || | | |  |     | |-' |  | 
 o-o  o   o   o       o   oo o o O-o  o--o  o-o  o-O 
                                 |                 | 
                                 o                 o 
 
 
##ONT-Ampseq Adapted Script (Sofia Olivier & Kathryn Morrissey - April 2025)

#########################################################################################

##STEP 01-03: DOWNLOAD, BASE CALL & QC REPORT OF RAW SEQUENCE FILES (Optional pre-processing of raw reads)
###Step 01&02 - Download raw sequences from sequencing device and base call using Dorado (Note: run on an HPC or on a local machine with ≥32GB RAM processing power)



###Step 03 - Create QC report using PycoQC



##STEP 04-07: SEQUENCE TRIMMING, QUALITY FILTERING, DEREPLICATION AND CHIMERA REMOVAL
###Step 04 - Trim adapters off the end of sequences 

for i in *.fastq;
	do
    porechop -i "$i" -o "${i%.fastq}_trimmed.fastq"
	done


###Step 05 - Quality and length filtering and conversion to .fasta format

for i in *_trimmed.fastq;
	do
    base="${i%_trimmed.fastq}"
    NanoFilt -q 9 -l 320 --maxlength 480 < "$i" > "${base}_filtered.fastq"
	done

for i in *_filtered.fastq;
	do
    base="${i%_filtered.fastq}"
    seqtk seq -a "$i" > "${base}.fasta"
	done

###Step 06 - Dereplication of sequences (grouping identical reads)

for i in *.fasta;
	do
    base="${i%.fasta}"
    vsearch --fastx_uniques "$i" --fastaout "${base}_uniques.fasta" --sizeout
	done
	
###Step 07 - Chimera identification and removal 

for i in *.fasta;
	do
    base="${i%.fasta}"
    vsearch --uchime2_denovo "$i" --uchimeout "${base}_nonchimeras.fasta" 
	done


######################################################################################### 
 _     _ ______   ______ 
| |   | (_____ \ / _____)
| |__ | |_____) ) /      
|  __)| |  ____/| |      
| |   | | |     | \_____ 
|_|   |_|_|      \______)

#########################################################################################


##THE FOLLOWING CAN BE RUN ON AN HPC OR ON A LOCAL MACHINE WITH ≥32GB RAM 
###Install vsearch, minimap2 and racon before continuing with the next steps and, for ease of workflow, make directories for each step 

cd /path/to/analysis/
mkdir 7_fastafiles         ##Transfer your quality control processed, non-chimeric sequence files from Step 07 into this directory
mkdir 8_centroids
mkdir 9_minimap
mkdir 10_racon
mkdir 11_relabel
mkdir 12_cluster_final
mkdir 13_tax_align


##STEP 08-10: SEQUENCE CLUSTERING, ALIGNMENT AND POLISHING
###Step 8 -  Sequence clustering and concatenating resulting centroids 

cd /path/to/7_fastafiles/

for i in *;
	do
	vsearch --cluster_unoise  /path/to/7_fastafiles/$i -id 0.90 --minsize 1 --centroids  /path/to/8_centroids/centroids_$i
	done

cat *fasta > cat_centroids.fa


###STEP 9 - Alignment
###SPECIFICATIONS CAN BE FOUND HERE: https://github.com/lh3/minimap2

cd /path/to/7_fastafiles/

for i in *;
	do 
	minimap2 -ax lr:hq /path/to/8_centroids/cat_centroids.fa /path/to/7_fastafiles/$i  > /path/to/9_minimap/minimap_$i.sam
	done
	

###STEP 10 - Polishing
###SPECIFICATIONS CAN BE FOUND HERE: https://github.com/lbcb-sci/racon

cd /path/to/7_fastafiles/

for i in *;
	do 
	racon -f -q 9 /path/to/8_centroids/cat_centroids.fa /path/to/9_minimap/minimap_$i.sam /path/to/7_fastafiles/$i  > /path/to/10_racon/polished_$i
	done


##STEP 11-13: RELABELLING, OTU GENERATION AND TAXONOMIC CLASSIFICATION
###STEP 11 - Relabel (this is used to relabel the reads in a barcode so that the reads can be linked back to the barcode)

cd /path/to/10_racon/

for i in polished_*;
	do
	base=$(basename "$i")
    bc_num="${base#polished_}"
    vsearch --sortbysize "$i" --relabel "BC${bc_num}:" --output /path/to/11_relabel/polished_bc${bc_num}_relabel.fasta
	done
	
cat *fasta > allbarcodes_polished.fa

###STEP 12 - OTU Generation (97% clustering similarity)

cd /path/to/12_cluster_final/

vsearch --cluster_fast /path/to/11_relabel/allbarcodes_polished.fa -id 0.97 --relabel OTU_ --sizeout --otutabout /path/to/12_cluster_final/polished_otutable.0.97.tsv --centroids /path/to/12_cluster_final/polished_otus_97.fa

###STEP 13 - Taxonomic classification 

cd /path/to/12_cluster_final/

vsearch --sintax polished_otus_97.fa --db sintaxformatted_db.fasta --tabbedout /path/to/13_tax_align/polished_taxonomic_classification.tsv 

