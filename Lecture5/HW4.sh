# BMMB 852 HW4
# Set the accession number and dataset directory
ACCESSION="GCF_000013425.1"
DATASET="ncbi_dataset/data/$ACCESSION"
ZIP_FILE="ncbi_dataset.zip"
GENOME_FILE="GCF_000013425.1_ASM1342v1_genomic.fna"
# Download the genome with additional feature types
datasets download genome accession ${ACCESSION} --include gff3,cds,protein,rna,genome

# Unzip the file
unzip $ZIP_FILE

# Extract features of type "gene"
awk '$3 == "gene"' ${DATASET}/genomic.gff > genes.gff

# Extract features of type "CDS"
awk '$3 == "CDS"' ${DATASET}/genomic.gff > cds.gff

# Link genome file under new name but first remove any file that might have the same name
rm -f genome.fa
ln -s ${GENOME_FILE} genome.fa

# Build index for the genome. This is where I could not get the code to reproduce.
samtools faidx genome.fa
