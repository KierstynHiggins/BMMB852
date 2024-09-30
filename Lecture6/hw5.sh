# Set the accession number and dataset directory
ACCESSION="GCF_000013425.1"
ZIP_FILE="ncbi_dataset.zip"
FILE="genome.fa" 
# Download the genome with additional feature types
datasets download genome accession ${ACCESSION} --include gff3,cds,protein,rna,genome

# Unzip file
unzip -o $ZIP_FILE

# Save under new name
ln -sf ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna genome.fa 

# Size of the file
FILE_SIZE=$(wc -c < "$FILE")
echo "Size of the file: $FILE_SIZE bytes"

# Total size
TOTAL_GENOME_SIZE=$(grep -v '^>' "$FILE" | tr -d '\n' | wc -c)
echo "Total size of the genome: $TOTAL_GENOME_SIZE bases"

# Number of chromosomes
NUM_CHROMOSOMES=$(grep -c '^>' "$FILE")
echo "Number of chromosomes: $NUM_CHROMOSOMES"

# Name and length of each chromosome
echo "Chromosome lengths:"
grep '^>' "$FILE" | while read -r line; do
    CHR_NAME=$(echo "$line" | sed 's/>//')
    CHR_LENGTH=$(awk -v chr="$CHR_NAME" '/^>/{if (seq) print length(seq); seq=""; if ($0 ~ chr) seq=""; next} {seq=seq$0} END{print length(seq)}' "$FILE")
    echo "$CHR_NAME: $CHR_LENGTH bases"
done


# Read length (assuming 150 bp)
READ_LENGTH=150

# Calculate genome length
GENOME_LENGTH=$(grep -v '^>' "$FILE" | tr -d '\n' | wc -c)

# Target coverage
TARGET_COVERAGE=10

# Calculate number of reads needed
NUM_READS=$(( (TARGET_COVERAGE * GENOME_LENGTH) / READ_LENGTH ))

# Generate simulated FASTQ files
wgsim -N $NUM_READS -1 $READ_LENGTH -2 $READ_LENGTH -e 0.02 -r 0.05 -R 0.15 -X 0.3 "$REFERENCE" simulated_reads_1.fq simulated_reads_2.fq

echo "Simulated FASTQ files generated: simulated_reads_1.fq and simulated_reads_2.fq"
echo "Number of reads generated: $NUM_READS"
echo "Average read length: $READ_LENGTH bp"

# Get the size of the FASTQ files
SIZE_BEFORE=$(du -sh simulated_reads_1.fq simulated_reads_2.fq)
echo "Size of the FASTQ files before compression:"
echo "$SIZE_BEFORE"

# Compress the FASTQ files
gzip simulated_reads_1.fq simulated_reads_2.fq

# Get the size of the compressed files
SIZE_AFTER=$(du -sh simulated_reads_1.fq.gz simulated_reads_2.fq.gz)
echo "Size of the FASTQ files after compression:"
echo "$SIZE_AFTER"



