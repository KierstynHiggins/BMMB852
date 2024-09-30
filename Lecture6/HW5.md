# BMMB 852: HW5
## Kierstyn Higgins

### Part 1
Set the accession number and dataset directory.
```bash
ACCESSION="GCF_000013425.1"
ZIP_FILE="ncbi_dataset.zip"
FILE="genome.fa" 
```
Download the genome with additional feature types.
```bash
datasets download genome accession ${ACCESSION} --include gff3,cds,protein,rna,genome
```
Save under new name.
```bash
ln -sf ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna genome.fa 
```
Size of the file.
```bash
FILE_SIZE=$(wc -c < "$FILE")
echo "Size of the file: $FILE_SIZE bytes"
```
This will print:
```bash
echo "Size of the file: $FILE_SIZE bytes"
Size of the file:  2856716 bytes
```
Total size of the genome.
```bash
TOTAL_GENOME_SIZE=$(grep -v '^>' "$FILE" | tr -d '\n' | wc -c)
echo "Total size of the genome: $TOTAL_GENOME_SIZE bases"
```
This will print:
```bash
echo "Total size of the genome: $TOTAL_GENOME_SIZE bases"
Total size of the genome:  2821361 bases
```
Number of chromosomes in the genome.
```bash
UM_CHROMOSOMES=$(grep -c '^>' "$FILE")
echo "Number of chromosomes: $NUM_CHROMOSOMES"
```
This will print;
```bash
echo "Number of chromosomes: $NUM_CHROMOSOMES"
Number of chromosomes: 1
```
Name (id) and length of each chromosome.
```bash
echo "Chromosome lengths:"
grep '^>' "$FILE" | while read -r line; do
    CHR_NAME=$(echo "$line" | sed 's/>//')
    CHR_LENGTH=$(awk -v chr="$CHR_NAME" '/^>/{if (seq) print length(seq); seq=""; if ($0 ~ chr) seq=""; next} {seq=seq$0} END{print length(seq)}' "$FILE")
    echo "$CHR
```
This will print:
```bash
NC_007795.1 Staphylococcus aureus subsp. aureus NCTC 8325 chromosome, complete genome: 2821361 bases
```
The file is 2856716 bytes.
The genome is 2821361 bases.
THere is 1 chromosome, NC_007795.1: 2821361 bases.

### Part 2

Set read length.
```bash
READ_LENGTH=150
```
Calculate genome length.
```bash
GENOME_LENGTH=$(grep -v '^>' "$FILE" | tr -d '\n' | wc -c)
```
Set target coverage.
```bash
TARGET_COVERAGE=10
```
Calculate number of reads needed.
```bash
NUM_READS=$(( (TARGET_COVERAGE * GENOME_LENGTH) / READ_LENGTH ))
```
Generate simulated FASTQ files.
```bash
wgsim -N $NUM_READS -1 $READ_LENGTH -2 $READ_LENGTH -e 0.02 -r 0.05 -R 0.15 -X 0.3 "$REFERENCE" simulated_reads_1.fq simulated_reads_2.fq

echo "Simulated FASTQ files generated: simulated_reads_1.fq and simulated_reads_2.fq"
echo "Number of reads generated: $NUM_READS"
echo "Average read length: $READ_LENGTH bp"
```
This will print:
```bash
echo "Number of reads generated: $NUM_READS"
Number of reads generated: 188090
```
and
```bash
echo "Average read length: $READ_LENGTH bp"
Average read length: 150 bp
```
Get the size of the FASTQ files.
```bash
SIZE_BEFORE=$(du -sh simulated_reads_1.fq simulated_reads_2.fq)
echo "Size of the FASTQ files before compression:"
echo "$SIZE_BEFORE"
```
This will print:
```bash
echo "$SIZE_BEFORE"
 64M	simulated_reads_1.fq
 64M	simulated_reads_2.fq
 ```
Compress the FASTQ files.
```bash
gzip simulated_reads_1.fq simulated_reads_2.fq
```
Get the size of the compressed files.
```bash
SIZE_AFTER=$(du -sh simulated_reads_1.fq.gz simulated_reads_2.fq.gz)
echo "Size of the FASTQ files after compression:"
echo "$SIZE_AFTER"
```
This will print:
```bash
echo "$SIZE_AFTER"
 12M	simulated_reads_1.fq.gz
 12M	simulated_reads_2.fq.gz
 ```
 It prints 188090 reads with and average read length of 150 bp. The size before compression is 64 M down to 12 M, saving 52 M of space. If you decrease read length you need more reads to keep the same coverage, and increasing read length will require fewer reads.

 ### Part 3
 Approximate genome size for:
 Yeast - 12 Mb
 Drosophila - 180 Mb
 Human - 3.1 Gb

 If we want to calculate 30x coverage, we can just multiply genome size by 30.

 Yeast - 360 Mb
 Drosophila - 5.4 Gb
 Human - 93 Gb

Yeast:
FASTA file size: 12 Mb
Number of reads: 2,400,000
FASTQ size before compression: approx 900 Mb
FASTQ size after compression: approx 450 Mb

Drosophila:
FASTA file size: 180 Mb
Number of reads: 36,000,000
FASTQ size before compression: approx 13.5 Gb
FASTQ size after compression: approx 6.75 Gb

Human:
FASTA file size: 3.1 Gb
Number of reads: 620,000,000
FASTQ size before compression: approx 232.5 Gb
FASTQ size after compression: approx 116.25 Gb



