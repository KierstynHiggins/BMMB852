# RNA seq Makefile
## Kierstyn Higgins
### 12/08/2024

#### Set variables here.

# The project name for reads.
PRJ=PRJNA699466

# SRR NUMBER.
SRR = SRR13628155

# Accession number of the Staphylococcus aureus genome.
ACC=GCF_000013425.1

# Genome reference.
REF=ncbi_dataset/refs/staph-2006.fa

# Genome annotation file.
GTF=ncbi_dataset/refs/staph-2006.gtf

# The counts in tab delimited format.
COUNTS_TXT = res/counts-hisat.txt

# Final combinted counts in CSV format.
COUNTS = res/counts-hisat.csv

# The design file
DESIGN = design.csv

# Sample alias.
SAMPLE = GSM5064589

# Extract the sample names from the first 3 rows of design.csv (adjust head if needed)
SAMPLES=$(head -3 design.csv | cut -d',' -f1)

# Output csv.
OUTPUT_CSV="combined_data.csv"

# Flags passed to parallel.
FLAGS = --eta --lb --header : --colsep ,

# Number of CPUS to use
NCPU = 4
# Number of reads.
N = 5000

# The GFF file.
GFF=refs/staph-2006.gff

# Directories
BAM_DIR=bam

# The path to read 1
R1=reads/${SRR}_1.fastq

# The path to read 2
R2=reads/${SRR}_2.fastq

# The resulting BAM file.
BAM=bam/${SRR}.bam
#_______________________________________________________________________________________
# Makefile settings
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory
# ____________________________________________________________________________________
# Print usage
usage:
	
	@echo "# make design (generates the ${DESIGN} file)"
	@echo "# make bam (aligns reads into BAM files)"
	@echo "# make count  (counts the reads in the BAM files)"
	@echo "# make csv (converts txt file to csv format)
	@echo "# make all    (runs all steps)"

# Generate the design file.
design: 
	bio search ${PRJ} -H --csv > ${DESIGN}


# Create the BAM alignment file.
bam: ${DESIGN}
	# Get the reference genome and the annotations.
	make -f src/run/genbank.mk ACC=${ACC} REF=${REF} GFF=${GFF} fasta gff

	# Index the reference genome.
	make -f src/run/hisat2.mk REF=${REF} index

	# Download the sequence data using fastq-dump.
	make fastq-dump

	# Align the reads to the reference genome. 
	# Use a sample name in the readgroup.
	make -f src/run/hisat2.mk SM=${SAMPLE} REF=${REF} R1=${R1} R2=${R2} BAM=${BAM} run 

# Download the sequence data using fastq-dump
fastq-dump:
	# Check if the 'reads' directory exists, if not create it
	mkdir -p reads
	
	# Run fastq-dump to fetch the paired-end reads and split them
	# --split-files will generate two separate fastq files: R1 and R2
	# --outdir specifies the directory where the fastq files will be stored
	echo "Downloading ${SRR} using fastq-dump..."
	fastq-dump --split-files --outdir reads ${SRR}

	# Optionally, limit the number of reads if needed (N = 50000)
	if [ -n "$(N)" ]; then
		head -n $(N) reads/${SRR}_1.fastq > reads/${SRR}_1_subset.fastq
		head -n $(N) reads/${SRR}_2.fastq > reads/${SRR}_2_subset.fastq
		mv reads/${SRR}_1_subset.fastq reads/${SRR}_1.fastq
		mv reads/${SRR}_2_subset.fastq reads/${SRR}_2.fastq
	fi

	# Verify that the files have been downloaded
	ls -lh reads/${SRR}_1.fastq reads/${SRR}_2.fastq

# Generate count data.
count: ${DESIGN}
	# Create the output directory
	mkdir -p res/

	# Count features for each sample
	cat design.csv | grep -w ${SAMPLE} | cut -d',' -f1 | while read sample; do \
		# Print the sample name being processed
		echo "Processing sample: ${sample}"; \
		featureCounts -a ${GTF} -o ${COUNTS_TXT} ${BAM_DIR}/${SRR}.bam; \
	done

# The final counts in CSV format.
csv: ${COUNTS_TXT}
	micromamba run -n stats Rscript src/r/format_featurecounts.r -c ${COUNTS_TXT} -o ${COUNTS}

# Run all.
all: design bam count csv




















