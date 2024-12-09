# RNA seq Makefile
## Kierstyn Higgins

# The URL for the data.
DATA_URL = https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download?filename=ncbi_dataset.zip&ncbi_phid=322CB9919F68EEE500003635566DC9E4.1.m_3.020

# The downloaded data file name.
DATA_FILE = $(notdir ${DATA_URL})

# The project name for reads.
PRJ=PRJNA976947

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

# Makefile settings
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Flags passed to parallel.
FLAGS = --eta --lb --header : --colsep ,

# Number of CPUS to use
NCPU = 4
# ____________________________________________________________________________________
# Print usage
usage:
	
	@echo "# REF=${REF}"
	@echo "# GTF=${GTF}"
	@echo "# DESIGN=${DESIGN}"
	@echo "# COUNTS=${COUNTS}"
	@echo "#"

	@echo "# make data   (downloads the FASTQ, GTF and FASTA files)"
	@echo "# make index  (generates the hisat2 index)"
	@echo "# make design (generates the ${DESIGN} file)"
	@echo "# make align  (aligns reads into BAM files)"
	@echo "# make count  (counts the reads in the BAM files)"
	@echo "#"

	@echo "# make all    (runs all steps)"
	@echo "#"

# Generate the design file.
design: 
	bio search ${PRJ} -H --csv > ${DESIGN}

# Download the sequencing data and references.
data: 
	curl -L "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000013425.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" -o ncbi_dataset.zip
unzip:
	unzip ncbi_dataset.zip -d ncbi_dataset

# Generate the HISAT2 index
index:
	make -f src/run/hisat2.mk index REF=${REF}

# Run the alignment
align: ${DESIGN}
	# Runs the alignment in parallel over all samples.
	cat ${DESIGN} | parallel ${FLAGS} \
		hisat2-build refs/staph-2006.fa refs/staph-2006.idx \
		NCPU=${NCPU} \
		REF=${REF} \
		R1=reads/{sample}_R1.fq \
		BAM=bam/{sample}.bam \
		run

# The counts as textfile produced by featurecounts.
${COUNTS_TXT}:
	# Make the directory name for the counts
	mkdir -p $(dir $@)

	# Count the features
	cat ${DESIGN} | \
		parallel --header : --colsep , -k echo bam/{sample}.bam | \
		parallel -u --xargs featureCounts -a ${GTF} -o ${COUNTS_TXT} {}

# The final counts in CSV format.
${COUNTS}: ${COUNTS_TXT}
	micromamba run -n stats Rscript src/r/format_featurecounts.r -c ${COUNTS_TXT} -o ${COUNTS}

# Trigger the counting explicitly
count: ${COUNTS}
	@ls -lh ${COUNTS_TXT}
	@ls -lh ${COUNTS}

all: data index align count

.PHONY: usage design data index align count all

