# Makefile for Genome and SRA Processing

# Variables
GENOME = GCF_000013425.1
FILE = $(GENOME)_ASM1342V1_genomic.fna
LENGTH = 150
COVERAGE = 10
N = 200000
SRA = SRR5837600
REF = brachy-2017.fa
R1 = $(SRA)_1.fastq
R2 = $(SRA)_2.fastq
R1_CUT = $(SRA)_1.trimmed.fastq
R2_CUT = $(SRA)_2.trimmed.fastq
BAM = align.bam
SAM = align.sam
R1sim = simulated_reads_1.fq
R2sim = simulated_reads_1.fq
SAMsim = align_sim.sam
BAMsim = align_sim.bam
#------------------------------------------------------------------------------
#
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# Target: usage
usage:
	@echo "Usage:"
	@echo "  make genome      - Download the genome"
	@echo "  make simulate    - Simulate reads for the genome"
	@echo "  make download    - Download reads from SRA"
	@echo "  make trim        - Trim SRA reads"
	@echo "  make clean       - Remove all generated files"
	@echo "  make index       - Index reference genome
	@echo "  make align       - Align reads to ref genome
	@echo "  make BAM         - Convert SAM to BAM
	@echo "  make stat        - Generate alignment stats

# Target: genome
genome:
	@echo "Downloading genome..."
	datasets download genome accession $(GENOME) --include gff3,cds,protein,rna,genome && \
	unzip -o ncbi_dataset.zip && \
	ln -sf ncbi_dataset/data/$(GENOME)/$(FILE) genome.fa
	@echo "Genome downloaded as genome.fa"
	@echo "File size: $$(wc -c < genome.fa) bytes"
	@echo "Total genome size: $$(grep -v '^>' genome.fa | tr -d '\n' | wc -c) bases"
	@echo "Number of chromosomes: $$(grep -c '^>' genome.fa)"
	@echo "Chromosome lengths:"
	@grep '^>' genome.fa | awk '{name=$$1; getline; print name, length($$0)}' | sed 's/>//g'

# Target: simulate
simulate: 
	@echo "Simulated FASTQ files: simulated_reads_1.fq and simulated_reads_2.fq"
	wgsim -N ${N} -1 ${LENGTH} -2 ${LENGTH} -e 0 -r 0 -R 0 -X 0 genome.fa ${R1sim} ${R2sim}
	seqkit stats ${R1sim} ${R2sim}
# Target: download
download:
	@echo "Downloading reads from SRA using fastq-dump..."
	fastq-dump -X 10000 --split-files $(SRA)

# Target: trim
trim: 
	@echo "Trimming FASTQ files..."
	fastqc $(R1) $(R2)
	fastp --cut_tail -i $(R1) -o $(R1_CUT) -I $(R2) -O $(R2_CUT)
	fastqc $(R1_CUT)


# Target: index
index:
	@echo "Indexing the reference genome..."
	bwa index ${REF} 

# Target: align
align:
	@echo "Aligning reads to the reference genome..."
	bwa mem ${REF} ${R1} ${R2} > ${SAM}
	@echo "Aligning simulated reads to reference genome..."
	bwa mem ${REF} ${R1sim} ${R2sim} > ${SAMsim}
# Target: BAM
BAM:
	@echo "Converting SAM to BAM..."
	cat ${SAM} | samtools sort > ${BAM}
	samtools index ${BAM}
	@echo "Converting SAMsim to BAMsim..."
	cat ${SAMsim} | samtools sort > ${BAMsim}
	samtools index ${BAMsim}
# Target: stat
stat:
	@echo "Generating alignment statistics... "
	samtools flagstat ${BAM}
	samtools flagstat ${BAMsim}

