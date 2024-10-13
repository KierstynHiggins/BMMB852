# BMMB 852: Makefile for Genome and SRA Processing Demo

# Variables
GENOME = GCF_000013425.1
FILE = $(GENOME)_ASM1342V1_genomic.fna
LENGTH = 150
COVERAGE = 10
SRA = SRR5837600
R1 = $(SRA)_1.fastq
R2 = $(SRA)_2.fastq
R1_CUT = $(SRA)_1.trimmed.fastq
R2_CUT = $(SRA)_2.trimmed.fastq

# Target: usage
.PHONY: usage
usage:
	@echo "Usage:"
	@echo "  make genome      - Download the genome"
	@echo "  make simulate    - Simulate reads for the genome"
	@echo "  make download     - Download reads from SRA"
	@echo "  make trim        - Trim SRA reads"
	@echo "  make clean       - Remove all generated files"

# Target: genome
.PHONY: genome
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
.PHONY: simulate
simulate: genome
	@echo "Simulating reads for the genome..."
	GENOME_LENGTH=$$(grep -v '^>' genome.fa | tr -d '\n' | wc -c); \
	NUM_READS=$$(((COVERAGE * GENOME_LENGTH) / $(LENGTH))); \
	wgsim -N $$NUM_READS -1 $(LENGTH) -2 $(LENGTH) -e 0.02 -r 0.05 -R 0.15 -X 0.3 genome.fa simulated_reads_1.fq simulated_reads_2.fq; \
	@echo "Simulated FASTQ files: simulated_reads_1.fq and simulated_reads_2.fq"; \
	@echo "Number of reads generated: $$NUM_READS"; \
	gzip simulated_reads_1.fq simulated_reads_2.fq; \
	@echo "Files compressed to: simulated_reads_1.fq.gz, simulated_reads_2.fq.gz"

# Target: download
.PHONY: download
download:
	@echo "Downloading reads from SRA using fastq-dump..."
	fastq-dump -X 10000 --split-files $(SRA)

# Target: trim
.PHONY: trim
trim: download
	@echo "Trimming FASTQ files..."
	fastqc $(R1) $(R2)
	fastp --cut_tail -i $(R1) -o $(R1_CUT) -I $(R2) -O $(R2_CUT)
	fastqc $(R1_CUT)

# Target: clean
.PHONY: clean
clean:
	@echo "Cleaning up..."
	rm -f genome.fa simulated_reads_1.fq.gz simulated_reads_2.fq.gz ncbi_dataset.zip $(R1) $(R2) $(R1_CUT) $(R2_CUT)

