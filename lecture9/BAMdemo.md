# BMMB 852: Makefiles and BAM files
## Kierstyn Higgins

**Variables**: First we will set the variables. These can be interchanged with the accession numbers of your choice.
```bash
GENOME = GCF_000013425.1
FILE = $(GENOME)_ASM1342V1_genomic.fna
LENGTH = 150
COVERAGE = 10
SRA = SRR5837600
REF = brachy-2017.fa
R1 = $(SRA)_1.fastq
R2 = $(SRA)_2.fastq
R1_CUT = $(SRA)_1.trimmed.fastq
R2_CUT = $(SRA)_2.trimmed.fastq
BAM = align.bam
SAM = align.sam
```
`Usage`: This target displays the available commands and their descriptions for the user. It outputs a help message, listing all the other targets that can be executed with the Makefile, along with a brief description of what each command does.
```bash
.PHONY: usage
usage:
	@echo "Usage:"
	@echo "  make genome      - Download the genome"
	@echo "  make simulate    - Simulate reads for the genome"
	@echo "  make download     - Download reads from SRA"
	@echo "  make trim        - Trim SRA reads"
	@echo "  make clean       - Remove all generated files"
	@echo "  make index       - Index reference genome
	@echo "  make align       - Align reads to ref genome
	@echo "  make BAM         - Convert SAM to BAM
	@echo "  make stat        - Generate alignment stats
```
`Genome`: This target downloads the specified genome and prepares the genomic file for use. It uses the datasets download command to download the genome associated with the specified accession number. It then unzips the downloaded dataset and creates a symbolic link to the genomic FASTA file, renaming it to genome.fa. It also outputs the size of the file, the total number of bases in the genome, the number of chromosomes, and the lengths of each chromosome.
```bash
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
```
`simulate`: This target simulates reads based on the downloaded genome. It uses the genome target to ensure the genome is downloaded and prepared. It calculates the number of reads needed based on the genome length, desired coverage, and read length. It uses wgsim to simulate paired-end reads and generates two FASTQ files (simulated_reads_1.fq and simulated_reads_2.fq). Finally, it compresses these files using gzip.
```bash
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
```
`download`: This target downloads reads from the SRA based on the specified accession. It uses the fastq-dump command to download a specified number of reads (10,000 in this case) and splits them into two files (SRR5837600_1.fastq and SRR5837600_2.fastq), which correspond to paired-end sequencing data.
```bash
.PHONY: download
download:
	@echo "Downloading reads from SRA using fastq-dump..."
	fastq-dump -X 10000 --split-files $(SRA)
```
`trim`: This trims the downloaded SRA reads for quality control. It depends on the download target to ensure that the reads are downloaded first. It uses fastqc to assess the quality of the raw FASTQ files and then employs fastp to trim the reads, producing two new trimmed FASTQ files (R1_cut and R2_cut). It performs a final quality check on the trimmed files.
```bash
.PHONY: trim
trim: download
	@echo "Trimming FASTQ files..."
	fastqc $(R1) $(R2)
	fastp --cut_tail -i $(R1) -o $(R1_CUT) -I $(R2) -O $(R2_CUT)
	fastqc $(R1_CUT)
```
`clean`: This cleans up all generated files and restores the working directory. This removes all files generated during the execution of the other targets, including the downloaded genome, simulated reads, and trimmed FASTQ files.
```bash
.PHONY: clean
clean:
	@echo "Cleaning up..."
	rm -f genome.fa simulated_reads_1.fq.gz simulated_reads_2.fq.gz ncbi_dataset.zip $(R1) $(R2) $(R1_CUT) $(R2_CUT)
```
`index`: This indexes the reference genome using bwa.
```bash
.PHONY: index
index:	
	@echo "Indexing the reference genome..."
	bwa index ${REF}
```
`align`: This aligns reads to the reference genome using `bwa mem`.
```bash
.PHONY: align
align:
	@echo "Aligning reads to the reference genome..."
	bwa mem ${REF} ${R1} ${R2} > ${SAM}
```
`BAM`: This generates a BAM file from the SAM file.
```bash
.PHONY: BAM
BAM:
	@echo "Converting SAM to BAM..."
	cat ${SAM} | samtools sort > ${BAM}
	samtools index ${BAM}
```
You can then visualize the BAM file in IGV. 

`stat`: This will generate alignment statistics.
```bash
.PHONY: stat
stat:
	@echo "Generating alignment statistics... "
	samtools flagstat ${BAM}
```
It will output...
```bash
$ samtools flagstat align.bam
20004 + 0 in total (QC-passed reads + QC-failed reads)
20000 + 0 primary
0 + 0 secondary
4 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
8 + 0 mapped (0.04% : N/A)
4 + 0 primary mapped (0.02% : N/A)
20000 + 0 paired in sequencing
10000 + 0 read1
10000 + 0 read2
0 + 0 properly paired (0.00% : N/A)
4 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
