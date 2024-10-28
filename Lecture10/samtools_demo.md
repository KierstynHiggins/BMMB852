# BMMB 852: Makefiles, BAM, and samtools
## Kierstyn Higgins
## 10/27/2024

**Variables**: First we will set the variables. These can be interchanged with the accession numbers of your choice.
```bash
GENOME = GCF_000013425.1
FILE = $(GENOME)_ASM1342V1_genomic.fna
LENGTH = 150
COVERAGE = 10
N = 200000
SRA = ERR13649352
REF = genome.fa
R1 = $(SRA)_1.fastq
R2 = $(SRA)_2.fastq
R1_CUT = $(SRA)_1.trimmed.fastq
R2_CUT = $(SRA)_2.trimmed.fastq
BAM = align.bam
FILTERED_BAM = filtered.bam
SAM = align.sam
R1sim = simulated_reads_1.fq
R2sim = simulated_reads_1.fq
SAMsim = align_sim.sam
BAMsim = align_sim.bam

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
```
`Usage`: This target displays the available commands and their descriptions for the user. It outputs a help message, listing all the other targets that can be executed with the Makefile, along with a brief description of what each command does.
```bash
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
	@echo "  make BAM_sim     - Converts SAMsim to BAMsim
	@echo "  make stat        - Generate alignment stats
	@echo "  make counts      - Counts reads for various specifications
	@echo "  make filter      - Creates filtered BAM file
	@echo "  make flagstats   - Gives stats for original and filtered BAM file
```
`Genome`: This target downloads the specified genome and prepares the genomic file for use. It uses the datasets download command to download the genome associated with the specified accession number. It then unzips the downloaded dataset and creates a symbolic link to the genomic FASTA file, renaming it to genome.fa. It also outputs the size of the file, the total number of bases in the genome, the number of chromosomes, and the lengths of each chromosome.
```bash
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
simulate: 
	@echo "Simulated FASTQ files: simulated_reads_1.fq and simulated_reads_2.fq"
	wgsim -N ${N} -1 ${LENGTH} -2 ${LENGTH} -e 0 -r 0 -R 0 -X 0 genome.fa ${R1sim} ${R2sim}
	seqkit stats ${R1sim} ${R2sim}
```
`download`: This target downloads reads from the SRA based on the specified accession. It uses the fastq-dump command to download a specified number of reads (10,000 in this case) and splits them into two files (SRR5837600_1.fastq and SRR5837600_2.fastq), which correspond to paired-end sequencing data.
```bash
download:
	@echo "Downloading reads from SRA using fastq-dump..."
	fastq-dump -X 10000 --split-files $(SRA)
```
`trim`: This trims the downloaded SRA reads for quality control. It depends on the download target to ensure that the reads are downloaded first. It uses fastqc to assess the quality of the raw FASTQ files and then employs fastp to trim the reads, producing two new trimmed FASTQ files (R1_cut and R2_cut). It performs a final quality check on the trimmed files.
```bash
trim: download
	@echo "Trimming FASTQ files..."
	fastqc $(R1) $(R2)
	fastp --cut_tail -i $(R1) -o $(R1_CUT) -I $(R2) -O $(R2_CUT)
	fastqc $(R1_CUT)
```
`index`: This indexes the reference genome using bwa.
```bash
index:	
	@echo "Indexing the reference genome..."
	bwa index ${REF}
```
`align`: This aligns reads to the reference genome using `bwa mem`.
```bash
align:
	@echo "Aligning reads to the reference genome..."
	bwa mem ${REF} ${R1} ${R2} > ${SAM}
	@echo "Aligning simulated reads to reference genome..."
	bwa mem ${REF} ${R1sim} ${R2sim} > ${SAMsim}
```
`BAM`: This generates a BAM file from the SAM file.
```bash
BAM:
	@echo "Converting SAM to BAM..."
	cat ${SAM} | samtools sort > ${BAM}
	samtools index ${BAM}
```
You can then visualize the BAM file in IGV. 
`BAM_sim`: This generates a BAM file from SAMsim.
```bash
	@echo "Converting SAMsim to BAMsim..."
	cat ${SAMsim} | samtools sort > ${BAMsim}
	samtools index ${BAMsim}
```

`stat`: This will generate alignment statistics.
```bash
stat:
	@echo "Generating alignment statistics... "
	samtools flagstat ${BAM}
	samtools flagstat ${BAMsim}
```
It will output...
```bash
$ make stat
Generating alignment statistics... 
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
404482 + 0 in total (QC-passed reads + QC-failed reads)
400000 + 0 primary
0 + 0 secondary
4482 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
404482 + 0 mapped (100.00% : N/A)
400000 + 0 primary mapped (100.00% : N/A)
400000 + 0 paired in sequencing
200000 + 0 read1
200000 + 0 read2
0 + 0 properly paired (0.00% : N/A)
400000 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
The first summary shows that only a small number of reads (8) are actually mapped, indicating that the original BAM file might have had significant issues with read quality or alignment.The second summary shows a much larger number of reads (404,482) and indicates that they are all mapped, suggesting this might be from a better-aligned dataset or after successful filtering.

I tried using a different reference genome, a different SRA number, different fastqc methods and nothing changed. I really can't figure out why.

`counts`: This target will count the alignments that we selected for.
```bash
counts:
	@echo "Counting reads that did not align..."
	@echo "Number of unmapped reads: $$(samtools view -c -f 4 $(BAM))"
	@echo "Primary alignments: $$(samtools view -c -F 256 $(BAM))"
	@echo "Secondary alignments: $$(samtools view -c -f 256 $(BAM))"
	@echo "Supplementary alignments: $$(samtools view -c -f 2048 $(BAM))"
	@echo "Properly paired alignments on the reverse strand (read1): $$(samtools view -c -f 2 -f 16 $(BAM))"
```
This will output...
```bash
$ make counts
Counting reads that did not align...
Number of unmapped reads: 19996
Primary alignments: 20004
Secondary alignments: 0
Supplementary alignments: 4
Properly paired alignments on the reverse strand (read1): 0
```
`filter`: This will create the filtered BAM file with properly paired alignments.
```bash
filter:
	@echo "Creating filtered BAM file with properly paired primary alignments (MQ > 10)..."
	samtools view -b -q 10 -f 2 -F 256 $(BAM) > $(FILTERED_BAM)
```
`flagstats`: This will compare stats for the filtered and unfiltered BAM files.
```bash
flagstats: filter
	@echo "Flagstats for original BAM file:"
	samtools flagstat $(BAM)
	@echo "Flagstats for filtered BAM file:"
	samtools flagstat $(FILTERED_BAM)
```
This will output....
```bash
$ make flagstats
Creating filtered BAM file with properly paired primary alignments (MQ > 10)...
Flagstats for original BAM file:
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
Flagstats for filtered BAM file:
0 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
0 + 0 mapped (N/A : N/A)
0 + 0 primary mapped (N/A : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
Here I tried changing the MQ threshold but that still didn't change the results.
