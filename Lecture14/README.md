## BMMB 852: RNA Seq Assignment
### Kierstyn Higgins
### 12/08/2024
#### *Note: I really struggled with this assignment and spent way more hours than I'd like to admit, so if you have any feedback please let me know!

First we need to set the variables:
```bash
# The URL for the data.
DATA_URL = https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download?filename=ncbi_dataset.zip&ncbi_phid=322CB9919F68EEE500003635566DC9E4.1.m_3.020

# The downloaded data file name.
DATA_FILE = $(notdir ${DATA_URL})

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

# Directories
BAM_DIR=bam
SORTED_BAM_DIR=sorted_bam
INDEX_DIR=index
```
`usage`: will print usage for each target.
```bash
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
```

`design`: This generates the design file from the project number.
(I had no issues with this step)
```bash
design: 
	bio search ${PRJ} -H --csv > ${DESIGN}
```
`data`: Download the data.
(I had to use curl instead of wget because it didn't like wget for some reason)
```
data: 
	curl -L ${DATA_URL}
```
`unzip`: Unzip the data.
```
unzip ncbi_dataset.zip -d ncbi_dataset
```
`index`: Generate HISAT2 index. 


First I tried:
```bash
make -f src/run/hisat2.mk index REF=${REF}
```
But that gave me this error:
```bash
ake -f src/run/hisat2.mk index REF=ncbi_dataset/refs/staph-2006.fa
# file not found: REF=ncbi_dataset/refs/staph-2006.fa
make[1]: *** [src/run/hisat2.mk:79: ncbi_dataset/refs//idx/staph-2006.fa.1.ht2] Error 255
make: *** [Makefile:85: index] Error 2
(bioinfo) 
```

Then I tried something else which seemed to work?
```bash
hisat2-build refs/staph-2006.fa refs/staph-2006.idx
```

`align`: Run the alignment. (I think this worked but I got an insanely long output.)
```bash
align: ${DESIGN}
	# Runs the alignment in parallel over all samples and converts SAM to BAM, sorts, and indexes.
	cat ${DESIGN} | parallel ${FLAGS} \
		hisat2 -x ${REF} -1 reads/{sample}_R1.fq -2 reads/{sample}_R2.fq -S ${BAM_DIR}/{sample}.sam \
		&& samtools view -bS ${BAM_DIR}/{sample}.sam > ${BAM_DIR}/{sample}.bam \
		&& samtools sort ${BAM_DIR}/{sample}.bam -o ${SORTED_BAM_DIR}/{sample}_sorted.bam \
		&& samtools index ${SORTED_BAM_DIR}/{sample}_sorted.bam
```
`count`: Once we got to the counts is where I ran into trouble...
```bash
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
```
This would output...
```bash
# Count the features
cat design.csv | \
	parallel --header : --colsep , -k echo bam/{sample}.bam | \
	parallel -u --xargs featureCounts -a ncbi_dataset/refs/staph-2006.gtf -o res/counts-hisat.txt {}
featureCounts(9504,0x204a99200) malloc: Incorrect checksum for freed object 0x7fa159704b90: probably modified after being freed.
Corrupt value: 0xf0003b7365736162
featureCounts(9504,0x204a99200) malloc: *** set a breakpoint in malloc_error_break to debug
make: *** [Makefile:100: res/counts-hisat.txt] Error 1
(bioinfo)
```
At this point I was out of solutions, so if anyone reading this has one please share!
