### BMMB 852: Variant Effect Prediction
####  Kierstyn Higgins
##### 11/17/24

We will generate a BAM file for the organism *Staphylococcus aurelis* and generate a variant call file. Then we will generate a variant effect prediction report using the tool snpEff.

First we need to set the variables.
```bash
# Accession number of the Staphylococcus aurelis genome.
ACC=GCF_000013425.1

# The reference file.
REF=refs/staph-2006.fa

# The GFF file.
GFF=refs/staph-2006.gff

# The sequencing read accession number.
SRR=ERR13649352

# The number of reads to get
N=5000

# The name of the sample (see: bio search ERR13649352)
SAMPLE=SAMEA115926234

# The path to read 1
R1=reads/${SAMPLE}_1.fastq

# The path to read 2
R2=reads/${SAMPLE}_2.fastq

# The resulting BAM file.
BAM=bam/${SAMPLE}.bam

# The resulting variant VCF file (compressed!).
VCF=vcf/${SAMPLE}.vcf.gz
```
We also need to set the makefile settings.
```bash
# Custom makefile settings.
SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
```
`usage`: prints the usage of the makefile.
```bash
usage:
	@echo "#"
	@echo "# SNP call demonstration"
	@echo "#"
	@echo "# ACC=${ACC}"
	@echo "# SRR=${SRR}"
	@echo "# SAMPLE=${SAMPLE}"
	@echo "# BAM=${BAM}"
	@echo "# VCF=${VCF}"
	@echo "#"
	@echo "# make bam|vcf|all"
	@echo "#"
```
`bam`: create the BAM alignment file.
```bash
bam: ${CHECK_FILE}
	# Get the reference genome and the annotations.
	make -f src/run/genbank.mk ACC=${ACC} REF=${REF} GFF=${GFF} fasta gff

	# Index the reference genome.
	make -f src/run/bwa.mk REF=${REF} index

	# Download the sequence data.
	make -f src/run/sra.mk SRR=${SRR} R1=${R1} R2=${R2} N=${N} run

	# Align the reads to the reference genome. 
	# Use a sample name in the readgroup.
	make -f src/run/bwa.mk SM=${SAMPLE} REF=${REF} R1=${R1} R2=${R2} BAM=${BAM} run stats
```
`vcf`: call the SNPs in the BAM file.
```bash
vcf:
	make -f src/run/bcftools.mk REF=${REF} BAM=${BAM} VCF=${VCF} run
```
`all`: run all steps.
```bash
all: bam vcf snpeff
```
`clean`: remove generate files.
```bash
clean:
	rm -rf ${REF} ${GFF} ${R1} ${R2} ${BAM} ${VCF}
```
`snpeff`: we can run the variant effect prediction target using snpEff.
```bash
snpeff:
	# Build the custom database
	make -f src/run/snpeff.mk GFF=${GFF} REF=${REF} build

	# Run snpEff using the custom database
	make -f src/run/snpeff.mk GFF=${GFF} REF=${REF} VCF=${VCF} run
# These targets do not correspond to files.
.PHONY: bam vcf all usage clean snpeff	
```
This will generate an html report as shown. We can see...
  - Number of variants by type, where INS = insertion and DEL = deletion.
   <img width="224" alt="Screenshot 2024-11-18 at 11 58 51 AM" src="https://github.com/user-attachments/assets/297e96f4-91b6-466d-a119-d2d7a4ebfc74">

   
  - Number of effects by impact.
   <img width="380" alt="Screenshot 2024-11-18 at 12 02 32 PM" src="https://github.com/user-attachments/assets/21fcdd45-fa1e-4c4d-b9ec-73c6442ff0c3">

   
  - Number of effects by functional class.
    
   <img width="360" alt="Screenshot 2024-11-18 at 12 05 05 PM" src="https://github.com/user-attachments/assets/762ba687-c5ef-42f2-8e54-8702c35d504c">
    
  Missense: a point mutation (single nucleotide change) in the DNA that results in a change in the amino acid sequence of the protein.
  
  Nonsense: a point mutation that changes a codon encoding an amino acid into a stop codon.
  
  Silent: a point mutation in the DNA that does not result in a change to the amino acid sequence of the protein.

  - And lastly, the number of effects by type and region. Here we see most of the effects are in downstream and upstream gene variants. These variants may not directly affect the protein sequence but can still have significant effects on gene expression and regulatory processes.
    
   <img width="794" alt="Screenshot 2024-11-18 at 12 11 18 PM" src="https://github.com/user-attachments/assets/a5d4ece4-626c-45dd-af57-07a05111fc5c">

