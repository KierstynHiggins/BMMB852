### BMMB 852: Variant Calling
####  Kierstyn Higgins
##### 11/10/24

We will generate a BAM file for the organism *Staphylococcus aurelis* and generate a variant call file.

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
all: bam vcf
```
`clean`: remove generate files.
```bash
clean:
	rm -rf ${REF} ${GFF} ${R1} ${R2} ${BAM} ${VCF}
```
We can then visualize the variant call file in IGV.
<img width="1143" alt="Screenshot 2024-11-10 at 5 06 02 PM" src="https://github.com/user-attachments/assets/0a8c4fb5-7baa-4696-8143-60fcc7b50fd9">

To investigate the vcf further, we can run the following code...
```bash
bcftools stats vcf/SAMEA115926234.vcf.gz
```
This will output:
```bash
bcftools stats vcf/SAMEA115926234.vcf.gz
# This file was produced by bcftools stats (1.21+htslib-1.21) and can be plotted using plot-vcfstats.
# The command line was:	bcftools stats  vcf/SAMEA115926234.vcf.gz
#
# Definition of sets:
# ID	[2]id	[3]tab-separated file names
ID	0	vcf/SAMEA115926234.vcf.gz
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN	[2]id	[3]key	[4]value
SN	0	number of samples:	1
SN	0	number of records:	1113
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	1108
SN	0	number of MNPs:	0
SN	0	number of indels:	5
SN	0	number of others:	0
SN	0	number of multiallelic sites:	0
SN	0	number of multiallelic SNP sites:	0
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV	[2]id	[3]ts	[4]tv	[5]ts/tv	[6]ts (1st ALT)	[7]tv (1st ALT)[8]ts/tv (1st ALT)
TSTV	0	70	1038	0.07	70	1038	0.07
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS	[2]id	[3]allele count	[4]number of SNPs	[5]number of transition[6]number of transversions	[7]number of indels	[8]repeat-consistent	[9]repeat-inconsistent	[10]not applicable
SiS	0	1	64	13	51	1	0	0	1
# AF, Stats by non-reference allele frequency:
# AF	[2]id	[3]allele frequency	[4]number of SNPs	[5]number of transitions	[6]number of transversions	[7]number of indels	[8]repeat-consistent	[9]repeat-inconsistent	[10]not applicable
AF	0	0.000000	64	13	51	1	0	0	1
AF	0	0.990000	1044	57	987	4	0	0	4
# QUAL, Stats by quality
# QUAL	[2]id	[3]Quality	[4]number of SNPs	[5]number of transitions (1st ALT)	[6]number of transversions (1st ALT)	[7]number of indels
QUAL	0	3.0	1	0	1	0
QUAL	0	3.2	5	0	5	0
QUAL	0	3.6	5	0	5	0
QUAL	0	3.7	4	1	3	0
QUAL	0	4.0	1	1	0	0
QUAL	0	4.3	1	0	1	0
QUAL	0	5.0	1	0	1	0
QUAL	0	5.6	9	0	9	0
QUAL	0	5.7	54	3	51	0
QUAL	0	6.0	1	0	1	0
QUAL	0	6.4	1	0	1	0
QUAL	0	7.2	5	0	5	0
QUAL	0	7.3	6	0	6	0
QUAL	0	7.4	1	0	1	0
QUAL	0	8.1	916	41	875	0
QUAL	0	8.8	1	0	1	0
QUAL	0	8.9	1	0	1	0
QUAL	0	9.2	3	0	3	0
QUAL	0	9.5	1	0	1	0
QUAL	0	10.7	1	1	0	0
QUAL	0	13.5	2	0	2	0
QUAL	0	13.6	1	0	1	0
QUAL	0	15.5	4	0	4	0
QUAL	0	16.0	2	0	2	0
QUAL	0	16.2	2	1	1	0
QUAL	0	16.5	22	1	21	0
QUAL	0	17.4	0	0	0	1
QUAL	0	18.4	5	0	5	0
QUAL	0	20.1	1	0	1	0
QUAL	0	20.2	1	1	0	0
QUAL	0	20.7	2	1	1	0
QUAL	0	21.4	1	0	1	0
QUAL	0	23.0	1	0	1	0
QUAL	0	23.2	2	0	2	0
QUAL	0	23.3	1	0	1	0
QUAL	0	23.4	2	1	1	0
QUAL	0	24.6	0	0	0	1
QUAL	0	25.6	1	1	0	0
QUAL	0	25.9	3	2	1	0
QUAL	0	26.3	1	1	0	0
QUAL	0	27.4	6	2	4	0
QUAL	0	30.4	1	1	0	0
QUAL	0	32.4	2	0	2	1
QUAL	0	33.9	1	1	0	0
QUAL	0	37.4	7	2	5	0
QUAL	0	40.4	1	1	0	0
QUAL	0	44.4	2	1	1	0
QUAL	0	49.4	1	0	1	0
QUAL	0	54.4	2	0	2	0
QUAL	0	57.7	1	0	1	0
QUAL	0	59.2	1	1	0	0
QUAL	0	62.4	1	0	1	0
QUAL	0	64.4	1	0	1	0
QUAL	0	67.4	1	1	0	0
QUAL	0	68.3	0	0	0	1
QUAL	0	71.3	1	0	1	0
QUAL	0	75.4	5	2	3	0
QUAL	0	79.4	0	0	0	1
QUAL	0	83.4	2	2	0	0
# IDD, InDel distribution:
# IDD	[2]id	[3]length (deletions negative)	[4]number of sites	[5]number of genotypes	[6]mean VAF
IDD	0	-1	2	0	.
IDD	0	1	2	0	.
IDD	0	2	1	0	.
# ST, Substitution types:
# ST	[2]id	[3]type	[4]count
ST	0	A>C	111
ST	0	A>G	22
ST	0	A>T	105
ST	0	C>A	286
ST	0	C>G	20
ST	0	C>T	21
ST	0	G>A	12
ST	0	G>C	18
ST	0	G>T	290
ST	0	T>A	112
ST	0	T>C	15
ST	0	T>G	96
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP	[2]id	[3]bin	[4]number of genotypes	[5]fraction of genotypes (%)	[6]number of sites	[7]fraction of sites (%)
DP	0	1	0	0.000000	920	82.659479
DP	0	2	0	0.000000	104	9.344115
DP	0	3	0	0.000000	53	4.761905
DP	0	4	0	0.000000	25	2.246181
DP	0	5	0	0.000000	9	0.808625
DP	0	6	0	0.000000	1	0.089847
DP	0	7	0	0.000000	1	0.089847
```
This tells us that the VCF file contains data for a single sample, with a total of 1113 variants. The majority of these variants are SNPs (1108), with a very small number of indels (5). There are no multiallelic sites or MNPs.


