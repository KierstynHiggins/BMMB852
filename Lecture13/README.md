# BMMB 852: How to Automate a VCF
## Kierstyn Higgins
### 12/1/2024
___________________________________________________________________________________
** Make sure the toolbox is installed using:
```bash
bio code
```
1. Download the Automation.mk Makefile to your directory.
```bash
# Define the GitHub repository and the path to the Makefile
REPO="https://github.com/KierstynHiggins/BMMB852/blob/main/Lecture13/Automation.mk"

# Download the Makefile using wget
wget $REPO
```
2. Run the Makefile to make sure it works.
```bash
make all
```
3. Now we need a list of samples. I used the samples associated with the project PRJEB75336- "WGS of S. aureus isolates serially passaged with either IMP-1700 or zoliflodicin".
   You can download them to a csv using the following code.
```bash
bio search PRJEB75336 -H --csv > design.csv
```
4. This code will run the Makefile with the list of samples.
```bash
cat design.csv | head -25 | \
  parallel --lb -j 4 --colsep , --header : \
  make all SRR={run_accession} SAMPLE={sample_alias}
```
5. Now we can merge the vcf files that we generated into one.
```bash
# Merge VCF files into a single one.
bcftools merge -0 vcf/*.vcf.gz -O z > merged.vcf.gz

# Index the merged VCF file
bcftools index merged.vcf.gz
```
6. Create a .tbi index file so that it will be compatible with IGV software.
```bash
tabix -p vcf merged.vcf.gz
```
When we load this merged VCF into IGV, it will look like this...
<img width="1136" alt="Screenshot 2024-12-01 at 9 57 52 PM" src="https://github.com/user-attachments/assets/26320c0b-ef00-4705-93ab-2cd2c8a18028">

We can also convert stats from the merged VCF to a .txt file using the following code:
```bash
bcftools stats merged.vcf > vcf_stats.txt
```
This will output something like this:
```bash
# This file was produced by bcftools stats (1.21+htslib-1.21) and can be plotted using plot-vcfstats.
# The command line was:	bcftools stats  merged.vcf.gz
#
# Definition of sets:
# ID	[2]id	[3]tab-separated file names
ID	0	merged.vcf.gz
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
SN	0	number of samples:	12
SN	0	number of records:	10333
SN	0	number of no-ALTs:	0
SN	0	number of SNPs:	10258
SN	0	number of MNPs:	0
SN	0	number of indels:	75
SN	0	number of others:	0
SN	0	number of multiallelic sites:	19
SN	0	number of multiallelic SNP sites:	19
# TSTV, transitions/transversions
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)
#   - transversions, see https://en.wikipedia.org/wiki/Transversion
# TSTV	[2]id	[3]ts	[4]tv	[5]ts/tv	[6]ts (1st ALT)	[7]tv (1st ALT)	[8]ts/tv (1st ALT)
TSTV	0	2484	7793	0.32	2473	7785	0.32
# SiS, Singleton stats:
#   - allele count, i.e. the number of singleton genotypes (AC=1)
#   - number of transitions, see above
#   - number of transversions, see above
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]
# SiS	[2]id	[3]allele count	[4]number of SNPs	[5]number of transitions	[6]number of transversions	[7]number of indels	[8]repeat-consistent	[9]repeat-inconsistent	[10]not applicable
SiS	0	1	455	137	318	11	0	0	11
# AF, Stats by non-reference allele frequency:
# AF	[2]id	[3]allele frequency	[4]number of SNPs	[5]number of transitions	[6]number of transversions	[7]number of indels	[8]repeat-consistent	[9]repeat-inconsistent	[10]not applicable
AF	0	0.000000	455	137	318	11	0	0	11
AF	0	0.080000	9163	1935	7228	52	0	0	52
AF	0	0.120000	54	30	24	1	0	0	1
AF	0	0.160000	448	289	159	5	0	0	5
AF	0	0.200000	13	6	7	0	0	0	0
AF	0	0.240000	103	67	36	5	0	0	5
AF	0	0.280000	1	1	0	0	0	0	0
AF	0	0.330000	17	10	7	1	0	0	1
AF	0	0.410000	9	5	4	0	0	0	0
AF	0	0.490000	8	3	5	0	0	0	0
AF	0	0.530000	1	1	0	0	0	0	0
AF	0	0.570000	4	0	4	0	0	0	0
AF	0	0.660000	1	0	1	0	0	0	0
# QUAL, Stats by quality
# QUAL	[2]id	[3]Quality	[4]number of SNPs	[5]number of transitions (1st ALT)	[6]number of transversions (1st ALT)	[7]number of indels
QUAL	0	3.0	2	0	2	0
QUAL	0	3.1	14	3	11	0
QUAL	0	3.2	83	28	55	2
QUAL	0	3.6	12	0	12	0
QUAL	0	3.7	26	5	21	0
QUAL	0	3.9	1	0	1	0
QUAL	0	4.0	1	1	0	0
QUAL	0	4.2	18	2	16	0
QUAL	0	4.3	3991	907	3084	0
QUAL	0	4.4	1	0	1	0
QUAL	0	4.9	5	0	5	0
QUAL	0	5.0	17	7	10	2
QUAL	0	5.3	0	0	0	1
QUAL	0	5.4	6	2	4	0
QUAL	0	5.5	6	2	4	0
QUAL	0	5.6	13	0	13	0
QUAL	0	5.7	170	15	155	0
QUAL	0	6.0	2	0	2	0
QUAL	0	6.3	2	1	1	0
QUAL	0	6.4	11	5	6	1
QUAL	0	6.5	3	2	1	0
QUAL	0	7.2	12	1	11	0
QUAL	0	7.3	61	11	50	1
QUAL	0	7.4	3	0	3	0
QUAL	0	7.5	5	3	2	0
QUAL	0	7.7	2	2	0	0
QUAL	0	7.8	1	1	0	0
QUAL	0	8.1	3601	368	3233	1
QUAL	0	8.2	3	3	0	0
QUAL	0	8.6	1	0	1	0
QUAL	0	8.8	3	1	2	0
QUAL	0	8.9	8	3	5	0
QUAL	0	9.0	1	0	1	0
QUAL	0	9.2	12	0	12	0
QUAL	0	9.5	1	0	1	0
QUAL	0	9.6	2	2	0	0
QUAL	0	9.8	23	7	16	1
QUAL	0	10.6	2	0	2	1
QUAL	0	10.7	18	10	8	1
QUAL	0	11.1	1	0	1	0
QUAL	0	11.5	1	0	1	0
QUAL	0	11.7	716	365	351	2
QUAL	0	11.9	1	1	0	0
QUAL	0	12.2	1	1	0	0
QUAL	0	12.6	4	0	4	1
QUAL	0	13.1	1	1	0	0
QUAL	0	13.4	1	0	1	0
QUAL	0	13.5	6	1	5	0
QUAL	0	13.6	2	0	2	1
QUAL	0	14.1	2	1	1	0
QUAL	0	14.3	1	0	1	0
QUAL	0	14.4	9	2	7	0
QUAL	0	14.5	16	8	8	2
QUAL	0	14.8	1	1	0	0
QUAL	0	14.9	1	1	0	0
QUAL	0	15.0	1	0	1	0
QUAL	0	15.2	1	1	0	0
QUAL	0	15.3	2	2	0	0
QUAL	0	15.4	3	1	2	0
QUAL	0	15.5	13	1	12	2
QUAL	0	16.0	2	0	2	0
QUAL	0	16.2	7	0	7	0
QUAL	0	16.5	59	5	54	1
QUAL	0	16.9	2	1	1	0
QUAL	0	17.0	1	0	1	0
QUAL	0	17.1	1	0	1	0
QUAL	0	17.2	1	0	1	0
QUAL	0	17.3	1	0	1	0
QUAL	0	17.4	4	2	2	1
QUAL	0	18.0	1	1	0	0
QUAL	0	18.1	1	0	1	0
QUAL	0	18.2	4	3	1	0
QUAL	0	18.3	2	1	1	0
QUAL	0	18.4	14	2	12	0
QUAL	0	18.7	1	1	0	0
QUAL	0	19.2	1	1	0	0
QUAL	0	19.4	27	15	12	3
QUAL	0	19.8	2	0	2	0
QUAL	0	19.9	1	0	1	0
QUAL	0	20.0	1	0	1	0
QUAL	0	20.1	3	0	3	0
QUAL	0	20.2	2	2	0	0
QUAL	0	20.3	1	0	1	0
QUAL	0	20.4	9	3	6	1
QUAL	0	20.7	5	2	3	0
QUAL	0	20.8	8	4	4	0
QUAL	0	21.0	6	1	5	0
QUAL	0	21.1	6	1	5	0
QUAL	0	21.2	16	4	12	0
QUAL	0	21.3	2	1	1	0
QUAL	0	21.4	39	22	17	2
QUAL	0	21.7	2	0	2	0
QUAL	0	21.8	1	0	1	0
QUAL	0	22.4	3	1	2	1
QUAL	0	22.7	1	1	0	0
QUAL	0	22.9	9	4	5	0
QUAL	0	23.2	3	0	3	0
QUAL	0	23.3	2	0	2	0
QUAL	0	23.4	44	18	26	0
QUAL	0	24.2	1	1	0	0
QUAL	0	24.3	1	0	1	0
QUAL	0	24.4	26	9	17	1
QUAL	0	24.6	0	0	0	1
QUAL	0	24.8	1	0	1	0
QUAL	0	25.3	1	0	1	0
QUAL	0	25.4	5	2	3	0
QUAL	0	25.6	3	3	0	0
QUAL	0	25.8	2	0	2	0
QUAL	0	25.9	12	7	5	0
QUAL	0	26.3	1	1	0	0
QUAL	0	26.4	1	1	0	1
QUAL	0	27.4	16	5	11	1
QUAL	0	28.3	1	1	0	0
QUAL	0	28.4	8	4	4	0
QUAL	0	29.4	31	17	14	0
QUAL	0	29.9	3	3	0	0
QUAL	0	30.3	1	0	1	0
QUAL	0	30.4	414	238	176	0
QUAL	0	30.6	3	0	3	0
QUAL	0	30.7	1	1	0	0
QUAL	0	31.4	1	0	1	0
QUAL	0	31.6	1	1	0	0
QUAL	0	31.9	17	9	8	0
QUAL	0	32.4	16	5	11	4
QUAL	0	32.5	2	2	0	0
QUAL	0	32.9	1	0	1	0
QUAL	0	33.4	25	14	11	0
QUAL	0	33.9	1	1	0	0
QUAL	0	34.4	7	3	4	1
QUAL	0	34.6	1	1	0	0
QUAL	0	35.3	1	1	0	0
QUAL	0	35.4	3	3	0	0
QUAL	0	36.4	2	1	1	0
QUAL	0	37.4	49	24	25	2
QUAL	0	38.4	5	4	1	1
QUAL	0	39.4	5	4	1	3
QUAL	0	39.6	1	1	0	0
QUAL	0	40.4	1	1	0	0
QUAL	0	41.4	27	17	10	0
QUAL	0	42.2	0	0	0	1
QUAL	0	42.3	1	1	0	0
QUAL	0	42.4	5	3	2	0
QUAL	0	43.4	8	4	4	2
QUAL	0	44.2	0	0	0	1
QUAL	0	44.4	11	7	4	1
QUAL	0	45.4	43	28	15	0
QUAL	0	46.4	3	2	1	1
QUAL	0	47.2	2	1	1	0
QUAL	0	47.4	1	1	0	0
QUAL	0	48.0	0	0	0	1
QUAL	0	48.4	17	11	6	0
QUAL	0	48.5	1	1	0	0
QUAL	0	49.3	0	0	0	1
QUAL	0	49.4	3	2	1	0
QUAL	0	50.4	1	1	0	0
QUAL	0	51.3	2	1	1	0
QUAL	0	51.4	2	2	0	0
QUAL	0	51.5	1	0	1	0
QUAL	0	52.4	15	11	4	0
QUAL	0	53.3	1	1	0	0
QUAL	0	53.4	2	0	2	1
QUAL	0	54.4	2	0	2	0
QUAL	0	55.0	0	0	0	1
QUAL	0	55.4	1	1	0	0
QUAL	0	56.4	3	2	1	0
QUAL	0	57.4	5	4	1	0
QUAL	0	57.7	1	0	1	0
QUAL	0	58.4	0	0	0	1
QUAL	0	59.2	1	1	0	0
QUAL	0	59.3	0	0	0	1
QUAL	0	59.4	2	2	0	0
QUAL	0	60.0	1	0	1	0
QUAL	0	60.2	1	1	0	0
QUAL	0	60.4	11	6	5	0
QUAL	0	61.0	1	0	1	0
QUAL	0	61.3	1	1	0	0
QUAL	0	61.4	0	0	0	1
QUAL	0	62.4	5	4	1	1
QUAL	0	63.3	0	0	0	1
QUAL	0	63.4	2	1	1	0
QUAL	0	64.4	16	10	6	1
QUAL	0	65.4	0	0	0	1
QUAL	0	65.8	1	1	0	0
QUAL	0	66.4	9	6	3	0
QUAL	0	67.1	1	1	0	0
QUAL	0	67.4	12	8	4	0
QUAL	0	68.4	1	0	1	0
QUAL	0	69.4	1	0	1	1
QUAL	0	70.1	1	1	0	0
QUAL	0	70.4	1	0	1	0
QUAL	0	71.1	0	0	0	1
QUAL	0	71.4	22	15	7	0
QUAL	0	72.4	3	2	1	0
QUAL	0	73.1	1	0	1	0
QUAL	0	73.4	2	2	0	1
QUAL	0	74.4	11	4	7	1
QUAL	0	74.7	1	0	1	0
QUAL	0	75.4	9	4	5	0
QUAL	0	77.2	1	0	1	0
QUAL	0	77.4	1	1	0	0
QUAL	0	78.4	3	2	1	0
QUAL	0	79.4	2	2	0	1
QUAL	0	80.4	7	3	4	2
QUAL	0	82.4	4	4	0	0
QUAL	0	83.4	3	2	1	0
QUAL	0	84.4	1	1	0	0
QUAL	0	85.4	2	2	0	0
QUAL	0	86.4	7	4	3	0
QUAL	0	87.4	2	2	0	0
QUAL	0	88.4	1	0	1	0
QUAL	0	88.8	1	0	1	0
QUAL	0	89.4	1	1	0	0
QUAL	0	90.4	13	9	4	4
QUAL	0	93.4	0	0	0	1
QUAL	0	95.4	1	0	1	0
QUAL	0	97.4	1	0	1	0
QUAL	0	99.4	1	1	0	0
QUAL	0	100.4	2	0	2	0
QUAL	0	101.4	6	4	2	0
QUAL	0	102.8	1	1	0	0
QUAL	0	103.4	0	0	0	1
QUAL	0	104.4	1	1	0	0
QUAL	0	105.4	4	2	2	0
QUAL	0	106.4	2	2	0	1
QUAL	0	107.4	2	2	0	0
QUAL	0	109.4	2	1	1	0
QUAL	0	117.4	1	0	1	1
QUAL	0	121.4	1	1	0	0
QUAL	0	122.4	1	1	0	0
QUAL	0	123.4	3	1	2	0
QUAL	0	126.4	1	1	0	0
QUAL	0	129.4	1	1	0	0
QUAL	0	133.4	1	0	1	0
QUAL	0	137.4	1	1	0	0
QUAL	0	138.4	1	1	0	0
QUAL	0	141.4	1	1	0	0
QUAL	0	144.4	1	1	0	0
QUAL	0	151.4	1	0	1	0
QUAL	0	157.4	1	1	0	0
QUAL	0	170.4	1	1	0	0
QUAL	0	225.4	2	1	1	0
# IDD, InDel distribution:
# IDD	[2]id	[3]length (deletions negative)	[4]number of sites	[5]number of genotypes	[6]mean VAF
IDD	0	-60	1	0	.
IDD	0	-16	1	0	.
IDD	0	-12	2	0	.
IDD	0	-6	2	0	.
IDD	0	-3	2	0	.
IDD	0	-2	4	0	.
IDD	0	-1	39	0	.
IDD	0	1	21	0	.
IDD	0	2	2	0	.
IDD	0	5	1	0	.
# ST, Substitution types:
# ST	[2]id	[3]type	[4]count
ST	0	A>C	1712
ST	0	A>G	896
ST	0	A>T	984
ST	0	C>A	1063
ST	0	C>G	106
ST	0	C>T	362
ST	0	G>A	378
ST	0	G>C	117
ST	0	G>T	1092
ST	0	T>A	950
ST	0	T>C	848
ST	0	T>G	1769
# DP, depth:
#   - set id, see above
#   - the depth bin, corresponds to the depth (unless --depth was given)
#   - number of genotypes with this depth (zero unless -s/-S was given)
#   - fraction of genotypes with this depth (zero unless -s/-S was given)
#   - number of sites with this depth
#   - fraction of sites with this depth
# DP, Depth distribution
# DP	[2]id	[3]bin	[4]number of genotypes	[5]fraction of genotypes (%)	[6]number of sites	[7]fraction of sites (%)
DP	0	1	0	0.000000	8015	77.567018
DP	0	2	0	0.000000	1219	11.797155
DP	0	3	0	0.000000	534	5.167909
DP	0	4	0	0.000000	270	2.612988
DP	0	5	0	0.000000	123	1.190361
DP	0	6	0	0.000000	64	0.619375
DP	0	7	0	0.000000	27	0.261299
DP	0	8	0	0.000000	25	0.241943
DP	0	9	0	0.000000	26	0.251621
DP	0	10	0	0.000000	6	0.058066
DP	0	11	0	0.000000	7	0.067744
DP	0	12	0	0.000000	8	0.077422
DP	0	13	0	0.000000	5	0.048389
DP	0	15	0	0.000000	2	0.019355
DP	0	24	0	0.000000	1	0.009678
DP	0	26	0	0.000000	1	0.009678
```
In summary...

-The majority of variants are SNPs, and a smaller proportion are insertions/deletions.

-There is a transition-biased mutation spectrum (transition to transversion ratio of 0.32), meaning more transversions are occurring than transitions. 

-Most variants in the dataset are of good quality (with high-quality scores).

-The coverage depth across variants is low, with a large proportion of sites having just 1x depth (single reads), which may indicate low-quality sequencing data for some sites.
