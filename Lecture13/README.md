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

