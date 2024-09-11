# **BMMB 852 HW2**
### Kierstyn Higgins
1. *Cyprinodon variegatus*, aka the sheephead minnow or sheephead pupfish, is a species of ray-finned fish. It can be found in the salt marshes on the east coasts of North and Central America.
   It can grow up to 3 inches and is best known for its ability to adapt to varying salinity levels.
2. 821648 features 
3. I am not sure why when I got to this step it gave me an error (see line 98).
4. There are 23019 genes listed.
5. The top ten most common feature types are exon, CDS, mRNA, five_prime_UTR, gene, three_prime_UTR, biological_region, region, ncRNA_gene, and snoRNA.
6. Yes this seems to be a well-annotated organism.

## Code 

Download the file from emsebl.org.

```bash
wget https://ftp.ensembl.org/pub/current_gff3/cyprinodon_variegatus/Cyprinodon_variegatus.C_variegatus-1.0.112.gff3.gz
 ```
Activate the environment.
```bash
conda activate bioinfo
```
List contents of the directory.
```bash
ls
```
This will print...
`
Cyprinodon_variegatus.C_variegatus-1.0.112.gff3.gz
bin
bioinformatics
bioniformatics
edirect
edu
micromamba
parent_dir
studies
work
(bioinfo)
`

Unzip the file.
```bash
gunzip Cyprinodon_variegatus.C_variegatus-1.0.112.gff3.gz
```
Let's see how big the file is..
```bash
ls -l
```
This prints...
`total 124M
-rw-r--r-- 1 ecouser ecouser 124M Mar  5  2024 Cyprinodon_variegatus.C_variegatus-1.0.112.gff3
drwxr-xr-x 3 ecouser ecouser 4.0K Aug 27 16:05 bin
drwxr-xr-x 4 ecouser ecouser 4.0K Aug 28 17:39 bioinformatics
drwxr-xr-x 3 ecouser ecouser 4.0K Aug 28 18:18 bioniformatics
drwxr-xr-x 7 ecouser ecouser 4.0K Aug 27 16:05 edirect
drwxr-xr-x 3 ecouser ecouser 4.0K Aug 28 11:57 edu
drwxr-xr-x 6 ecouser ecouser 4.0K Aug 27 16:04 micromamba
drwxr-xr-x 3 ecouser ecouser 4.0K Aug 28 16:56 parent_dir
drwxr-xr-x 3 ecouser ecouser 4.0K Aug 28 17:03 studies
drwxr-xr-x 7 ecouser ecouser 4.0K Aug 27 16:15 work
(bioinfo)`

Let's see how many lines of code this file is...
```bash
cat Cyprinodon_variegatus.C_variegatus-1.0.112.gff3 | wc -l
```
This tells us...
`821648
(bioinfo)`

We can print the head...
```bash
cat Cyprinodon_variegatus.C_variegatus-1.0.112.gff3 | head
```

Or we can print the tail!
```bash
cat Cyprinodon_variegatus.C_variegatus-1.0.112.gff3  | tail
```
We will see that the sequence regions start with `JPKM` and a series of numbers.

`##sequence-region   JPKM01095229.1 1 78823
##sequence-region   JPKM01097702.1 1 55185
##sequence-region   JPKM01097741.1 1 55057
##sequence-region   JPKM01100000.1 1 33715
##sequence-region   JPKM01101084.1 1 24462
##sequence-region   JPKM01101207.1 1 23389
##sequence-region   JPKM01101290.1 1 22677
##sequence-region   JPKM01101291.1 1 22661
##sequence-region   JPKM01102216.1 1 15975`

How many sequence regions does this file have?
```bash
cat Cyprinodon_variegatus.C_variegatus-1.0.1
12.gff3 | grep JPKM | wc -l
```
However, here is where I got an error. 
`cat: Cyprinodon_variegatus.C_variegatus-1.0.1: No such file or directory
12.gff3: command not found
0
(bioinfo)`

Here we can remove sequence regions and shorten the gff3 file name.
 ```bash
 cat Cyprinodon_variegatus.C_variegatus-1.0.112.gff3 | grep
-v '#' > cypr.gff3

cat cypr.gff3 | cut -f 1,2 | head
```
Here we can print the feature types.
```bash
cat cypr.gff3 | cut -f 3 | head
```
Which will print
`
region
gene
mRNA
exon
CDS
exon
CDS
gene
mRNA
exon
`
This will sort the feature types and put them in order of largest to smallest.
```bash
cat cypr.gff3 | cut -f 3 | sort | uniq -c | sort -rn | head
```
Which should look like this:
`
 333767 exon
 326230 CDS
  32898 mRNA
  24017 five_prime_UTR
  23019 gene
  17908 three_prime_UTR
  11664 biological_region
   9259 region
    404 ncRNA_gene
    175 snoRNA
`

