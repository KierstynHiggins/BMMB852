# BMMB 852 HW3
## Kierstyn Higgins

Download genome for *Staphylococcus aureus*.
```bash
datasets download genome accession GCF_000013425.1
```
Unzip the file.
```bash
unzip ncbi_dataset.zip
```
Include other feature types.
```bash
datasets download genome accession GCF_000013425.1 --include gff3,cds,protein,rna,genome
```
Unzip again...
```bash
 unzip ncbi_dataset.zip
 ```
Build a FASTA index for the genome.
```bash
samtools samtools faidx ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna
```
Extract features of type "gene".
```bash
 cat ncbi_dataset/data/GCF_000013425.1/genomic.gff | awk '$3 == "gene"
' > genes.gff
 ```
 Extract features of type "CDS".
 ```bash
  cat ncbi_dataset/data/GCF_000013425.1/genomic.gff | awk '$3 == "CDS"
' > cds.gff
```
 Link genome file under new name.
 ```bash
  ln -s ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna genome.fa
  ```
  Build index for the genome.
  ```bash
  samtools faidx genome.fa
  ```
My next step was to visualize the genome in IGV. You can see that the coding sequences line up with the translation table and start with a start codon and end with a stop codon.
![Screenshot 2024-09-15 210801](https://github.com/user-attachments/assets/b9cfbf9b-806f-48c0-81ee-4d9e395aab34)

Next I created my own demo gff file using:
```bash
code demo.gff
```
I formatted the gff to look like this, using tabs instead of spaces.
![Screenshot 2024-09-15 213753](https://github.com/user-attachments/assets/ac35f17b-96ea-4cfb-8ebf-3a12c18ec17d)

Then I visualized my gff file in IGV.
![Screenshot 2024-09-15 213544](https://github.com/user-attachments/assets/5112a8fc-0dd3-4d3d-8bb4-dd8506c65db5)

