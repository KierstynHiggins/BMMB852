# BMMB 852: Fastqc demo
## Kierstyn Higgins

Set accession and directories. I chose the genome for **Brachycephalus brunneu**, also known as the pumpkin toadlet, native to the Brazilian Atlantic Forest.
```bash 
ACCESSION=SRR5837600
FASTQ1=SRR5837600_1.fastq
FASTQ2=SRR5837600_2.fastq
FASTQ1_TRIM=SRR5837600_1.trimmed.fastq
FASTQ2_TRIM=SRR5837600_2.trimmed.fastq
```
Activate environment.
```bash
conda activate bioinfo
```
Use fast-dump to download dataset.
```bash
fastq-dump -X 10000 --split-files ${ACCESSION}
```
Check quality of fastq file.
```bash
fastqc ${FASTQ1}
```
This file appears to have per base sequence content, overrepresented sequences, and adaptor content.
![Screenshot 6-10-2024](https://github.com/KierstynHiggins/BMMB852/blob/main/Lecture7/Screenshot_6-10-2024_22236_.jpeg?raw=true)

Now, Use fastp to trim files.
```bash
fastp --cut_tail -i ${FASTQ1} -o ${FASTQ1_TRIM} -I ${FASTQ2} -O ${FASTQ2_TRIM}
```
Check quality of fastq file post-trim.
```bash
fastqc ${FASTQ1_TRIM}
```
Now the file appears much cleaner.
![Screenshot](https://github.com/KierstynHiggins/BMMB852/blob/main/Lecture7/Screenshot_6-10-2024_22236_.jpeg?raw=true)
