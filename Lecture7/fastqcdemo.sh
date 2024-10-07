# BMMB:852 fastqc demo
# Kierstyn Higgins

# Set accession and directories.
ACCESSION=SRR5837600
FASTQ1=SRR5837600_1.fastq
FASTQ2=SRR5837600_2.fastq
FASTQ1_TRIM=SRR5837600_1.trimmed.fastq
FASTQ2_TRIM=SRR5837600_2.trimmed.fastq
# Activate environment.
conda activate bioinfo

# Use fast-dump to download dataset.
fastq-dump -X 10000 --split-files ${ACCESSION}

# Check quality of fastq file.
fastqc ${FASTQ1}

# Use fastp to trim files.
fastp --cut_tail -i ${FASTQ1} -o ${FASTQ1_TRIM} -I ${FASTQ2} -O ${FASTQ2_TRIM}

# Check quality of fastq file post-trim.
fastqc ${FASTQ1_TRIM}
