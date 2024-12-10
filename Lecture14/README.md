## BMMB 852: RNA Seq Assignment
### Kierstyn Higgins
### 12/08/2024

Using the makefile "rnaseq.mk", you will be able to perform an RNA-Seq analysis on each dataset and generate a ccount matrix that
summarizes read counts for each dataset.

1. Download the makefile.
   ```bash
   curl -L -o Makefile 
   https://github.com/KierstynHiggins/BMMB852/blob/main/Lecture14/rnaseq.mk
   ```
2. Update the variables for the dataset of your choice.
   
3. Print usage for each target.
   ```bash
   make usage
   ```
4. Run the makefile.
   ```bash
   make all
   ```
5. Run the makefile in parallel with the design file.
   ```
   cat design.csv | head -3 | \
    parallel --lb -j 1 --colsep , --header : \
    make all SRR={run_accession} SAMPLE={sample_alias}
   ```
   This should generate a file called res/counts-hisat.csv. However when I do this, it still only lists one sample in the csv.
   ![image](https://github.com/user-attachments/assets/345e792a-c6fc-4f49-907f-d02e4e315eb4)

6. Lastly you can visualize in IGV. However, only one of my bam files is showing up so I'm not sure what to do here...
   ![image](https://github.com/user-attachments/assets/1b5bad70-558e-4824-a4e4-fe790716f997)
