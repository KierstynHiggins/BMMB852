## BMMB 852: Perform a differential expression analysis
### Kierstyn Higgins
#### 12/15/2024

This makefile is intended to simulate reads and generate a design and counts csv in order to perform a differential expression analysis.

1. Download the makefile.
   ```bash
   curl -o Makefile https://github.com/KierstynHiggins/BMMB852/blob/main/Lecture15/diff.mk
   ```

2. Print usage to see specific target purposes.
   ```bash
   make usage
   ```
3. Activate and prepare stats environment.
   ```bash
   make activate
   make tools
   make stat
   ```

4. Simulate reads and create edger file.
   ```bash
   make simulate
   make edger
   ```

5. Generate plots.
   ```bash
   make pca
   make heat
   ```

   This will create both a pca plot and a heat map as shown below.
   <img width="764" alt="Screenshot 2024-12-15 at 9 16 05 PM" src="https://github.com/user-attachments/assets/989e3c68-c1b0-4eaa-89c2-8afe6fa041b8" />
   <img width="603" alt="Screenshot 2024-12-15 at 9 15 57 PM" src="https://github.com/user-attachments/assets/e7f555fa-c931-401a-9cde-2c9b8750da6c" />

We found 61 genes. It appears that approximately half of the genes are expressed in group A and the other half in group B (if I'm interpreting it correctly?). The PCA shows well-defined clusters which suggests reliable data.
