## BMMB 852: HW4 Part 2
###  09/25/24
### Kierstyn Higgins

#### Sequence Ontology
The sequence feature type I chose is **tmrna**. A tmRNA liberates a mRNA from a stalled ribosome. To
accomplish this part of the tmRNA is used as a reading frame
that ends in a translation stop signal. The broken mRNA is
replaced in the ribosome by the tmRNA and translation of the
tmRNA leads to addition of a proteolysis tag to the
incomplete protein enabling recognition by a protease.
Recently a number of permuted tmRNAs genes have been found
encoded in two parts. TmRNAs have been identified in
eubacteria and some chloroplasts but are absent from archeal
and Eukaryote nuclear genomes. The parent of this feature type is **small_regulatory_ncrna** and the child is **tmrna_region**. To investigate the feature type I did the following steps.

1. Activate *bioinfo* environment.
```bash
conda activate bioinfo
```
2. Install bio command.
```bash
pip install bio --upgrade
```
3. Download bio command.
```bash
bio --download
```
4. Choose GFF feature type and explain.
```bash
bio explain tmrna
```
This will output:
```bash
## tmrna (SO:0000584)

A tmRNA liberates a mRNA from a stalled ribosome. To
accomplish this part of the tmRNA is used as a reading frame
that ends in a translation stop signal. The broken mRNA is
replaced in the ribosome by the tmRNA and translation of the
tmRNA leads to addition of a proteolysis tag to the
incomplete protein enabling recognition by a protease.
Recently a number of permuted tmRNAs genes have been found
encoded in two parts. TmRNAs have been identified in
eubacteria and some chloroplasts but are absent from archeal
and Eukaryote nuclear genomes.

Parents:
- small_regulatory_ncrna 

Children:
- tmrna_region (part_of)
```

#### Gene Ontology
Using bio explain to investigate gene ontologies of mice.

```bash
conda activate bioinfo
```
*Choose a biological process to explain.*
```bash
bio explain GO:0009058
```
**DNA biosynthetic process:** The cellular DNA metabolic process resulting in the
formation of DNA, deoxyribonucleic acid, one of the two main
types of nucleic acid, consisting of a long unbranched
macromolecule formed from one or two strands of linked
deoxyribonucleotides, the 3'-phosphate group of each
constituent deoxyribonucleotide being joined in
3',5'-phosphodiester linkage to the 5'-hydroxyl group of the
deoxyribose moiety of the next one.

Parents:
- dna metabolic process 
- cellular macromolecule biosynthetic process 
- nucleobase-containing compound biosynthetic process 

Children:
- dna synthesis involved in dna repair 
- dna replication (has_part)
- dna amplification 
- rna-dependent dna biosynthetic process 
- dna strand elongation (has_part)
- dna polymerase activity (part_of)
- dna synthesis involved in dna replication 
- negative regulation of dna biosynthetic process (negatively_regulates)
- positive regulation of dna biosynthetic process (positively_regulates)
 
Genes: POLK

*Choose a molecular function to explain.*
```bash
bio explain GO:0003677
```
**DNA binding:** Any molecular function by which a gene product interacts
selectively and non-covalently with DNA (deoxyribonucleic
acid).

Parents:
- nucleic acid binding 

Children:
- dna secondary structure binding 
- dna template activity 
- bent dna binding 
- damaged dna binding 
- dna clamp loader activity (has_part)
- double-stranded dna binding 
- single-stranded dna binding 
- dna topoisomerase activity (has_part)
- dna binding, bending 
- recombination hotspot binding 
- dna translocase activity (has_part)
- chromatin dna binding 
- positive regulation of dna binding (positively_regulates)
- negative regulation of dna binding (negatively_regulates)
- sequence-specific dna binding 
- dna end binding 
- triplex dna binding 
- g-quadruplex dna binding 
- dna clamp unloader activity (has_part)
- flap-structured dna binding 
- histone-dependent dna binding 
- atp-dependent chromatin remodeler activity (has_part)

Genes: NEUROD2, ONECUT2, POLK

*Choose a cellular component to explain.*

```bash
bio explain GO:0016020
```

**Membrane:** A lipid bilayer along with all the proteins and protein
complexes embedded in it an attached to it.

Parents:
- cellular anatomical entity 

Children:
- prospore membrane 
- annulate lamellae 
- plasma membrane 
- outer membrane 
- extrinsic component of membrane (part_of)
- organelle membrane 
- intrinsic component of membrane (part_of)
- protein transport within lipid bilayer (occurs_in)
- phagophore assembly site membrane 
- photosynthetic membrane 
- ascus membrane 
- nuclear outer membrane-endoplasmic reticulum membrane network 
- membrane-bounded organelle (has_part)
- coated membrane 
- prospore membrane leading edge 
- respirasome (part_of)
- leaflet of membrane bilayer (part_of)
- side of membrane (part_of)
- plasma membrane region 
- membrane protein complex (part_of)
- membrane microdomain 
- pathogen-containing vacuole membrane 
- spore inner membrane

Genes: CEACAM20, INPP5D

It seems that the genome is well annotated with specific terms.
