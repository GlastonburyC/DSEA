# FalseTransFilter
Many *trans*-eQTLs can be 'phantom *cis*-eQTL effects' driven by multimapping reads due to gene families, paralogues and regions of high sequence similarity across the genome (*CBWD1* below). As multi-mapping is inevitable across RNA-seq experiments - These false *trans*-eQTLs are replicable across studies.

We have devised a simple post eQTL analysis filter package to mark / remove eQTLs that are potentially false, or non-interpretable.

```
python get_dosages.py dir/to/dosages/ dir/to/trans/eqtls/eQTLs.sorted output.dosages
```

```
R CMD BATCH FalseTransFilter.R gencode.v19.annotation.gtf trans.eQTLs.txt.sorted /path/to/BAMs/ output.dosages snplocations.txt listOfQuantifiedGenes.txt
```

# Functionality:

1. Overall test for psuedogene enrichment given expression background.

2. Mark trans-eQTLs that are pseudogenes

3. eQTLs that are members of gene families are also marked. If gene families are all in *cis*, multi-gene families can also create multiple spurious *cis*-eQTL signals in which no one 'causal' gene can be selected.

4. Mark phantom *cis*-eQTLs. These are *trans*-eQTLs driven by multimapping. We have observed many false *trans*-eQTL signals that are driven by shared sequence homology with regions in *cis*. Whilst the QTL effect is real (differential expression due to genotype) - It should act in *cis* rather than trans (*CBWD1* below). We automate the process of checking mean nucleotide % read coverage across all individuals in a study. If many reads align to very few basepairs, the overall coverage across the gene will be non-uniformly distributed. In addition we provide diagnostic plots of each spurious association and mark it for filtering.

# Example: *CBWD1*

An example of this effect is the following *trans*-eQTL rs2034278-CBWD1 (*P* =  1.96 x 10<sup>-19</sup>) Gtex replication ( *P* = 0.0078, Subcutaneous Adipose tissue), demonstrating these spurious signals are replicable across cohorts. When scanning for paralogous genes in cis of rs2034278, it becomes apparent that there is 1 gene, *CBWD2*. Mean basepair coverage calculated across 720 individuals clearly demonstrates that the signal is spurious, with only 3.27% of the gene having >= 1 read per basepair, even though the gene is expressed at 1CPM > 90% of individuals. It is likely that this *trans*-eQTL is actually a phantom *cis*-eQTL effect on *CBWD2*, and due to their high sequence similarity (98%, ENSEMBL) the reads are erronously assigned to a region of *CBWD1*.


![alt tag](https://raw.githubusercontent.com/GlastonburyC/FalseTransFilter/master/CBWD2.png?token=AEA_S5x_uaLXAiGYa97Xvk1hY_3UQ4Fuks5W6vIUwA%3D%3D)


This spurious signal can be compared to a *trans*-eQTL for *NME3* that did not have a paralogue within 1MB of its associated SNP and was not found to have poor/non-uniform read coverage (below).


![alt tag](https://raw.githubusercontent.com/GlastonburyC/FalseTransFilter/master/NME.png?token=AEA_S53_l77mLA8tmlymtquCJAjbvowkks5W6vM5wA%3D%3D)
