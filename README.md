# FalseTransFilter
Many *trans*-eQTLs can be 'phantom *cis*-eQTL effects' driven by multimapping reads due to gene families, paralogues and regions of high sequence similarity across the genome (figure 1). As multi-mapping is inevitable across RNA-seq experiments - These false *trans*-eQTLs are replicable across studies.

We have devised a simple post eQTL analysis filter package to mark / remove eQTLs that are potentially false, or non-interpretable.

```
python get_dosages.py dir/to/dosages/ dir/to/trans/eqtls/eQTLs.sorted output.dosages
```

```
R CMD BATCH FalseTransFilter.R /path/to/annotation /path/to/transeQTLresults/ /path/to/BAMs/ path/to/dosages/
```

# Filters:

1. Overall test for psuedogene enrichment given expression background.

2. Pseudogenes originate from retrotransposition and gene duplication events. Many pseudogenes are highly similar in sequence and therefore extremely difficult to map reads to correctly/uniquely. All pseudogene-trans associations are marked to allow users to choose whether they should take them forward.

3. eQTLs that are members of gene families are also marked. If gene families are all in cis, this can also create multiple spurious cis-eQTL signals in which no one 'causal' gene can be selected.

4. Mark phantom *cis*-eQTLs. These are *trans*-eQTLs driven by multimapping. Read coverage tests are performed. We have observed many false *trans*-eQTL signals that are driven by shared sequence homology with regions in *cis*. Whilst the QTL effect is real (differential expression due to genotype) - It should act in cis rather than trans (figure X). We utilise a sliding window across the gene to assess read coverage uniformity. If many reads align to very few windows, the overall coverage across the gene will be non-uniformly distributed. In addition we provide diagnostic plots of each offender and mark it for filtering by the user.
