# FalseTransFilter
Many *trans*-eQTLs can be 'phantom *cis*-eQTL effects' driven by multimapping reads due to gene families, paralogues and regions of high sequence similarity across the genome (figure 1). As multi-mapping is inevitable across RNA-seq experiments - These false *trans*-eQTLs are replicable across studies.

We have devised a simple post eQTL analysis filter package to mark / remove eQTLs that are potentially false, or non-interpretable.
```
R CMD BATCH FalseTransFilter.R /path/to/annotation /path/to/transeQTLresults/ /path/to/BAMs/ path/to/dosages/
```

# Filters:
1. Pseudogenes originate from retrotransposition and gene duplication events. Many pseudogenes are highly similar in sequence and therefore extremely difficult to map genes to. All pseudogene-trans associations are marked to allow users to choose to take them forward.

2. eQTLs that are members of gene families are also marked. If gene families are all in cis, this can also create multiple spurious cis-eQTL signals in which no one 'causal' gene can be selected.

3. Mark *trans*-eQTLs driven by multimapping. Read coverage tests are performed. We have observed many false *trans*-eQTL signals that are driven by shared sequence homology. Whilst the QTL effect is real (differential expression due to genotype) - It should act in cis rather than trans (figure X). We utilising a sliding window of read coverage to assess uniformity. If many reads align to very few windows, the overall coverage across the gene will be poisson distributed. If the read coverage for a given trans-gene approximates a poisson, we provide diagnostic plots of each offender and also mark it for filtering by the user.
