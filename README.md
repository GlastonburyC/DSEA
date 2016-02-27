# FalseTransFilter
Many trans-eQTLs can be 'phantom cis-eQTL effects' driven by multimapping reads due to gene families, paralogues and regions of high sequence similarity across the genome (figure 1). As multi-mapping is common across RNA-seq experiments - These false trans-eQTLs are replicable across studies.

False trans-eQTL effects are readily detected visually by inspecting gene family occurrence around the trans-SNP. In addition, it is expected that reads should map uniformly to the trans-gene of interest, but given shared domains (e.g. Zinc Finger Domains) across large protein families, some reads map to a small portion of the expressed gene and these reads can preferentially map to a distant trans gene, whilst they actually originate from a gene much closer to the SNP (cis-effect).

Multimapping of RNA-seq reads, particularly short reads (50-75bp) is extremely common. Researchers attempt to reduce multi-mapping by primarily filtering on mapping quality (MAPQ SAM flag). However, whilst this can minimise the effects of multi-mapping it cannot remove it outright. In addition, completely removing reads that do not map uniquely can bias gene expression estimates, as multimaps make up a considerable portion of the total reads. Even in the case of uniquely mappable reads, sequencing error, PCR artifacts and natural genetic variation can cause reads to preferentially map to a highly similar sequence.

We have devised some simple post eQTL analysis filters to mark / remove eQTLs that are potentially false, or non-interpretable. As eQTLs are readily taken forward for functional genomic followup - it is essential to eliminate false positives.

# Filters:
1. Pseudogenes originate from retrotransposition and gene duplication events. Many pseudogenes are highly similar in sequence and therefore extremely difficult to map genes to. All pseudogene-trans associations are marked to allow users to choose to take them forward.

2. eQTLs that are members of gene families are also marked. If gene families are all in cis, this can also create multiple spurious cis-eQTL signals in which no one 'causal' gene can be selected.

3. Mark trans-eQTLs driven by multimapping. Read coverage tests are performed. We have observed many false trans-eQTL signals that are driven by shared sequence homology. Whilst the QTL effect is real (differential expression due to genotype) - It should act in cis rather than trans (figure X). We utilising a sliding window of read coverage to assess uniformity. If many reads align to very few windows, the overall coverage across the gene will be poisson distributed. If the read coverage for a given trans-gene approximates a poisson, we provide diagnostic plots of each offender and also mark it for filtering by the user.
