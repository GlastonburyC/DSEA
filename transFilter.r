args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(plyr)
library(XML)

annotation.file = "/home/glastonc/new_bams/gencode.v19.proteins.RMretainedIntrons.gtf" # args[1]
# FDR5 significant trans
trans_eQTL_file = "/home/glastonc/new_bams/counts/eQTLs/trans.eQTLs.TMM.F.txt.sorted"  # args[2]
dosage.dir = " "
SNPlocation_file = "/home/glastonc/FileMatrixQTL/FilesMatrixQTL/all_snps.txt"
BAM.dir = " "
quantified_genes = " "
ENSEMBL_paralogs = "paralog_DB.txt"

####### load in files ########
trans_eQTL       <- fread(trans_eQTL_file,head=T,sep="\t")
trans_eQTL       <- as.data.frame(trans_eQTL)

gencode          <- fread(annotation.file,head=F,sep="\t")
gencode          <- as.data.frame(gencode)

SNPlocations     <- fread(SNPlocation_file,head=T,sep="\t")
SNPlocations     <- as.data.frame(SNPlocations)

ENSEMBL_paralogs <- fread(ENSEMBL_paralogs,head=T,sep="\t")
ENSEMBL_paralogs <- as.data.frame(ENSEMBL_paralogs)
cat("1) Processing Gencode annotation file..","\n")

# obtain the gene name and biotype as a seperate columns in the annotation table.
gene.split <- lapply(gencode[,9],function(x)strsplit(x,";")[[4]][1])

geneN   <- as.data.frame(unlist(lapply(gene.split,
	                     function(x)strsplit(x,"\"")[[1]][2])))

biotype <- as.data.frame(unlist(lapply(gencode$V9,
	                     function(x)strsplit(strsplit(x,";")[[1]][3],
	                     	      "\"")[[1]][2])))



colnames(biotype)[1]="biotype"

colnames(geneN) <- "gene"
gencodeM        <- cbind(biotype,gencode)
gencodeM        <- cbind(geneN,gencodeM)

gencode         <- gencodeM
rm(gencodeM)

gencode_genes <- subset(gencode,V3 == "gene")
gencode.sub   <- gencode[gencode[,1] %in% trans_eQTL[,2],]
pseudo.genes  <- unique(gencode.sub[which(gencode.sub$biotype == "pseudo_gene"),
	                    'gene'])

trans_eQTL$pseudo.gene=rep('-')	

### annotate Pseudogenes that are trans-eQTL hits ###
#####################################################

trans_eQTL$pseudo.gene[which(trans_eQTL[,2] == pseudo.genes)] = "Psuedogene"

gencode_pseudo=subset(gencode,biotype =="pseudo_gene")
total_pseudo=quantified_genes[quantified_genes[,1] %in% gencode[,1],]
No_pseudotrans         <- length(which(trans_eQTL$pseudo.gene == "Psuedogene"))
trans_length           <- length(unique(trans_eQTL[,2]))
No_of_genes_expressed  <- dim(quantified_genes)[1]

# number of trans pseudo genes, total number of psuedo genes 
# (should be background list psuedo genes) and total N trans genes)
# 1-phyper(trans.pseudo,totalPseudoGenes,1-totalPsuedoGenes,trans.N)

# Hypergeometric test for enrichment. 
# Does your trans-eQTL list have more pseudo genes than expected by chance?
enrich=1-phyper(No_pseudotrans,total_pseudo,
	     No_of_genes_expressed-total_pseudo,
	     trans_length)

cat("Hypergeometric enrichment of pseudogenes P = ",
	 enrich,"\n",sep="")

cat("2) ",paste(length(which(trans_eQTL$pseudo.gene != '-')),
	    " Psuedogenes detected in trans-eQTL results..",
	    "\n","\n",sep=""))

# extract paralogous genes.
ensemblParalogs <- function(i,gene_of_interest) {

	paralogs <- ENSEMBL_paralogs[which(ENSEMBL_paralogs[,1]==gene_of_interest),'Human_Paralog_Ensembl_Gene_ID']
	
	return(paralogs)
}


# This function checks whether the SNP has any paralogous genes in cis (5MB) 
# outputs all trans genes with annotation of ENSEMBL curated paralogs
# gsub("[[:digit:]]","",my.data) 
cis_window_scan <- function(gencode_genes,trans_eQTL) {

	for(i in 1:dim(trans_eQTL)[1]){

		gene_of_interest          <- substring(trans_eQTL[i,2],1,15)
		para_genes                <- tryCatch(ensemblParalogs(i,gene_of_interest),
									 error=function(e) cat('Gene has no paralogs, skipping..','\n'),silent=TRUE)
		SNP_idx                   <- which(SNPlocations[,1]==trans_eQTL[i,1])
		SNP_pos                   <- SNPlocations[,3][SNP_idx]
		cis_genes_idx    	      <- which(gencode_genes$V4 >= SNP_pos-5e6 & gencode_genes$V4 <= SNP_pos+5e6)
		cis_genes        	      <- gencode_genes$gene[cis_genes_idx]
		para_in_cis      	      <- cis_genes[substr(cis_genes,1,15) %in% para_genes]
		para_in_cis_out        	  <- paste(para_in_cis, collapse=", ")
		trans_eQTL$para_in_cis[i] <- para_in_cis_out

}
return(trans_eQTL)
}

annotated_trans_eQTLS <- cis_window_scan(gencode_genes,trans_eQTL)

# mpileup function - obtain reads from BAMS for all possible trans genes
# provide diagnostic read coverage plot per dodgy gene if Perplexity score 


# For given cis-eQTL results. Compute enrichment from the reference.
# Annotate which SNP is REF/ALT. Are there more + effect REF eQTLs than compared to ALT + effect eQTLs.





