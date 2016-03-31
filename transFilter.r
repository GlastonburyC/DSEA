args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(plyr)
library(zoo)
library(ggplot2)
library(GenomicRanges)
library(gtools)
library(parallel)
library(GenomicFeatures)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

annotation.file  = "/home/glastonc/new_bams/gencode.v19.annotation.gtf" # args[1]
# FDR5 significant trans
trans_eQTL_file  = "/home/glastonc/new_bams/counts/eQTLs/trans.eQTLs.TMM.F.txt.sorted"  # args[2]
SNP_file         = "/home/glastonc/new_bams/counts/eQTLs/FDR/cisDosages.txt"
BAM.dir          = " "
quantified_genes = "/home/glastonc/new_bams/counts/Adipose.genes.txt"
ENSEMBL_paralogs = "/home/glastonc/paralog_DB.txt"
dir2BAMs         = "/home/glastonc/new_bams/bams_with_genotypes/"
dir2dosages      = "/home/glastonc/FileMatrixQTL/FilesMatrixQTL/F/dosages/"
SNPlocation_file = "/home/glastonc/FileMatrixQTL/FilesMatrixQTL/all_snps.txt"
####### load in files ########
trans_eQTL       <- fread(trans_eQTL_file,head=T,sep="\t")
trans_eQTL       <- as.data.frame(trans_eQTL)

gencode          <- fread(annotation.file,head=F,sep="\t")
gencode          <- as.data.frame(gencode)

SNPlocations     <- fread(SNPlocation_file,head=T,sep="\t")
SNPlocations     <- as.data.frame(SNPlocations)

ENSEMBL_paralogs <- fread(ENSEMBL_paralogs,head=T,sep="\t")
ENSEMBL_paralogs <- as.data.frame(ENSEMBL_paralogs)

dosages   	 <- read.table(SNP_file,head=F,sep="\t")

expressed_genes  <- read.table(quantified_genes,head=T,sep="\t")


# format gencode annotation as a GRange object - reduce to meta-exon structure.
txdb                <- makeTxDbFromGFF("/home/glastonc/new_bams/gencode.v19.annotation.gtf"
				      ,format="gtf")

exons.list.per.gene <- exonsBy(txdb,by="gene")


# replace with mclapply (parrelelise)
union_exons         <-mclapply(exons.list.per.gene,function(x){reduce(x)},mc.cores=no_cores)

cat("1) Processing Gencode annotation file..","\n")

# obtain the gene name, ensembl ID and biotype as a seperate columns in the annotation table.
gene.split <- mclapply(gencode[,9],
	      function(x)strsplit(x,";")[[1]][1],mc.cores=no_cores)

geneN      <- as.data.frame(unlist(mclapply(gene.split,
	                function(x)strsplit(x,"\"")[[1]][2],mc.cores=no_cores)))

biotype    <- as.data.frame(unlist(mclapply(gencode$V9,
	                   function(x)strsplit(strsplit(x,";")[[1]][3],
	                   "\"")[[1]][2],mc.cores=no_cores)))

gene_id    <- as.data.frame(unlist(mclapply(gencode$V9,
			   function(x)strsplit(strsplit(x,";")[[1]][5],
			   "\"")[[1]][2],mc.cores=no_cores)))


colnames(biotype)[1]="biotype"
colnames(gene_id)[1]="gene_id"
colnames(geneN)[1]="ensembl_id"

colnames(geneN) <- "gene"
gencodeM        <- cbind(biotype,gencode)
gencodeM        <- cbind(geneN,gencodeM)
gencodeM        <- cbind(gene_id,gencodeM)
gencode         <- gencodeM
rm(gencodeM)

gencode_genes <- subset(gencode,V3 == "gene")
gencode.sub   <- gencode[gencode[,2] %in% trans_eQTL[,2],]
pseudo.genes  <- unique(gencode.sub[which(gencode.sub$biotype == "pseudogene"),
	                    'biotype'])

trans_eQTL$pseudo.gene=rep('-')	

### annotate Pseudogenes that are trans-eQTL hits ###
#####################################################

trans_eQTL$pseudo.gene[which(trans_eQTL[,2] == pseudo.genes)] <- "psuedogene"

gencode_pseudo         <- subset(gencode_genes,biotype =="pseudogene")
total_pseudo           <- expressed_genes[expressed_genes[,1] %in% gencode_pseudo[,2],]
No_pseudotrans         <- length(which(trans_eQTL$pseudo.gene == "psuedogene"))
trans_length           <- length(unique(trans_eQTL[,2]))
No_of_genes_expressed  <- dim(expressed_genes)[1]

# number of trans pseudo genes, total number of psuedo genes 
# (should be background list psuedo genes) and total N trans genes)
# 1-phyper(trans.pseudo,totalPseudoGenes,1-totalPsuedoGenes,trans.N)

# Hypergeometric test for enrichment. 
# Does your trans-eQTL list have more pseudo genes than expected by chance?
enrich=1-phyper(No_pseudotrans,total_pseudo,
	     No_of_genes_expressed-total_pseudo,
	     trans_length)


cat("2) ",paste(length(which(trans_eQTL$pseudo.gene != '-')),
	    " Psuedogenes detected in trans-eQTL results..",
	    "\n","\n",sep=""))

cat("Hypergeometric enrichment of pseudogenes P = ",
	 signif(enrich,digits=3),"\n",sep="")


# extract paralogous genes.
ensemblParalogs <- function(i,gene_of_interest) {

	paralogs <- ENSEMBL_paralogs[which(ENSEMBL_paralogs[,1] == gene_of_interest)
	                           ,'Human_Paralog_Ensembl_Gene_ID']
	
	return(paralogs)
}


# This function checks whether the SNP has any paralogous genes in cis (5MB) 
# outputs all trans genes with annotation of ENSEMBL curated paralogs
# gsub("[[:digit:]]","",my.data) 
cis_window_scan <- function(gencode_genes,trans_eQTL) {

	for(i in 1:dim(trans_eQTL)[1]){

		gene_of_interest          <- substring(trans_eQTL[i,2],1,15)
		gene_id 		  <- gencode_genes$gene_id[which(substr(gencode_genes$gene,1,15)==gene_of_interest)]
		para_genes                <- tryCatch(ensemblParalogs(i,gene_of_interest),
						      error=function(e) cat('Gene has no paralogs, skipping..','\n'),
						      silent=TRUE)
		chr        		  <- SNPlocations[which(SNPlocations[,1]==trans_eQTL$SNP[i]),2]
		SNPlocations_tmp          <- subset(SNPlocations,chrm_snp==chr)
		SNP_idx                   <- which(SNPlocations_tmp[,1]==trans_eQTL[i,1])
		SNP_pos                   <- SNPlocations_tmp[,3][SNP_idx]
		gencode_genes_tmp         <- subset(gencode_genes,V1 ==chr)
		cis_genes_idx    	  <- which(gencode_genes_tmp$V4 >= SNP_pos-1e6 & gencode_genes_tmp$V4 <= SNP_pos+1e6)
		cis_genes        	  <- gencode_genes_tmp$gene[cis_genes_idx]
		para_in_cis      	  <- cis_genes[substr(cis_genes,1,15) %in% para_genes]
				if (gene_of_interest %in% substring(para_in_cis,1,15)){
			para_in_cis=para_in_cis[-which(para_in_cis==gene_of_interest)]
		}
		para_in_cis_out        	  <- paste(para_in_cis, collapse=", ")
		trans_eQTL$para_in_cis[i] <- para_in_cis_out
		cat(paste(trans_eQTL[i,2]," has ",length(para_in_cis)," paralogues",sep=""),"\n")
}
return(trans_eQTL)
}

annotated_trans_eQTLS <- cis_window_scan(gencode_genes,trans_eQTL)

# Obtain BAM file names if they contain the sample ID.
system(paste("cd ",dir2BAMs,";"," ls *.bam > bamIDs.txt;",sep=""))
bamIDs        <- read.table(paste(dir2BAMs,'bamIDs.txt',sep=''),head=F,stringsAsFactors=F)
system(paste("rm ",paste(dir2BAMs,"bamIDs.txt",sep=""),sep=""))
bamIDs        <- unlist(lapply(bamIDs$V1,function(x)strsplit(x,".bam")))

annotated_trans_eQTLS$coverage=rep(NA)

pdf("trans_eQTL.coverage.test.pdf",width=15,height=5)


read_coverage <- function(bamIDs,dosages, dir2BAMs, annotated_trans_eQTLS, gencode_genes,union_exons){

for(i in 1:nrow(annotated_trans_eQTLS)){

	# format to mpileup co-ordinates (chrX:pos1-pos2)
	cord       <- gencode_genes[which(gencode_genes[,2]==annotated_trans_eQTLS[i,2]),c(1,4,7,8)]
	cord_form  <- paste(cord[,2],":",cord[,3],"-",cord[,4],sep="")
	# extract reads from all BAMs from cord positions
	cat("Extracting reads covering gene position ",
	     cord_form," corresponding to gene: ", as.character(cord$gene_id),"\n",sep="")
	
	system(paste("module load samtools; cd ",dir2BAMs,"; samtools mpileup -r "
		   ,cord_form," *.bam > "
		   ,cord[1,1],".pileup.txt",sep=""))

	file.name     <- paste(dir2BAMs,cord[1,1],".pileup.txt",sep="")
	read.cov      <- fread(file.name,head=F,sep="\t")
	read.cov      <- as.data.frame(read.cov)
	# remove unwanted data from mpileup output
	nums          <- sapply(read.cov, is.numeric)
	read.cov      <- read.cov[,nums]
	# label columns with sample IDs
	colnames(read.cov)[2:ncol(read.cov)] <- as.character(bamIDs)
	read.cov      <- read.cov[,mixedorder(names(read.cov))]
	read.cov      <- read.cov[,c(ncol(read.cov),1:ncol(read.cov)-1)]
	x             <- paste(annotated_trans_eQTLS[i,2],sep="")
	# select gene from collapsed union exon representation
	targetGene    <- union_exons[x]
	# convert gene coordinates to GRanges object for overlap function
	bpCoverage    <- GRanges(Rle(paste(cord$V1,sep="")),IRanges(start=read.cov$V2,width=1))
	targetGene    <- GRangesList(targetGene)
	my_index      <- overlapsAny(bpCoverage, targetGene, ignore.strand=TRUE)
	my_exon       <- findOverlaps(bpCoverage, unlist(targetGene),type = "within", select = "first", ignore.strand=TRUE)
	# keep only positions that correspond to exons
	# exon.cov      <- read.cov[my_index,]
	exon.cov      <- read.cov[!is.na(my_exon),]
	exon.cov$exon <- my_exon[!is.na(my_exon)]
	exon.cov      <- exon.cov[,c(1,ncol(exon.cov),2:(ncol(exon.cov)-1))]
        exon_member   <- exon.cov[,2]

	# extract dosage for this trans-eQTL SNP
	SNP           <- dosages[which(dosages$V1==annotated_trans_eQTLS[i,1]),2:ncol(dosages)]
	SNP           <- round(SNP)
	# order samples and subset to individuals with genotypes
	colnames(SNP) <- colnames(exon.cov[,3:ncol(exon.cov)])
	# calculate the number of exon basepairs that have <= 1 read covering them.
	bpwith1       <- length(which(rowMeans(exon.cov[1:nrow(exon.cov),2:ncol(exon.cov)]) <= 1))
    gene_percent      <- (bpwith1/sum(width(targetGene)))*100
    # Width of exons for percent coverage!!
    cat(paste(round(gene_percent,digits=3),"% of nucleotides from ",cord$gene_id
    	                                  ," exons have < 1 read coverage (mean)",sep=""),"\n")

    annotated_trans_eQTLS$coverage[i]=paste(signif((1-gene_percent/100)*100,digits=3),"%",sep="")
    
    # annotate individuals with their SNP dosage and stratify by genotype.
    temp=rbind(rep(NA,ncol(exon.cov)),exon.cov)
    temp[1,2:ncol(temp)]=SNP
    colnames(exon.cov)[1]="bptmp"
    #subset read coverage by genotype calculate the mean coverage per base pair
    AA=exon.cov[,which(temp[1,]=="0")]
    AA=rowMeans(AA)
    AA=rollapply(AA, width = 1, by =1, FUN = mean, align = "left")
    AT=exon.cov[,which(temp[1,]=="1")]
    AT=rowMeans(AT)
    AT=rollapply(AT, width = 1, by =1, FUN = mean, align = "left")
    TT=exon.cov[,which(temp[1,]=="2")]
    TT=rowMeans(TT)
    TT=rollapply(TT, width = 1, by =1, FUN = mean, align = "left")

    # form a dataframe with each alleles coverage
    final<-data.frame(AA=AA,AT=AT,TT=TT,position=rollapply(exon.cov$bptmp, width = 1
    	            , by = 1, FUN = mean, align = "left"),exon=exon_member)
    colnames(final)=c('0','1','2','Position',"exon")
    dt=melt(final,id.vars=c("Position","exon"))
    colnames(dt)=c('Position','exon','dosage','value')
    recs <- data.frame(xmin=start(targetGene), xmax=end(targetGene),ymin=0,ymax=(max(dt$value)+1))

    print(ggplot() + geom_line(data=dt,aes(x=Position,y=value,colour=dosage,group=interaction(exon,dosage)),alpha=0.8,size=0.7) + 
				xlab(paste(cord$gene_id,"basepair position",sep=" ")) +
				ylab("Per nucleotide read coverage") +
				ggtitle(paste(cord$gene_id,"Read coverage",sep=" ")) + 
				theme_bw() + 
				geom_rect(data=recs, aes(xmin=xmin.value, xmax=xmax.value, ymin=ymin, ymax=ymax),
				fill='gray80', alpha=0.3))


	#system(paste("cd ",dir2dosages,"; awk '$1==\"",annotated_trans_eQTLS[i,1],"\"' all.dosages.csv > ",annotated_trans_eQTLS[i,1],".txt",sep=""))
	#SNP.file.name=paste(annotated_trans_eQTLS[i,1],".txt",sep="")
	#SNP=read.table(SNP.file.name,head=F,sep="\t")

}
dev.off()

}

read_coverage(bamIDs,dosages,dir2BAMs,annotated_trans_eQTLS,gencode_genes,union_exons)


write.table(annotated_trans_eQTLS,"/home/glastonc/new_bams/counts/eQTLs/FDR/DSEA.cis.out.txt",col.names=T,row.names=F,sep="\t",quote=F)
# mpileup function - obtain reads from BAMS for all possible trans genes
# provide diagnostic read coverage plot per dodgy gene if Perplexity score 



