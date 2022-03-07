####################################################
##
##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
##   Version 1.0 08.07.2016
##   Needs packages ggplot2,grid,gridExtra
##
####################################################

library(ggplot2)
library(grid)
library(gridExtra)

## calculate Immunophenoscore
#ipsmap<- function (x) {
#	if (x<=0) {
#		ips<-0
#	} else {
#		if (x>=3) {
#		 ips<-10
#		} else {
#			ips<-round(x*10/3, digits=0)
#		}
#	}
#	return(ips)
#}


## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
IPSCore <- function(data, file = '../0.data/IPS_genes.txt') {
	gene_expression<-data %>% as.data.frame
	sample_names<-names(gene_expression)

	## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
	# For different 
	IPSG<-read.table(file,header=TRUE, sep="\t", dec = ".",check.names=FALSE)
	unique_ips_genes<-as.vector(unique(IPSG$NAME))

	IPS<-NULL
	MHC<-NULL
	CP<-NULL
	EC<-NULL
	SC<-NULL
	AZ<-NULL

	# Gene names in expression file
	GVEC<-row.names(gene_expression)
	# Genes names in IPS genes file
	VEC<-as.vector(IPSG$GENE)
	# Match IPS genes with genes in expression file
	ind<-which(is.na(match(VEC,GVEC)))
	# List genes missing or differently named
	MISSING_GENES<-VEC[ind]
	dat<-IPSG[ind,]
	if (length(MISSING_GENES)>0) {
		cat("differently named or missing genes: ",MISSING_GENES,"\n")
	}

	for (i in 1:length(sample_names)) {	
		GE<-gene_expression[[i]]
		mGE<-mean(GE)
		sGE<-sd(GE)
		Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
		W1<-IPSG$WEIGHT
		WEIGHT<-NULL
		MIG<-NULL
		k<-1
		for (gen in unique_ips_genes) {
			MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
			WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
			k<-k+1
		}
		WG<-MIG*WEIGHT
		MHC[i]<-mean(WG[1:10])
		CP[i]<-mean(WG[11:20])
		EC[i]<-mean(WG[21:24])
		SC[i]<-mean(WG[25:26])
		AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
		#IPS[i]<-ipsmap(AZ[i])
		IPS[i]<- AZ[i]
	}
	DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ)
	return(DF)
}
