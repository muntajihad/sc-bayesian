# The data of the three omics are analyzed. Each by its own pipeline. 
# The results of each data were annotated with gene names.
# Then, only mutual genes that have been detected in 
# the three omics data were chosen for the downstream analysis.
# AWK, paste, SAMtools and Bedtools were used for this purpose.

##### In R #####

# load required packages
library(bnlearn)
library(tidyr)
library(dplyr)

# set working directory, and load the data.
setwd("~/omics_data/")
meth <- read.table("methylation.txt",header = TRUE)
exp <- read.table("expression.txt",header = TRUE)
cnv <- read.table("cnv.txt",header = TRUE)


### Data Preparation for BNlearn ###

## FOR DNA METHYLATION:
# categorize (as factor) methylation data
meth <- meth/100
meth <- as.data.frame(sapply(meth, function(x) cut(x, breaks = c(-1,0.20,0.80,Inf), 
                                                   labels = c(0,1,2)))) %>% mutate_if(is.character,as.factor)
# keep only genes that have shown variance through samples:
nmeth <- meth[,sapply(meth, function(x) length(unique(x)))>1]

## FOR GENE EXPRESSION:
# categorize (as factor) expression data
exp <- as.data.frame(sapply(exp, function(x) cut(x, breaks = c(-1,10,1000,Inf), labels = c(0,1,2)))) %>% mutate_if(is.character,as.factor)
# keep only genes that have shown variation through samples
nexp <- exp[,sapply(exp, function(x) length(unique(x)))>1]

## FOR CNV:
# convert cnv from integer to factor
cnv <- cnv %>% mutate_if(is.integer,as.factor)

# AFTER DATA FILTERATION KEEP ONLY GENES THAT HAVE 
# BEEN DETECTED IN THE THREE OMICS WILL BE KEPT (1901 GENES)
nnmeth <- nmeth[intersect(colnames(nmeth),colnames(nexp))]
nnexp <- nexp[intersect(colnames(nmeth),colnames(nexp))]
nncnv <- cnv[intersect(colnames(cnv),colnames(nnexp))]

########################################
# BN MODEL SCORING FUNCTION
options(scipen = 999)
scoring <- function(methy, exp, cnv,name){
  scores<-matrix(data=NA,nrow=ncol(methy),ncol=6)
  for (i in 1:ncol(methy)){
    data<-data.frame(cnv[,i],methy[,i],exp[,i])
    colnames(data)<-c("C","M","E")
    
    s1 = empty.graph(c("C","M","E"))
    s2 = empty.graph(c("C","M","E"))
    s3 = empty.graph(c("C","M","E"))
    
    # Arcs for S1 METH <-- CNV --> EXP
    arcs1 = matrix(c("C", "M", "C", "E"),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("from", "to")))
    arcs(s1) = arcs1
    
    # Arcs for S2 CNV-->METH--> EXP
    arcs2 = matrix(c("C", "M", "M", "E"),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("from", "to")))
    arcs(s2) = arcs2
    
    # Arcs for S3 CNV-->EXP--> METH
    arcs3 = matrix(c("C", "E", "E", "M"),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("from", "to")))
    arcs(s3) = arcs3
    
    
    
    fit1 = bn.fit(s1,data)
    fit2 = bn.fit(s2,data)
    fit3 = bn.fit(s3,data)
    
    sc1 = score(s1,data,type=name)
    sc2 = score(s2,data,type=name)
    sc3 = score(s3,data,type=name)
    
    rowData <- c(sc1,sc2,sc3)
    min.model=sort(rowData,decreasing = FALSE)[1]
    second.min=sort(rowData,decreasing = FALSE)[2]
    difference=second.min-min.model
    rowData <- c(rowData,min.model,second.min,difference)
    scores[i,]<-rowData
    
  }
  
  colnames(scores)<-c("INDEP","CME","CEM","Min.Model","Second.Min","Diff")
  class(scores) 
  scores=data.frame(scores)
  scores["Gene"]=colnames(cnv)
  scores$rltv.lik=exp(abs(scores$Second.Min-scores$Min.Model)/2)
  new_scores <- scores %>%
    select("Gene", everything())
  new_scores <- new_scores[order(new_scores$rltv.lik,decreasing = TRUE),]
  return(new_scores)
}

######################################################

# APPLY THE FUNCTION TO THE OMICS DATA
aic <- scoring(nnmeth,nnexp,nncnv,"aic")
bic <- scoring(nnmeth,nnexp,nncnv,"bic")

# NUMBER OF GENES FOR EACH SPECIFIC SCORE..

aic_RL10 <- aic[aic$rltv.lik>=10,] #AIC, relative likelihood >=10, 
nrow(aic_RL10) # should be 23 genes
bic_RL10 <- bic[bic$rltv.lik>=10,] #BIC, relative likelihood >=10, 
nrow(bic_RL10) # should be 27
bic_RL10 %>% inner_join(aic_RL10,by = "Gene") %>% nrow() # mutual should be 21
nrow(bic[bic$Diff>=2,]) # delta BIC >=2, should be 295 genes
nrow(bic[bic$Diff>=6,]) # delta BIC >=2, should be 13 genes
nrow(bic[bic$Diff>=10,]) # delta BIC >=2, should be 2 genes
bic_d2 <- bic[bic$Diff>=2,]

## add the name of the best model column
bic_d2$model <- names(bic_d2[,2:4])[apply(X = bic_d2[,2:4],MARGIN = 1,FUN = which.min)]

