setwd("G:/Semester 2 Courses/Gene Expression/project/group4_T/group4_T/")


setRepositories(ind = c(1,2,3,4))

pkgs<- c("GEOquery", "biomaRt", "DESeq2", "BiocParallel",
         "openxlsx", "topGO", "goseq", "GO.db", "org.Hs.eg.db")

options(timeout = 1200)
install.packages(pkgs)

lapply(pkgs, require, character.only = TRUE)


##### Connecting to the GEO series and download data
library(GEOquery)
gse <- getGEO("GSE95640",GSEMatrix=FALSE)
#### Check contents of object
# check meta data for the object
head(Meta(gse))
# Check GSM lists (Data for individual sample) in gse
head(names(GSMList(gse)))
# See contents for 1st sample (GSM)
GSMList(gse)[[1]]
# Check meta data for a sample
names(Meta(GSMList(gse)[[1]]))
# Get sample characteristics within data
Meta(GSMList(gse)[[1]])$characteristics_ch1

# See the included platform (GPL)
names(GPLList(gse))
#### Building our analysis data
#Build sample annotation matrix
samAnno<- t(sapply(GSMList(gse), function(tmp)
  return(c(Meta(tmp)$characteristics_ch1, Meta(tmp)$title))))

# See first few rows
head(samAnno)
# See frequency of each annotation except for age
table(samAnno[,-c(3,5)])
# Make simplified sample annotation data frame relevant to
# the analysis
samAnnoDf<- data.frame(
  title=as.character(samAnno[,5]),
  gender=gsub("gender: ", "", samAnno[,2]),
  age=as.numeric(gsub("age: ", "", samAnno[,3])),
  time= gsub("without ", "NO", gsub("time: ", "",
                                    gsub("( is for the baseline | is after 8 weeks of )", "",
                                         samAnno[,4])))
)
# Check first few rows of annotation table
head(samAnnoDf)
# Get read-count data
download.file(url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95640/suppl/GSE95640%5Fcid12%5Fcng%5Frawcnts%5Fqced%2Etsv%2Egz",
              destfile="G:/Semester 2 Courses/Gene Expression/project/GSE95640racnts.txt.gz")

# Load read-count data
cnt<- read.table(gzfile("GSE95640racnts.txt.gz"), header=T,
                 sep=",", stringsAsFactors=F)
# See dimensions of count data
dim(cnt)
# See the first few rows of count data
head(colnames(cnt))
cntRowNames<- cnt[,1]
cnt<- cnt[,-1]
rownames(cnt)<- cntRowNames
#Reorder the annotation table to match cnt columns
rownames(samAnnoDf)<- as.character(samAnnoDf[,1])
samAnnoDf<- samAnnoDf[colnames(cnt),]

# Check if order of columns in count data frame is
# similar to that of sample annotation data frame
# This should return TRUE
all(colnames(cnt)== as.character(samAnnoDf[,1]))
#### Select your samples to work with
# Check the test and ctrl frequency
(tabTime<- table(samAnnoDf[,"time"]))
# Load the analysis samples
sampleIds<- scan(file="sampleIds.txt", what="character")
# select the analysis sub-samples
cntDat<- cnt[,sampleIds]
annoDat<- samAnnoDf[sampleIds,]




##question1:
length(which(rowSums(cntDat)==0))

##question2:
13159/nrow(cntDat)*100



cntDat<- cntDat[-which(rowSums(cntDat)==0),]


library(biomaRt)
ensembl<-useMart("ENSEMBL_MART_ENSEMBL",
                 host="http://grch37.ensembl.org", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( attributes=c("ensembl_gene_id",
                                       "external_gene_name", "gene_biotype",
                                       "chromosome_name", "start_position", "end_position"),
                          filters = "ensembl_gene_id",
                          values = unique(as.character(rownames(cntDat))),
                          mart=ensembl, uniqueRows = TRUE)
# Are there the same number of genes in the count and
# gene-annotation data ?
nrow(annoGen)== nrow(cntDat)
# Are the orders of the genes similar in the count and
# gene-annotation data ?
table(rownames(cntDat)==annoGen$"ensembl_gene_id")



##question3
length(which(annoGen$gene_biotype=="rRNA"))


##question4
length(which(annoGen$gene_biotype=="Mt_rRNA"))


##question5
cntDat_rRNA<- cntDat[which(annoGen$gene_biotype=="rRNA"),]
sum(cntDat_rRNA)


##question6
cntDat_MtrRNA<- cntDat[which(annoGen$gene_biotype=="Mt_rRNA"),]
sum(cntDat_MtrRNA)


rRNA<- c(which(annoGen$gene_biotype=="rRNA"),which(annoGen$gene_biotype=="Mt_rRNA"))


fil<- cntDat[-rRNA,]

save.image(file = "gene expreession project.RData")

load("G:/Semester 2 Courses/Gene Expression/project/group4_T/group4_T/gene expreession project.RData")

rownames(annoGen)<- annoGen$ensembl_gene_id
annoGenFil<- annoGen[rownames(fil),]



###Data normalization
library(DESeq2)
vstFil<- vst(as.matrix(fil))
plot(x=rowMeans(vstFil),y=apply(vstFil,1,sd),xlab="mean of the expression of the genes",
     ylab="standard deviation of the expression of the genes",main= "VST transformed data",
     col="blue")

plot(x=rowMeans(fil),y=apply(fil,1,sd),xlab="mean of the expression of the genes",
     ylab="standard deviation of the expression of the genes",main= "filtered untransformed
     data",col="blue")


pdf("VSTnormalization_sd_vs_mean.pdf")
par(mfrow=c(2,1))

plot(x=rowMeans(fil),y=apply(fil,1,sd),xlab="mean (filtered untransformed data)",
     ylab="sd",col="blue")

plot(x=rowMeans(vstFil),y=apply(vstFil,1,sd),xlab="mean (VST transformed data)",
     ylab="sd",col="blue")

dev.off()


# Run PCA on transpose of the
genVstFilPca<- prcomp(t(vstFil), scale=TRUE)
#Measure the variance and plot the variance of the data that
# is explained by each PC dimension
pcVar<- genVstFilPca$sdev^2
pdf("PC_Var.pdf")
barplot(pcVar, names=colnames(genVstFilPca$x), ylab="Var")
dev.off()
# Measure the variance percentage and
# plot the percentage of the variance within the data that
# is explained by each PC dimension
varPerc<- 100*pcVar/sum(pcVar)
pdf("PC_Var_Perc.pdf")
barplot(varPerc, names=colnames(genVstFilPca$x), ylab="Var %")
dev.off()



pdf("pairs_time.pdf")
pairs(genVstFilPca$x[,1:4], col=as.numeric(as.factor(annoDat$time)),oma=c(4,4,6,12))
par(xpd=T)
legend(0.805, 0.2, as.vector(c("CID1NOLCD","CID2LCD")),
       fill=c("black","red"))
dev.off()



pdf("pairs_gender.pdf")
pairs(genVstFilPca$x[,1:4], col=as.numeric(as.factor(annoDat$gender)),oma=c(4,4,6,12))
par(xpd=T)
legend(0.82, 0.2, as.vector(c("F","M")),
       fill=c("black","red"))
dev.off()


dds<- DESeq2::DESeqDataSetFromMatrix(countData = fil, colData =
                                       annoDat[,c("gender","time")], design = ~gender+time)
library(BiocParallel)
dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=2))
ddsDiff<- DESeq2::results(dds, name="time_CID2LCD_vs_CID1NOLCD", pAdjustMethod =
                            "BH", parallel = TRUE, BPPARAM = SnowParam(workers=2))
# Adjust and correct Log Fold change relative to mean expression
ddsAshr <- DESeq2::lfcShrink(dds, coef="time_CID2LCD_vs_CID1NOLCD", type="ashr",
                             parallel = TRUE, BPPARAM = SnowParam(workers=2))

save.image(file = "gene expreession project1.RData")

load("gene expreession project1.RData")


colnames(ddsAshr)

##question 9
padj<- ddsAshr$padj
padj<- padj[complete.cases(padj)]
padj_0.05<- padj[padj<0.05]
length(padj_0.05)

##question 10
pvalue<- ddsAshr$pvalue
pvalue<- pvalue[complete.cases(pvalue)]
pvalue_0.05<- pvalue[pvalue<0.05]
length(pvalue_0.05)

##question 11
diff<- data.frame(ddsAshr)
diff_upregulated<- subset(diff,diff$log2FoldChange>0 & diff$pvalue<0.05)
##or
diff_upregulated<- diff[which(ddsAshr$pvalue <0.05 & ddsAshr$log2FoldChange > 0),]

nrow(diff_upregulated)

##question 12

diff<- data.frame(ddsAshr)
diff_downregulated<- subset(diff,diff$log2FoldChange<0 & diff$pvalue<0.05)
##or
diff_downregulated<- diff[which(ddsAshr$pvalue <0.05 & ddsAshr$log2FoldChange < 0),]

nrow(diff_downregulated)

pdf("LogFC_vs_mean_expr_Unadjusted.pdf")
DESeq2::plotMA(ddsDiff)
dev.off()

pdf("LogFC_vs_mean_expr_ASHR_adjusted.pdf")
DESeq2::plotMA(ddsAshr)
dev.off()


# Bind two data frames by columns
outDf<- cbind(annoGenFil, as.data.frame(ddsAshr))
# Write data frame into a text file
write.table(outDf[which(ddsAshr$pvalue<0.05),],
            file="sigGenes.tsv", col.names=T, row.names=F, quote=F, sep="\t")


outDf<- cbind(annoGenFil, as.data.frame(ddsAshr))
write.table(outDf[which(ddsAshr$pvalue<0.05),],
            file="sigGenes.tsv",
            col.names=T, row.names=F, quote=F, sep="\t")
library(goseq)
bpAnnoList<-
  getgo(genes=unique(as.character(outDf[,"ensembl_gene_id"])),
        genome="hg19",
        id='ensGene', fetch.cats="GO:BP")
upGenes<- as.numeric(outDf$log2FoldChange>0)
names(upGenes)<- outDf[,1]
dnGenes<- as.numeric(outDf$log2FoldChange<0)
names(dnGenes)<- outDf[,1]
### First GO enrichent for Up-regulated genes
# Assigning weights based on the length of genes
pwfUp<- nullp(upGenes,"hg19","ensGene", plot.fit=FALSE)
# Using the default Wallenius approximation
goWallUp<-goseq(pwfUp,"hg19","ensGene")
# Adjust for multiple testing
padj<-p.adjust(goWallUp$over_represented_pvalue,
               method="BH")
goWallUp<- cbind(goWallUp, padj)
enrichedGoUp<- goWallUp[which(padj<.01),]


##question13
nrow(enrichedGoUp)

# See the most significant GO categories
head(enrichedGoUp)
enrichedGoCatUp<- enrichedGoUp$category
library(GO.db)
# Show descriptions and reports of 10 top GOs enriched in
# significantly upregulated genes
for(go in
    enrichedGoCatUp[1:min(10,length(enrichedGoCatUp))]){print(GOTERM[[
      go]])
  cat("--------------------------------------\n")}
save.image("gene expreession project2.RData")




pwfDn<- nullp(dnGenes,"hg19","ensGene", plot.fit=FALSE)
# Using the default Wallenius approximation
goWallDn<-goseq(pwfDn,"hg19","ensGene")
# Adjust for multiple testing
padj<-p.adjust(goWallDn$over_represented_pvalue,
               method="BH")
goWallDn<- cbind(goWallDn, padj)
enrichedGoDn<- goWallDn[which(padj<.01),]


##question14
nrow(enrichedGoDn)

# See the most significant GO categories
head(enrichedGoDn)
enrichedGoCatDn<- enrichedGoDn$category
library(GO.db)
# Show descriptions and reports of 10 top GOs enriched in
# significantly downregulated genes
for(go in
    enrichedGoCatDn[1:min(10,length(enrichedGoCatDn))]){print(GOTERM[[
      go]])
  cat("--------------------------------------\n")}


library(openxlsx)
openxlsx::write.xlsx(x=list(down_regulated_significant_GO=enrichedGoDn,up_regulated_significant_GO=enrichedGoUp),
                     file="./EnrichedGO.xlsx", colNames = TRUE, keepNA=TRUE)

save.image("gene expreession project3.RData")



for(go in
    c("GO:0006119","GO:0009055","GO:0042773","GO:0022904")){print(GOTERM[[
      go]])
  cat("--------------------------------------\n")}
