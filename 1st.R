# Fiest part : create correspondence file for each dataset

# 1. Set working directory
setwd("/home/dm20/Desktop/Recent_work/Bioinformatics/Scleroderma/Diseases_folder")

# 2. Load required libraries 
library(RCurl)
library(GEOquery)
library(limma)
library(topGO)
library(genefilter)

# Perform 3.1 for micro array data
# perform 3.2 for RNAseq data

# 3.1 For micro-array data perform below steps
gse <- getGEO(filename = "GSE104174_series_matrix.txt.gz",destdir=getwd())
d <- factor(c(rep('CTRL', 7),rep('SSc',10)))


mod <- model.matrix(~0+d)
# print(mod)
# dSSc dCTRL dCVCTRL  
# 1    0     0       0     1
# 2    0     0       1     0

fit_1 <- lmFit(gse, mod)
# #print(fit_1)
contr <- makeContrasts(dCTRL-dT2D,levels = mod)
#contr <- makeContrasts(dCD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)
# 
# #> colnames(fit_3)
# #[1] "dSSc - dCTRL"   "dTREAT - dSSc"  "dTREAT - dCTRL"
# #print(fit_3)
# #show(fit_3)
# #n=length(fit_3)
# #print(n)
# #topTable(fit_3, coef=1, adjust="BH") #coef=1 means dSSc-dCTRL
# #topTable(fit_3, coef=2, adjust="BH") #coef=1 means dTRAET-dSSc
# 
#table_result <- topTable(fit_3, coef=1,p.value=1e-2, sort.by = "logFC")
#table_result <- topTable(fit_3, coef=1,p.value=1e-2, sort.by = "logFC")
#table_result <- topTable(fit_3, coef=1, n=Inf,p.value=5e-2, sort.by = "logFC")
#dim(table_result)
#[1] 156  32
table_result <- topTable(fit_3, coef=1,n=Inf,adjust="BH", sort.by = "logFC")
#> dim(table_result)
#[1] 18981    32
#[1] 156  32
#dim(table_result)
# # #[1]  9 11


subtable_result <- subset(table_result, select=c("ID","Gene.Symbol","adj.P.Val","P.Value","t","B","logFC"))
subtable_result <- subtable_result[is.na(subtable_result$ID)==0,]
# End of micro-array data part

# 3.2 For RNAseq : Start from here 
# data from Biojupies / Grein
subtable_result <- read.table("Control vs Perturbation.txt", 
                 header = TRUE)

geneList <- subtable_result$logFC
#names(geneList) <- subtable_result$Gene.Symbol
names(geneList) <- subtable_result$gene_symbol

# 4. Save subtable
write.csv(subtable_result,"GSE104174_table.csv")
# nrow(subtable_result[subtable_result$P.Value<0.05,])

# 5. Apply condition of LogFC threshold 
topDiffGenes <- function(allScore){
  return(abs(allScore)>1.00)
}

# 6. Create topGO class with annotation
CD_GOdata <- new("topGOdata",
                 description = "cd study",
                 ontology = "BP", 
                 allGenes = geneList,
                 geneSel = topDiffGenes,
                 annot = annFUN.org, 
                 ID = "symbol", 
                 mapping = "org.Hs.eg.db",
                 nodeSize = 10)

# 7. Create genes and GO terms
n_sg <- sum(topDiffGenes(geneList))
sg <- sigGenes(CD_GOdata)
ug <- usedGO(CD_GOdata)

# 8. Fisher's Exact Test (and Kolmogorov-Smirnov for evaluation, not described)
resultFisher <- runTest(CD_GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(CD_GOdata, algorithm = "classic", statistic = "ks")


# 9. GO terms tree
allRes <- GenTable(CD_GOdata, classic = resultFisher, KS = resultKS, orderBy = "classic", topNodes = 30)
showSigOfNodes(CD_GOdata,score(resultFisher),firstSigNodes = 5, useInfo = "all")
printGraph(CD_GOdata,resultFisher,firstSigNodes = 5, fn.prefix = "CD1_GSE11501", useInfo = "all", pdfSW = TRUE)

# 10. Create text file for the correspondence GO terms - genes (this file is mandatory for the 2nd script)
terms <- allRes$GO.ID
genes <- genesInTerm(CD_GOdata,terms)
for (i in 1:length(terms))
{
  term <- terms[i]
  genes_term <- genes[term][[1]]
  # find the genes that are in the list of genes of interest
  fact <- genes_term %in% sg
  genes_term_2 <- genes_term[fact == TRUE]
  genes_term_2 <- paste(genes_term_2, collapse=',')
  cat(term,"genes:",genes_term_2,"\n", append = TRUE, file = "CD1_GSE11501_correspondence.txt" )
}
