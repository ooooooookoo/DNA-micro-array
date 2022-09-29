library('GEOquery')
library(Biobase)
library(limma)
library(affy)
library(oligo)
library(oligoClasses)
library(tidyverse)

gse56481 <- getGEO('GSE56481')
gse56481 <- gse56481[[1]]
methods(class=class(gse56481))

#get supplementary data and CEL files
oligoClasses::list.celfiles('GSE56481')
gse56481_cel <- read.celfiles(list.celfiles('D:/OneDrive/NUS/BL5631-R/CA1/GSE56481_RAW',full.names=TRUE,listGzipped=TRUE))
gse56481_cel
image(gse56481_cel[,1])

#metadata
pd <- pData(gse56481)
pd
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
pd['cel_file']
gse56481_cel <- read.celfiles(paste0('D:/OneDrive/NUS/BL5631-R/CA1/GSE56481_RAW/',pd$cel_file),phenoData = phenoData(gse56481))

#normalization
gse56481_cel_norm <- rma(gse56481_cel)


#specifying our model for differential gene expression analysis
library(limma)

#Using contrasts(removing the intercept)
design <- model.matrix(~ 0 + gse56481_cel_norm[['diagnosis:ch1']])
colnames(design)[1] <- 'GPA'
colnames(design)[2] <- 'health'
design

#Deciding the contrasts(when the GPA-health is positive, expression is higher in GPA than health)
contrast_matrix <- makeContrasts( GPA - health ,levels=design)
contrast_matrix

#moving ahead to fit
fit <- lmFit(gse56481_cel_norm,design)
fit2 <- contrasts.fit(fit,contrast=contrast_matrix)
fit2 <- eBayes(fit2)
fitted.ebayes <- eBayes(fit2)
topTable(fitted.ebayes)
summary(decideTests(fit2,lfc=1))
fitted.ebayes

#Annotation of up-regulation gene and down-regulating gene
BiocManager::install('hugene20sttranscriptcluster.db')
library('hugene20sttranscriptcluster.db')
ps <- rownames(topTable(fitted.ebayes))
ps
columns(hugene20sttranscriptcluster.db)
keytypes(hugene20sttranscriptcluster.db)
keys(hugene20stprobeset.db,keytype="PROBEID")
ps2 <- topTable(fitted.ebayes,number=Inf,p.value=0.05,lfc=1)
ps2_up <- rownames(ps2[ps2$logFC>0,])
df <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype = 'PROBEID')
df

ps2_down <- rownames(ps2[ps2$logFC<0,])
df2 <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype = 'PROBEID')
df2

#downstream processing(volcanoplot)
interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value=0.05,lfc=2)
# interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value=0.05,lfc=1) #for GSEA
head(interesting_genes)
volcanoplot(fitted.ebayes , coef=1 ,main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

#downstream processing(heatmap)
library(RColorBrewer)
gse_interesting <- gse56481_cel_norm[rownames(interesting_genes),]
heatmap(exprs(gse_interesting),
        labCol = gse56481_cel_norm[['diagnosis:ch1']],labRow=NA,
        col = rev(brewer.pal(10,"RdBu")),
        disfun = function(x)as.dist(1-cor(tx)))

#downstream processing(GSEA)


interesting_genes2 <- rownames_to_column(interesting_genes, var = "Gene")

library(clusterProfiler)
library(msigdbr)
library(magrittr)

  # Determine the gene sets we want to use e.g. Hallmark gene sets and C7: immunologic signature gene sets)
Gene_set_H <- msigdbr(
  species = "Homo sapiens", category = "H"
)

Gene_set_C7 <- msigdbr(
  species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB"
)

  M# Create mapped data frame from ProbeID to EntrezID
mapped_Interesting_genes2 <- data.frame(
  entrez_id = mapIds(
    hugene20sttranscriptcluster.db,
    keys = interesting_genes2$Gene,
    keytype = "PROBEID",
    column = "ENTREZID",
    multiVals = "first"
  )
) %>%
  # If doesn't match, drop the row
  dplyr::filter(!is.na(entrez_id)) %>%
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::inner_join(interesting_genes2, by = c("Ensembl" = "Gene"))

  # Check is there any duplicate in gene id 
any(duplicated(mapped_Interesting_genes2$entrez_id))

  # If TRUE, keep the only row with highest t-statistic values
filtered_mapped_Interesting_genes2 <- mapped_Interesting_genes2 %>%
  dplyr::arrange(dplyr::desc(abs(t))) %>% # Rank the t-statistic values
  dplyr::distinct(entrez_id, .keep_all = TRUE) #Filter out the first row with highest absolute value of the t-statistic


  # Create a vector ranked based on the t-statistic values
t_vector <- filtered_mapped_Interesting_genes2$t
names(t_vector) <- filtered_mapped_Interesting_genes2$entrez_id

  # Sort the t-statistic values in descending order here
t_vector <- sort(t_vector, decreasing = TRUE)


  #Running GSEA using Hallmark gene sets
set.seed(20220928)

gsea_results <- GSEA(
  geneList = t_vector,
  minGSSize = 25, maxGSSize = 500, 
  pvalueCutoff = 0.1, eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    Gene_set_H,
    gs_name,
    entrez_gene
  )
)

  #Review and store result
head(gsea_results@result)
gsea_result <- data.frame(gsea_results@result)

H_most_enrich <- gsea_result %>%
 dplyr::slice_max(n = 3, order_by = NES) # view 3 rows with the largest NES values

H_most_enrich


