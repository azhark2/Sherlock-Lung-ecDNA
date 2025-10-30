library(BiocManager)
library(tidyverse)
library(scales)
library(ggrepel)
library(broom)
library(data.table)
library(AnnotationHub)
library(org.Hs.eg.db)
library(edgeR)
library(dplyr)
library(limma)
library(ggbiplot)

matrix <- read.table("~/iCloud/dev/Sherlock-Lung/results2/Sherlock_EAGLE_RNAseq_tumor_only_matrix_filtered.txt", header = TRUE, row.names = 1, sep = "\t")

raw_counts <- read.table("~/iCloud/dev/Sherlock-Lung/data2/Sherlock1+2+TCGA-LUAD+EAGLE_RAW_RNA_COUNTS.txt", header = TRUE, sep = "\t")

data <- read.csv("~/iCloud/dev/Sherlock-Lung/metadata/rna2wgs.txt", sep = '\t', header=TRUE)

ecDNA <- read.delim("~/iCloud/dev/Sherlock-Lung/results2/ecDNA-annotations.txt")

metadata <- read.csv("~/iCloud/dev/Sherlock-Lung/results2/log-reg-table5.tsv", sep = '\t', header=TRUE)
metadata = metadata %>% dplyr::select(Tumor_Barcode, Gender, Age, Histology, Ancestry, Smoking, Tumor_Purity, ecDNA_status)

#select rows in data "Sample" column that are in the columns of matrix
data2 <- data[data$Sample %in% colnames(matrix),]
data3 <- merge(data2, ecDNA, by="Tumor_Barcode")

data4 <- data3 %>% dplyr::select(Tumor_Barcode, Sample)
covdata <- merge(metadata, data4, by="Tumor_Barcode")

#subset columns of matrix that are in the data "Sample" column
matrix2 <- matrix[,colnames(matrix) %in% data3$Sample]
rownames(covdata) <- covdata$Sample

#filter for nonsmokers
covdata <- covdata %>% filter(Smoking == "Smoker")

#remove the Sample and Tumor_Barcode columns from covdata
covdata <- covdata %>% dplyr::select(-Sample, -Tumor_Barcode, -Smoking)

#remove columns from matrix2 that are not in covdata
matrix2 <- matrix2[,rownames(covdata)]

#align the row names of covdata with the column names of matrix2
covdata2 <- covdata[match(colnames(matrix2), rownames(covdata)), , drop = FALSE]

#write.table(covdata2, "/Users/khandekara2/iCloud/dev/Sherlock-Lung/results2/Sherlock_EAGLE_RNAseq_covdata.txt", row.names = TRUE, quote = FALSE, sep = "\t")

covdata2$Gender = relevel(factor(covdata$Gender), ref = 'Female')
covdata2$Histology = relevel(factor(covdata$Histology), ref = 'Adenocarcinoma')
covdata2$Ancestry = relevel(factor(covdata$Ancestry), ref = 'EUR')
covdata2$ecDNA_status <- factor(covdata$ecDNA_status, levels = c("0", "1"))
#covdata2$Smoking = relevel(factor(covdata$Smoking), ref = 'Non-Smoker')

#subset the raw counts matrix to only include the samples that are in the matrix
raw_counts2 <- raw_counts[,colnames(raw_counts) %in% colnames(matrix2)]

#align the row names of covdata with the column names of matrix2

# Get gene annotations
gene_ids <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
genes <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_ids, 
                               columns = c("ENSEMBL", "SYMBOL", "GENETYPE"), 
                               keytype = "ENSEMBL")

genes <- data.frame(gene = rownames(raw_counts2)) %>%
  left_join(genes %>% as.data.frame, 
            by = c("gene"="SYMBOL")) %>%
  dplyr::distinct(gene, .keep_all=TRUE) 

rownames(genes) <- genes$gene

dge <- DGEList(counts = raw_counts2, 
               genes = genes,
               samples = covdata2,
               remove.zeros = TRUE)

# Perform PCA analysis:
pca_analysis <- prcomp(t(cpm(dge, log=TRUE)))
#summary(pca_analysis)

pca_plot <- ggbiplot::ggbiplot(pca_analysis, 
                               groups = dge$samples$ecDNA_status, 
                               ellipse = TRUE,
                               var.axes = FALSE)

print(pca_plot)

keepTheseGenes <- (rowSums(cpm(dge) > 2) >= 50) 

beforeFiltering_plot <- dge %>% 
  cpm(log = TRUE) %>% 
  melt %>% 
  dplyr::filter(is.finite(value)) %>% 
  ggplot(aes(x = value, colour = Var2)) +
  geom_density() + 
  guides(colour = FALSE) +
  ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
  labs(x = "logCPM", y = "Density")

afterFiltering_plot <- dge %>% 
  cpm(log = TRUE) %>% 
  magrittr::extract(keepTheseGenes,) %>%
  melt %>% 
  dplyr::filter(is.finite(value)) %>% 
  ggplot(aes(x = value, colour = Var2)) +
  geom_density() + 
  guides(colour = FALSE) +
  ggtitle("B. After filtering", subtitle = paste0(table(keepTheseGenes)[[2]], " genes"))+
  labs(x = "logCPM", y = "Density")

dge <- dge[keepTheseGenes,,keep.lib.sizes = FALSE] 

design <- model.matrix(~ 0 + ecDNA_status + Gender + Age + Histology + Tumor_Purity, data = covdata2)

voomData <- voom(dge, design = design, plot = FALSE)

# Make the column names of the design matrix syntactically valid
colnames(design) <- make.names(colnames(design))

contrasts <- makeContrasts(ecDNA_PosVsNeg = ecDNA_status1 - ecDNA_status0, levels = design)

fit <- lmFit(voomData, design) %>%
  contrasts.fit(contrasts) %>%
  treat(lfc = 1)

results <- decideTests(fit, 
                       p.value = 0.05,
                       adjust.method = "fdr") 


fit2 <- eBayes(fit)

results2 <- decideTests(fit2, 
                       p.value = 0.05,
                       adjust.method = "fdr") 

allDEresults <- topTable(fit2, 
                         coef = "ecDNA_PosVsNeg", 
                         number = Inf, 
                         adjust.method = "fdr") %>% as.data.frame()

allDEresults <- allDEresults %>%
  dplyr::mutate(isSignificant = case_when(
    abs(logFC) > 1.5 ~ TRUE, 
    TRUE ~ FALSE # If conditions in the line above are not met, gene is not DE. 
  ))

sigDEresults <- allDEresults %>%
  dplyr::filter(isSignificant == TRUE)

#order the results by logFC (highest to lowest)
allDEresults <- allDEresults[order(-allDEresults$logFC),]

#write out the results
write.table(allDEresults, "~/iCloud/dev/Sherlock-Lung/results2/log-reg-results2/smoker_RNAseq_DE_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")

volcano_plot <- allDEresults %>%
  ggplot(aes(x = logFC, 
             y = -log10(P.Value ),
             colour = isSignificant)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values = c("grey", "red")) +
  ggtitle("Never-Smokers RNA-seq Differential Expression Analysis",
          subtitle = "ecDNA Positive vs Negative") 

volcano_plot

smokers <- read.table("~/iCloud/dev/Sherlock-Lung/results2/log-reg-results2/smoker_RNAseq_DE_results.txt", header = TRUE, sep = "\t")

smoker_plot <- smokers %>%
  ggplot(aes(x = logFC, 
             y = -log10(P.Value ),
             colour = isSignificant)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_colour_manual(values = c("grey", "red")) +
  ggtitle("Smokers RNA-seq Differential Expression Analysis",
          subtitle = "ecDNA Positive vs Negative")


# Gene set enrichment analysis
library(clusterProfiler)
gene_list <- results$logFC
names(gene_list) <- rownames(results)
gene_list <- sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, ont="BP", keyType="ENSEMBL", minGSSize=3, maxGSSize=800, pvalueCutoff=0.05, verbose=FALSE, OrgDb=org.Hs.eg.db, pAdjustMethod="BH")

