# Assignment Script

# Davids Jalisevs, 12/12/2023

####################### Step 1
# Change working directory.
setwd("C:/Users/jalis/Downloads")
path  = "C:/Users/jalis/Downloads" # change this to your own directory

file_name = "brca_tcga_pan_can_atlas_2018.tar.gz"

file.exists("brca_tcga_pan_can_atlas_2018.tar.gz")

# We will first extract the files into folders.

####################### Step 2

untar(file_name)

# change directory to the extracted folders

setwd(paste(getwd() , "/brca_tcga_pan_can_atlas_2018", sep = ""))

####################### Step 3 and Step 4 and Step 5

# We will use the following files:

# Read clinical, RNASeq, and CNA data
clinical <- read.delim("data_clinical_patient.txt")
rnaseq <- read.delim("data_mrna_seq_v2_rsem.txt")
cna <- read.delim("data_cna.txt")

# Remove duplicated genes from RNASeq data
keep <- !duplicated(rnaseq[, 1])
rnaseq <- rnaseq[keep, ]
rownames(rnaseq) <- rnaseq[, 1]

# Find ERBB2 in CNA data and plot histogram
erbb2_indx <- which(cna[, 1] == 'ERBB2')
hist(as.numeric(cna[erbb2_indx, -c(1, 2)]))

####################### Step 6

# Match patients in RNASeq to patients in CNA
rna_cna_id <- which(is.element(colnames(rnaseq[, -c(1, 2)]), colnames(cna[, -c(1, 2)])))
rna_cna_sub <- rnaseq[, 2 + rna_cna_id]



# Check patients in RNASeq subset are in CNA
no_pats_in_rna_cna_sub_and_cna <- sum(is.element(colnames(rnaseq[, 2 + rna_cna_id]), colnames(cna[, -c(1, 2)])))
sanity_check <- no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

####################### Step 7

# Pre-allocate memory for ERBB2

meta_erbb2 = matrix(0,length(rna_cna_id),1)
for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
  
}



# simple checks to make sure. 

col_i = colnames(rna_cna_sub)[1]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.

pos_example = which(meta_erbb2==1)[1]


col_i = colnames(rna_cna_sub)[pos_example]

col_cna = which(colnames(cna)==col_i)

# sanity check

(cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# We will add a title to the metadata.

colnames(meta_erbb2) = 'ERBB2Amp'

# transform into integers

rna_cna_sub = round(rna_cna_sub)




library(DESeq2)
###########################################################################
# Step 8 : Normalising data 

# Create DESeqDataSet from count data
dds <- DESeqDataSetFromMatrix(countData = rna_cna_sub, colData = meta_erbb2, design = ~ ERBB2Amp)

# Run the DESeq normalization steps
dds <- DESeq(dds)

# Access normalized count data
normalized_counts <- counts(dds, normalized = TRUE)

#differential expression analysis
res <- results(dds)

#########################################################################
# Step 9: Obtaining Differentially expressed Genes

# Filter for significantly differentially expressed genes P value can be adjusted
sig_genes <- subset(res, padj < 0.05)

#top 10 differentially expressed genes ranked by fold change
top_genes <- head(sig_genes[order(-abs(sig_genes$log2FoldChange)), ], 10)

# Display the results
print("Differentially Expressed Genes:")
print(sig_genes)

print("Top 10 Differentially Expressed Genes:")
print(top_genes)

#######################################################################
# Step 10 :Pathaways Analysis
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

library(enrichplot)


# Assuming you have 'sig_genes' containing differentially expressed genes
de_genes <- sig_genes$gene_symbol

# Perform pathway enrichment analysis using clusterProfiler
library(clusterProfiler)

# Convert gene symbols to Entrez IDs using org.Hs.eg.db
gene_ids <- bitr(de_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform pathway enrichment analysis using KEGG
enrich_result <- enrichKEGG(gene         = gene_ids$ENTREZID,
                            organism     = 'hsa',
                            keyType      = 'kegg',
                            pAdjustMethod = 'BH', # You can choose the adjustment method
                            qvalueCutoff  = 0.05)  # Adjust this cutoff as needed

# Display the results
print(enrich_result)

# Plot the results
dotplot(enrich_result, showCategory = 15)  # Adjust 'showCategory' based on your preference



a########################################################################
# Step 11 :Get the variance stabilised transformed expression values.
#With the vst values obtain a PCA plot.

library(ggplot2)


# Get the variance-stabilized transformed expression values
vst_values <- vst(dds)

#Step 12 : With VST valuse create a PCA plot 
# Create a PCA plot
pca <- plotPCA(vst_values, intgroup = "ERBB2Amp", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
pca_plot <- ggplot(pca, aes(x = PC1, y = PC2, color = ERBB2Amp)) +
  geom_point(size = 3) +
  ggtitle("PCA Plot of VST Values") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"))

# Display the PCA plot
print(pca_plot)



####################################################
#
# Optional For EXTRA MARKS attempts
#
####################################################
# Perform clustering on vst values (using a subset of data)
subset_dist_values <- dist(t(assay(vst_values[, 1:10])))  # Adjust the subset size as needed
subset_cluster_values <- hclust(subset_dist_values)
subset_cluster_assignment <- cutree(subset_cluster_values, k = 4)  # Specify the number of clusters

# Create PCA plot with clustering information
pca_cluster_plot <- ggplot(pca, aes(x = PC1, y = PC2, color = factor(subset_cluster_assignment))) +
  geom_point(size = 3) +
  ggtitle("PCA Plot with Clustering") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"))

# Display the PCA plot with clustering
print(pca_cluster_plot)
####################################################
# Install the survival package if not already installed
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}
# Load the survival package
library(survival)

# Assuming Overall.Survival..Months. and Overall.Survival.Status are columns in your clinical data
surv_object <- with(clinical[5:nrow(clinical), ], {
  surv_time <- as.numeric(ifelse(is.na(as.character(Overall.Survival..Months.)), NA, as.character(Overall.Survival..Months.)))
  surv_event <- as.numeric(ifelse(as.character(Overall.Survival.Status) == "1", 1, 0))
  Surv(time = surv_time, event = surv_event)
})

# Ensure that there are non-zero observations
surv_object <- surv_object[surv_object$time > 0]

# Filter vst values for matching IDs
matching_ids <- rownames(vst_values)
vst_data <- assay(vst_values)[matching_ids %in% matching_ids, ]

# Ensure the number of rows matches
if (length(surv_object$time) != nrow(vst_data)) {
  stop("Number of rows in surv_object and vst_data do not match.")
}

# Combine survival object with vst values
surv_data <- data.frame(Survival = surv_object, vst_data)



