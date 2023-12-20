# Breast Cancer RNASeq Data Analysis

## Overview

This repository contains the code and analysis for the assignment on breast cancer RNASeq data. The goal of this assignment is to explore and analyze publicly available TCGA RNASeq data for breast cancer, specifically focusing on the ERBB2+ subtype.

## Dataset

The dataset used for this analysis can be found [here](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018). After downloading and extracting the files, we utilized three main data sources:

1. RNASeq Data: `data_mrna_seq_v2_rsem.txt`
2. Patient Data: `data_clinical_patient.txt`
3. Copy Number Aberrations Data: `data_cna.txt`

## Analysis Steps

1. **Data Preprocessing:**
   - Matching patient IDs between RNASeq, Patient Data, and Copy Number Aberrations Data.
   - Creating metadata based on ERBB2 amplification from CNA data.

2. **Normalization:**
   - Using DESeq2 to normalize the RNASeq data.

3. **Differential Expression Analysis:**
   - Identifying differentially expressed genes between ERBB2+ and other breast cancer tumors.

4. **Pathway Enrichment Analysis:**
   - Analyzing enriched pathways among the differentially expressed genes.

5. **PCA Plot:**
   - Generating a PCA plot using the variance-stabilized transformed expression values.

6. **Optional Analysis (Additional Marks):**
   - Gene Expression Cluster
   - Cox Regression Model with DE Genes
   - Lasso Cross-Validation for Gene Selection in Survival Prediction

## Code Structure

The primary code for this assignment is available in the `Assignment_Code.R` file. Additionally, the complete code is provided in a separate text document named `Assignment_Code_Full.txt`.

## Report

For a detailed analysis report, please refer to the report document.

## Instructions for Reproduction

1. Download the dataset from the provided link.
2. Untar the folder and extract the files.
3. Execute the R code in `Assignment_Code.R` for step-by-step analysis.

Feel free to explore the code and adapt it to your needs!

---

*Note: Detailed analysis results and additional discussions are available in the full report document.*
