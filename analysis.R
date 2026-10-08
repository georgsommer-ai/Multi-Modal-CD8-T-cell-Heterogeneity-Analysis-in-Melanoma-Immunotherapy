# ==============================================================================
# Workflow: Mahuron / Pauken et al. 2025 - Melanom
# Re-Implementation by Georg Sommer
# ==============================================================================


# Load libraries ----
# ------------------------------------------------------------------------------

library(Seurat)
library(SeuratData)
library(tidyverse)
library(DESeq2)
library(tximport)
library(fgsea)
library(TITAN)
library(pdist)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)
library(Azimuth)

set.seed(42)


# SETUP PROJECTFOLDER ----
# ------------------------------------------------------------------------------

dirs <- c("data/bulk_raw", "data/sc_raw", "plots")
for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}



# ==============================================================================
# A. Bulk RNA-Seq Workflow ----
# ==============================================================================

# download Bulk RNA-seq Dataset ---- 
bulk_dir <- "data/bulk_raw/GSE147620"

if (!dir.exists(bulk_dir)) {
  GEOquery::getGEOSuppFiles("GSE147620", baseDir = "data/bulk_raw", makeDirectory = TRUE)
  tars <- list.files(bulk_dir, pattern = "\\.tar$", full.names = TRUE)
  for (t in tars) {
    message("Unpacking ", t, " ...")
    untar(t, exdir = bulk_dir)
  }
} else {
  message("Bulk-Data exist in: ", bulk_dir)
}


# prepare bulk RNA-seq data ----

bulk_files <- list.files(bulk_dir, pattern = "\\.(txt|tsv|csv)(\\.gz)?$", full.names = TRUE)

dat <- readr::read_delim(bulk_files[1], show_col_types = TRUE)

dat <- as.data.frame(dat)

# Remove colname from index col for DESeq2
rownames(dat) <- dat[, 1]
counts <- dat[, -1, drop = FALSE]

num_samples <- ncol(counts)


# 4.1 Differential Expression with DESeq2 ----
# ------------------------------------------------------------------------------

# Prep Metadata for DESeq2 
correct_condition <- ifelse(grepl("PD-1hi", colnames(counts)), "CPHi", "CPLo")

correct_patient <- gsub("melanoma|_.*", "", colnames(counts))

sample_info <- data.frame(
  condition = factor(correct_condition, levels = c("CPLo", "CPHi")), 
  patient = as.factor(correct_patient),
  row.names = colnames(counts)
)
sample_info

# Create DESeq2-Object 
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = sample_info,
  design = ~ patient + condition
)

# DESeq2-analysis
dds <- DESeq(dds)

# Compare CPHi vs CPLo
res <- results(
  dds,
  contrast = c("condition", "CPHi", "CPLo")
)
res


# 4.2 Ranking and CP^Hi and CP^Lo Gene Signatures  ----
# ------------------------------------------------------------------------------

res_df <- as.data.frame(res) %>% drop_na()

# VOLCANO PLOT to show CPHi gene signature
p_bulk_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & log2FoldChange > 0), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Bulk RNA-seq (CPHi Signatur)", x = "log2 Fold Change", y = "-log10(padj)") +
  theme(legend.position = "none")

p_bulk_volcano

ggsave("plots/4_2_Bulk_Volcano_Plot.png", plot = p_bulk_volcano, width = 9, height = 9)

# copy for CPLo signature
res_df_cplo <- res_df

# Ranking & gene signatures CPHi
res_df <- res_df %>%
  filter(padj < 0.05, log2FoldChange > 0) %>% 
  mutate(
    rank_padj = rank(padj),
    rank_log2FC = rank(-log2FoldChange), # Minus:most upregulated Gen den Platz 1
    average_rank_CPHi = (rank_padj + rank_log2FC) / 2
  ) %>%
  arrange(average_rank_CPHi)

CPHi_sig <- head(rownames(res_df), 200)
CPHi_sig
# only 46 CPHi genes
length(CPHi_sig)

# Ranking & gene signatures CPLo:
res_df_cplo <- res_df_cplo %>%
  filter(padj < 0.05, log2FoldChange < 0) %>% 
  mutate(
    rank_padj = rank(padj),
    rank_log2FC = rank(log2FoldChange), # Kein Minus! Kleinste (negativste) Werte bekommen Rang 1
    average_rank_CPLo = (rank_padj + rank_log2FC) / 2
  ) %>%
  arrange(average_rank_CPLo)

CPLo_sig <- head(rownames(res_df_cplo), 200)
CPLo_sig
length(CPLo_sig)



# ==============================================================================
# B. SINGLE-CELL RNA-SEQ WORKFLOW ----
# ==============================================================================

# download Single-Cell RNA-seq data----
sc_dir <- "data/sc_raw/GSE148190"
if (!dir.exists(sc_dir)) {
  GEOquery::getGEOSuppFiles("GSE148190", baseDir = "data/sc_raw", makeDirectory = TRUE)
  
  tars <- list.files(sc_dir, pattern = "\\.tar$", full.names = TRUE)
  for (t in tars) {
    message("Unpacking ", t, " ...")
    untar(t, exdir = sc_dir)
  }
} else {
  message("Single-Cell-Data exist in: ", sc_dir)
}

# quick summary of SC-files:
list.files(sc_dir)

# 6 samples from 3 patients :

# Patient - Tissue:
# K383 lymph node (aka pilot2)

# K409 blood
# K409 lymph node
# K409 tumor
# K409 VDJ-files

# K411 blood
# K411 lymph node



# 4.5 Create Seurat object ----
# -------------------------------------------------

mtx_files <- list.files(sc_dir, pattern = "matrix.mtx(\\.gz)?$", full.names = TRUE, recursive = TRUE)
mtx_files

seurat_list <- list()

for (mtx in mtx_files) {
  
  message("1. current matrix-file (mtx): \n   ", mtx)
  pat_id <- gsub("_?matrix.mtx(\\.gz)?$", "", basename(mtx))
  
  # Fallback if ID empty
  if (pat_id == "") pat_id <- basename(dirname(mtx))
  message("\n2. Extract Patient-ID (pat_id): \n   ", pat_id)
  
  b_file <- list.files(dirname(mtx), pattern = paste0("^(", pat_id, "_)?barcodes.tsv(\\.gz)?$"), full.names = TRUE)[1]
  message("\n3. Found Barcode-file (b_file): \n   ", b_file)
  
  f_file <- list.files(dirname(mtx), pattern = paste0("^(", pat_id, "_)?(features|genes).tsv(\\.gz)?$"), full.names = TRUE)[1]
  message("\n4. Found Feature-file (f_file): \n   ", f_file)
  
  message("\n5. Creating Seurat-object for ", pat_id, "...")
  seurat_list[[pat_id]] <- CreateSeuratObject(counts = ReadMtx(mtx, b_file, f_file), project = pat_id)
  
  seurat_list[[pat_id]]$patient_id <- pat_id
  
  message("Sucessfully completed for : ", pat_id)
}

sc_obj <- if(length(seurat_list) > 1) merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list)) else seurat_list[[1]]

# EDA:
sc_obj
# meta data
head(sc_obj@meta.data, 5)
# countmatrix first 5 rows, 2 columns
LayerData(sc_obj[["RNA"]], layer = "counts.GSM4455931_pilot2_GEX")[1:5, 1:2]
      



# 4.6 Quality control ----
# ------------------------------------------------------------------------------
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
head(sc_obj@meta.data, 10)

# --- PLOT: EDA BEFORE QC-FILTERING ---
p_qc_vln <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave("plots/4_6_QC_VlnPlot_PreFilter.png", plot = p_qc_vln, width = 10, height = 5)
p_qc_vln

# Explanation for 5% cutoff: 
# 5% cutoff for mitochondrial RNA is what the authors did, and was kept to compare results here with theirs
# In solid tumor cells, the metabolism is massively upregulated and the mitochondrial proportion correlates with this, reaching even 20%. 
# In contrast, healthy T cells shouldn't get high mitochondrial RNA percentages so we can cut off at 5% without too much data loss
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & percent.mt < 5)

# --- PLOT: EDA AFTER QC-FILTERING ---
p_qc_vln2 <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave("plots/4_6_QC_VlnPlot_Post_Filter.png", plot = p_qc_vln2, width = 24, height = 13, limitsize = FALSE)
p_qc_vln2

# meta data
head(sc_obj@meta.data, 5)



# 4.7 SCT- Dim. Reduction - Clustering ----
# ------------------------------------------------------------------------------

sc_obj <- SCTransform(sc_obj, verbose = FALSE)
sc_obj <- RunPCA(sc_obj, verbose = FALSE)
sc_obj <- RunUMAP(sc_obj, dims = 1:30, verbose = FALSE) 
sc_obj <- FindNeighbors(sc_obj, dims = 1:30, verbose = FALSE)
# Adjusted resolution to get same nr of clusters as original study for comparable results 
# note: different UMAP due to bigger cohort of original study 
sc_obj <- FindClusters(sc_obj, resolution = 0.12, verbose = FALSE) # 0.8

# UMAP pre-integration ---
p_umap_pre <- DimPlot(sc_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP Clusters (Pre-Integration)")
ggsave("plots/4_7_UMAP_PreIntegration.png", plot = p_umap_pre, width = 5, height = 5, limitsize = FALSE)
p_umap_pre



# 4.8 Removing TCR-genes before integration ----
# ------------------------------------------------------------------------------
# clustering should be only influenced by the cell state 
tcr_genes <- grep("^TR[ABGD][VJC]", rownames(sc_obj), value = TRUE)
keep_genes <- setdiff(rownames(sc_obj), tcr_genes)
sc_obj <- subset(sc_obj, features = keep_genes)

# EDA
tcr_genes
head(keep_genes, 20)
length(tcr_genes)
length(keep_genes)



# 4.9 Integration (CCA) ----
# ------------------------------------------------------------------------------
obj_list <- SplitObject(sc_obj, split.by = "patient_id")
obj_list
obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:20)
sc_obj_int <- IntegrateData(anchorset = anchors, dims = 1:20)



# 4.10 Normalization & UMAP Comparison ----
# ------------------------------------------------------------------------------
DefaultAssay(sc_obj_int) <- "RNA"

sc_obj_int <- JoinLayers(sc_obj_int)
sc_obj_int <- NormalizeData(sc_obj_int, normalization.method = "LogNormalize", scale.factor = 10000)
sc_obj_int <- FindVariableFeatures(sc_obj_int, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc_obj_int <- ScaleData(sc_obj_int, verbose = FALSE)
sc_obj_int <- RunPCA(sc_obj_int, verbose = FALSE)
sc_obj_int <- RunUMAP(sc_obj_int, dims = 1:20, verbose = FALSE)

sc_obj_int <- FindNeighbors(sc_obj_int, dims = 1:20, verbose = FALSE)
sc_obj_int <- FindClusters(sc_obj_int, resolution = 0.19, verbose = FALSE)

# EDA: 
sc_obj_int
head(sc_obj_int@meta.data, 5)

# PLOT: UMAP POST-INTEGRATION
p_umap_int <- DimPlot(sc_obj_int, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP Clusters (Post-Integration)")
ggsave("plots/4_10_UMAP_PostIntegration.png", plot = p_umap_int, width = 7, height = 5)
p_umap_int

# UMAP cluster original study
orig_data <- readRDS(
   "data/sc_raw/GSE148190/orig_data/GSE148190_HumanTumorIntegrated_filtered.rds")
orig_data

# EDA
orig_data
head(orig_data@meta.data, 10)

# Samples:
unique(orig_data@meta.data$sample)

# PLOT:orig UMAP Post-INTEGRATION
p_umap_post_orig <- DimPlot(orig_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP Clusters (Post-Integration_orig)")
ggsave("plots/4_10_UMAP_Post_Integration_orig.png", plot = p_umap_post_orig, width = 7, height = 5)
p_umap_post_orig
p_umap_int | p_umap_post_orig

# Tissue / Patient Source on UMAP ----
# ------------------------------------------------------------------------------

# add source to metadata 
sc_obj_int$tissue <- "Unknown"
sc_obj_int$tissue[grepl("GSM4455931", sc_obj_int$orig.ident)] <- "Lymph Node (K383)"
sc_obj_int$tissue[grepl("GSM4455932", sc_obj_int$orig.ident)] <- "Blood (K409)"
sc_obj_int$tissue[grepl("GSM4455933", sc_obj_int$orig.ident)] <- "Lymph Node (K409)"
sc_obj_int$tissue[grepl("GSM4455935", sc_obj_int$orig.ident)] <- "Tumor (K409)"
sc_obj_int$tissue[grepl("GSM4455937", sc_obj_int$orig.ident)] <- "Blood (K411)"
sc_obj_int$tissue[grepl("GSM4455938", sc_obj_int$orig.ident)] <- "Lymph Node (K411)"

p_tissue <- DimPlot(sc_obj_int, reduction = "umap", group.by = "tissue") + 
  ggtitle("Tissue / Patient source")
p_int_check <- p_umap_int | p_tissue

ggsave("plots/4_10_p_int_check.png", plot = p_int_check, width = 16, height = 9, limitsize = FALSE)
p_int_check



# 4.11 Comparability of the results with the original study using a confusionmatrix and UMAPs of Commonalities
# ------------------------------------------------------------------------------
# which cluster of our UMAP corresponds with the cluster of the original study
# note: original study has bigger cohort 

# get cell barcodes from rownames
cells <- rownames(sc_obj_int@meta.data)
cells_orig <- rownames(orig_data@meta.data)
head(cells, 10)

# isolate barcodes
barcodes <- str_extract(rownames(sc_obj_int@meta.data), "[A-Z]{16}")
barcodes_orig <- str_extract(rownames(orig_data@meta.data), "[A-Z]{16}")

# match barcodes
match_idx <- match(barcodes, barcodes_orig)
match_idx

# prepare table for confusion matrix
comparison_table <- data.frame(
cell = cells,
barcode = barcodes,
cell_orig = cells_orig[match_idx], 
cluster_orig = orig_data$seurat_clusters[match_idx] 
)
# add my corresponding clusters as column my_clusters 
comparison_table$my_cluster <- sc_obj_int$seurat_clusters
head(comparison_table, 15)

# remove not matching cells (na)
comparison_table <- comparison_table[!is.na(comparison_table$cluster_orig), ]
head(comparison_table, 15)

confusion_matrix <- round(
prop.table(
table(
my_cluster   = comparison_table$my_cluster,
cluster_orig = comparison_table$cluster_orig
),
margin = 1
) * 100,
1
)
confusion_matrix

# Heatmap of confusion_matrix
conf_df <- as.data.frame(confusion_matrix)

p_heatmap <- ggplot(conf_df, aes(x = cluster_orig, y = my_cluster, fill = Freq)) +
geom_tile(color = "grey90", size = 0.5) +
geom_text(aes(label = ifelse(Freq > 0, Freq, "")), color = "black", size = 4) +
scale_fill_gradient(low = "white", high = "firebrick3", name = "Percent") +
theme_minimal() +
theme(
  axis.text = element_text(size = 12, face = "bold"),
  axis.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  panel.grid = element_blank() # Verhindert störende Gitterlinien hinter den Kacheln
) +
labs(
  title = "Confusion Matrix: UMAP clusters from this work vs original study",
  x = "orig. study clusters",
  y = "my clusters"
)
p_heatmap

# Inspect matching clusters 
Idents(sc_obj_int) <- "seurat_clusters"
Idents(orig_data) <- "seurat_clusters"

p_my_umap <- DimPlot(sc_obj_int, reduction = "umap", label = TRUE, label.size = 5) + 
  ggtitle("DEINE UMAP (sc_obj_int)") +
  theme(legend.position = "none") 

p_umap_orig <- DimPlot(orig_data, reduction = "umap", label = TRUE, label.size = 5) + 
  ggtitle("AUTOREN UMAP (orig_data)") +
  theme(legend.position = "none")

p_my_umap | p_umap_orig

sc_obj_int$Shared_Color <- "no clear match"
orig_data$Shared_Color <- "no clear match"

my_cluster <- as.character(sc_obj_int$seurat_clusters)
orig_cluster <- as.character(orig_data$seurat_clusters)

sc_obj_int$Shared_Color[my_cluster == "11"] <- "Match: 11 = 9"
orig_data$Shared_Color[orig_cluster == "9"]  <- "Match: 11 = 9"

sc_obj_int$Shared_Color[my_cluster == "0"] <- "Match: 0 = 1"
orig_data$Shared_Color[orig_cluster == "1"] <- "Match: 0 = 1"

sc_obj_int$Shared_Color[my_cluster == "8"] <- "Match: 8 = 3"
orig_data$Shared_Color[orig_cluster == "3"] <- "Match: 8 = 3"

sc_obj_int$Shared_Color[my_cluster == "2"] <- "Match: 2 = 0,2,6"
orig_data$Shared_Color[orig_cluster %in% c("0", "2", "6")] <- "Match: 2 = 0,2,6"

color_palette <- c(
  "Match: 11 = 9" = "#E31A1C", # Kräftiges Rot
  "Match: 0 = 1"  = "#1F78B4", # Kräftiges Blau
  "Match: 8 = 3"  = "#33A02C", # Kräftiges Grün
  "Match: 2 = 0,2,6"  = "#6A3D9A", # Kräftiges Lila
  "no clear match"   = "grey85"   # Dezentes Hellgrau
)

p_similarity_my <- DimPlot(sc_obj_int, group.by = "Shared_Color", pt.size = 0.5) +
  scale_color_manual(values = color_palette) +
  ggtitle("my UMAP") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) + NoLegend()

p_similarity_orig <- DimPlot(orig_data, group.by = "Shared_Color", pt.size = 0.5) +
  scale_color_manual(values = color_palette) +
  ggtitle("orig. study UMAP") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_main1 <- p_umap_int | p_umap_post_orig | p_similarity_my | p_similarity_orig
p_main1b <- p_umap_int | p_umap_post_orig
p_similarity <- p_heatmap | (p_similarity_my | p_similarity_orig + plot_layout(guides = "collect"))

ggsave("plots/main1.png", plot = p_main1, width = 32, height = 9, limitsize = FALSE)
ggsave("plots/p_main1b.png", plot = p_main1b, width = 16, height = 9, limitsize = FALSE)
ggsave("plots/p_similarity.png", plot = p_similarity, width = 24, height = 9, limitsize = FALSE)
p_main1
p_similarity



# 4.12. Gene-signature scoring w. ModuleScore / Differential expression ----
# ------------------------------------------------------------------------------

# convert Ensembl-IDs to Gene-Symbols
cphi_symbols <- mapIds(
org.Hs.eg.db,
keys = CPHi_sig,           
column = "SYMBOL",         
keytype = "ENSEMBL",       
multiVals = "first"        
)
head(cphi_symbols, 20)

# delete Ensembl-IDs with no Gene-Symbols (NA)
cphi_sig_ena <- as.character(na.omit(unname(cphi_symbols)))
cphi_sig_ena

message(sprintf("Successful conversions: %s of %s genes.", 
              length(cphi_sig_ena), length(CPHi_sig)))

sc_obj_int <- AddModuleScore(
  object = sc_obj_int, 
  features = list(cphi_sig_ena), 
  name = "CPHi_Score",
  ctrl = 100,  # Standard 
  nbin = 24    # Standard
)

# ModuleScore/CPHi_Score added
head(sc_obj_int@meta.data, 10)

# FeaturePlot 
p_score_feature <- FeaturePlot(sc_obj_int, features = "CPHi_Score1", cols = c("lightblue", "darkred")) +
                   ggtitle("CPHi Gene Signature score: CPHi_Score1")

ggsave("plots/4_12_CPHi_ModuleScore200_FeaturePlot.png", plot = p_score_feature, width = 7, height = 5)
p_score_feature



# 4.13 Methods for CPHi/CPLo-classifikation ----
# ------------------------------------------------------------------------------

# METHOD 1: CPHi signature Score based on Module Score (mean & SD) ---- 
mean_score <- mean(sc_obj_int$CPHi_Score1)
sd_score <- sd(sc_obj_int$CPHi_Score1)

# Method 1A: Strict Threshold (Mean + 1 Standard Deviation)
threshold_1_sd <- mean_score + (1 * sd_score)
sc_obj_int$CP_Class_Method2_Sig_1SD <- ifelse(
  sc_obj_int$CPHi_Score1 > threshold_1_sd, 
  "CPHi", 
  "CPLo"
)

# Method 1B: relaxed Threshold (Mean + 0.25 Standard Deviation)
threshold_025_sd <- mean_score + (0.25 * sd_score)
sc_obj_int$CP_Class_Method2_Sig_025SD <- ifelse(
  sc_obj_int$CPHi_Score1 > threshold_025_sd, 
  "CPHi", 
  "CPLo"
)


# METHOD 2: Classification by avg Module Score by Cluster ----
cluster_mean_scores <- tapply(sc_obj_int$CPHi_Score1, sc_obj_int$seurat_clusters, mean)

# CPHi_Score1 by cell
head(sc_obj_int$CPHi_Score1, 20)
# Cluster nr by cell
head(sc_obj_int$seurat_clusters, 20)
cluster_mean_scores

# Top 5( corresponds to a cutoff of -0.01)
cphi_clusters <- names(head(sort(cluster_mean_scores, decreasing = TRUE), 5))

cphi_clusters

sc_obj_int$CP_Class_Method3_Cluster <- ifelse(
  sc_obj_int$seurat_clusters %in% cphi_clusters, 
  "CPHi", 
  "CPLo"
)


# METHOD 3: PDCD1/CTLA4 Coexpression: only double-positive cells are "CPHi" ----
pdcd1_expr <- GetAssayData(sc_obj_int, assay="RNA", layer="data")["PDCD1", ]
ctla4_expr <- GetAssayData(sc_obj_int, assay="RNA", layer="data")["CTLA4", ]

# Visual Identification of Cutoffs (Natural Breaks)
# Natural Breaks lies at local minimum between the two peaks: 0.4
plot(density(pdcd1_expr), main = "Density Plot: PDCD1", col = "blue", lwd = 2)
plot(density(ctla4_expr), main = "Density Plot: CTLA4", col = "red", lwd = 2)

cutoff_pdcd1 <- 0.4
cutoff_ctla4 <- 0.4

sc_obj_int$CP_Class_Method1_Direct <- ifelse(
  pdcd1_expr > cutoff_pdcd1 & ctla4_expr > cutoff_ctla4, 
  "CPHi", 
  "CPLo"
)


p_class_m2 <- DimPlot(sc_obj_int, reduction = "umap", group.by = "CP_Class_Method2_Sig_1SD", cols = c("darkred", "lightblue")) +
  ggtitle("Method 1: Module Score \n (Basis: mean & SD)")

p_class_m3 <- DimPlot(sc_obj_int, reduction = "umap", group.by = "CP_Class_Method3_Cluster", cols = c("darkred", "lightblue")) +
  ggtitle("Method 2: avg Module Score \n on Cluster-level")

p_class_m1 <- DimPlot(sc_obj_int, reduction = "umap", group.by = "CP_Class_Method1_Direct", cols = c("darkred", "lightblue")) +
  ggtitle("Method 3: Direct PDCD1/CTLA4 \n Coexpression")

p_combined <-  p_class_m2 | p_class_m3 | p_class_m1
ggsave("plots/4_13_CP_Class_UMAP_All_Methods.png", plot = p_combined, width = 15, height = 5)
p_combined



# 4.14 Inhibitory-receptor coexpression: SUM of "PDCD1", "CTLA4", "HAVCR2", "LAG3" "TIGIT"----
# ------------------------------------------------------------------------------
ir_genes <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT")

# Extract expression matrix for ir_genes 
ir_matrix <- GetAssayData(sc_obj_int, assay="RNA", layer="data")[ir_genes, ]

# create score by summing (normalised) counts 
sc_obj_int$IR_Score <- colSums(ir_matrix)

# COEXPRESSION FEATURE PLOT
p19_feature <- FeaturePlot(sc_obj_int, features = "IR_Score", cols = c("lightgrey", "darkred"), pt.size = 1) +
               ggtitle("Inhibitory-Receptor Coexpression \n PDCD1, CTLA4, HAVCR2, LAG3 TIGIT")
p_4methods_cphi <- p_combined | p19_feature
p_4methods_cphi
# VIOLIN PLOT 
p19_vln <- VlnPlot(sc_obj_int, 
                   features = ir_genes, 
                   group.by = "CP_Class_Method2_Sig_1SD", 
                   pt.size = 0, # pt.size = 0 blendet die störenden schwarzen Punkte im Violinplot aus
                   ncol = 3) # Ordnet die 5 Plots schöner an

p19_vln

ggsave("plots/4_14_4methods_cphi.png", plot = p_4methods_cphi, width = 32, height = 9, limitsize = FALSE)
ggsave("plots/4_14_IR_Coexpression_VlnPlot.png", plot = p19_vln, width = 12, height = 8)



# 4.15 TCR clonotype analysis ----
# ------------------------------------------------------------------------------
# LOAD CONTIGS (Barcode <-> Klonotype-ID)
ln_contigs <- read.csv(file.path(sc_dir, "GSM4455934_K409_LN_VDJ_filtered_contig_annotations.csv.gz"))
tumor_contigs <- read.csv(file.path(sc_dir, "GSM4455936_K409_tumor_VDJ_filtered_contig_annotations.csv.gz"))
head(ln_contigs,10)
head(tumor_contigs,10)

# LOAD CLONOTYPES (Klonotyp-ID <-> Frequency/CloneSize)
ln_clonotypes <- read.csv(file.path(sc_dir, "GSM4455934_K409_LN_VDJ_clonotypes.csv.gz"))
tumor_clonotypes <- read.csv(file.path(sc_dir, "GSM4455936_K409_tumor_VDJ_clonotypes.csv.gz"))
head(ln_clonotypes,10)
head(tumor_clonotypes,10)
      
ln_cells <- unique(ln_contigs[, c("barcode", "raw_clonotype_id")])
tumor_cells <- unique(tumor_contigs[, c("barcode", "raw_clonotype_id")])

# look up Frequency (clonesize) 
ln_cells <- merge(ln_cells, ln_clonotypes[, c("clonotype_id", "frequency")], 
                  by.x = "raw_clonotype_id", by.y = "clonotype_id", all.x = TRUE)
tumor_cells <- merge(tumor_cells, tumor_clonotypes[, c("clonotype_id", "frequency")], 
                     by.x = "raw_clonotype_id", by.y = "clonotype_id", all.x = TRUE)
head(ln_cells,10)
head(tumor_cells,10)

# join LN & Tumor
all_tcr_cells <- rbind(ln_cells, tumor_cells)

# CREATE VEKTOR FOR SEURAT (cell barcode, frequency)
clone_sizes <- all_tcr_cells$frequency
names(clone_sizes) <- all_tcr_cells$barcode
head(clone_sizes, 20)

# BARCODE-MAPPING 
seurat_barcodes <- rownames(sc_obj_int@meta.data)
head(seurat_barcodes,10)

seurat_raw_barcodes <- sub(".*_", "", seurat_barcodes)
head(seurat_raw_barcodes, 10)

clone_sizes_aligned <- clone_sizes[match(seurat_raw_barcodes, names(clone_sizes))]
head(clone_sizes_aligned, 50)

names(clone_sizes_aligned) <- seurat_barcodes
head(clone_sizes_aligned, 50)

sc_obj_int <- AddMetaData(object = sc_obj_int, metadata = clone_sizes_aligned, col.name = "CloneSize")
head(sc_obj_int@meta.data, 20)

# Identify cloned cells: Singleton vs Expanded/Cloned
sc_obj_int$CloneSize[is.na(sc_obj_int$CloneSize)] <- 1
head(sc_obj_int@meta.data, 20)

sc_obj_int$Clonal_Expansion <- ifelse(sc_obj_int$CloneSize > 1, "Expanded/Cloned", "Singleton")

table(sc_obj_int$Clonal_Expansion)


p20_dim <- DimPlot(sc_obj_int, 
                   group.by = "Clonal_Expansion", 
                   reduction = "umap", 
                   cols = c("Expanded/Cloned" = "darkblue", "Singleton" = "lightgrey")) +
           ggtitle("TCR Clonal Expansion Mapping")

# CLONE SIZE BY EXHAUSTION STATE
p20_vln <- VlnPlot(sc_obj_int, 
                   features = "CloneSize", 
                   group.by = "CP_Class_Method2_Sig_1SD", 
                   pt.size = 0) + 
           ggtitle("Clonal Expansion by Exhaustion State")

p20_dim
p20_vln

# for readme
main2 <-p_class_m2 | p_class_m1 | p19_feature | p20_dim

dir.create("plots", showWarnings = FALSE)
ggsave("plots/4_15_TCR_Clonal_Expansion_UMAP.png", plot = p20_dim, width = 8, height = 9, limitsize = FALSE)
ggsave("plots/4_15_TCR_Clone_Size_Vln.png", plot = p20_vln, width = 6, height = 5)
ggsave("plots/main2.png", plot = main2, width = 32, height = 9, limitsize = FALSE)



# 4.16. TITAN Topic Modeling (Soft Clustering) ----
# ------------------------------------------------------------------------------

# 4.16.1: Compatibility Patch (Seurat V5 vs. TITAN)
# ------------------------------------------------------------------------------
# problem: in Seurat V5 "slot" was replaced with "layer" 
# solution: patch the TITAN function 'runLDA'

replace_slot_with_layer <- function(expr) {
  
  # Ist das aktuelle Textstück ein Befehl /ein 'call'
  if (!is.call(expr)) {
    return(expr)
  }
  
  # Ist es der Befehl 'GetAssayData', der den Fehler verursacht?
  if (identical(expr[[1]], as.name("GetAssayData"))) {
    
    # alle Argumente, die an diesen Befehl übergeben werden, in einer Variable speichern 
    argument_names <- names(expr)
    
    # Enthält der Befehl den fehlerhaften Begriff 'slot'?
    if (!is.null(argument_names) && "slot" %in% argument_names) {
      
      
      cat("\n Die Funktion 'GetAssayData' wurde gefunden!\n")
      cat("Status VOR der Reparatur (Die Argumente heißen):\n")
      print(argument_names) 
      
      # Reparatur: Wir suchen in der Liste der Argumente nach 'slot' und überschreiben es mit 'layer'.
      argument_names[argument_names == "slot"] <- "layer"
      
      cat("Status NACH der Reparatur (Wir haben 'slot' gelöscht und 'layer' eingefügt):\n")
      print(argument_names) # Hier drucken wir die reparierten Variablennamen aus
      cat("--------------------------------------------------\n")
      
      # Update code: reparierten Namen zurück in den eigentlichen Code-Befehl geben.
      names(expr) <- argument_names
    }
  }
  
  # Rekursion:
  for (i in seq_along(expr)) {
    
    # Die Funktion ruft sich selbst auf (Rekursion). Das bedeutet: Sie nimmt das
    # innere Element und schickt es wieder ganz nach oben zu SCHRITT 1.
    expr[[i]] <- replace_slot_with_layer(expr[[i]])
  }
  
  # Wir geben den komplett überprüften und (falls nötig) reparierten Code zurück.
  expr
}


# 4.16.2: Patch TITAN Function 
# ------------------------------------------------------------------------------

# create a copy  of TITAN-function 'runLDA' that is compatible with Seurat V5 (runLDA_v5).
runLDA_v5 <- TITAN::runLDA

body(runLDA_v5) <- replace_slot_with_layer(body(runLDA_v5))

environment(runLDA_v5) <- environment(TITAN::runLDA)


# 4.16.3: Run TITAN Topic Modeling 
# ------------------------------------------------------------------------------
LDA_model <- runLDA_v5(
  sc_obj_int, 
  ntopics = 3, 
  normalizationMethod = "CLR", 
  seed.number = 42, 
  assayName = "RNA" 
)

# Add Topic Weights to metadata
sc_obj_int <- TITAN::addTopicsToSeuratObject(model = LDA_model, Object = sc_obj_int)

# 4.16.4 Validate TITAN Topics
# ------------------------------------------------------------------------------
model_names <- names(LDA_model)
model_names
LDA_model$topics[,1:12]
topic_matrix <- LDA_model$topics

top_genes_list <- list()

num_topics<-nrow(LDA_model$topics)

# extract the Top120 gene counts by Topic
for (i in 1:num_topics) {
  
  current_weights <- topic_matrix[i, ]
  
  sorted_weights <- sort(current_weights, decreasing = TRUE)
  
  top_genes <- names(sorted_weights)[1:120]
  
  topic_name <- paste0("Topic_", i)
  
  top_genes_list[[topic_name]] <- top_genes
}

titan_top_genes <- as.data.frame(top_genes_list)

# the genes show:Topic 1: Progenitor, Topic 2: Glycolysis / Effector, Topic 3: Terminal Exhaustion
titan_top_genes

# TPEX marker: TCF7, SLAMF6 and others
marker <- c("TCF7", "SLAMF6", "IL7R", "LEF1",
            "BCL6", "ID3", "MYB", "XCL1")
idx <- match(marker, colnames(topic_matrix))

marker_check <- do.call(rbind, lapply(seq_len(nrow(topic_matrix)), function(i) {
  w <- topic_matrix[i, ]

  data.frame(
    gene = marker,
    Topic = paste0("Topic_", i),
    gene_in_data = !is.na(idx),
    rank = unname(rank(-w, ties.method = "min")[idx]),
    weight = unname((w / sum(w))[idx])
  )
}))

print(marker_check, row.names = FALSE, digits = 4)


# plots
p21_topics <- FeaturePlot(
  sc_obj_int, 
  features = c("Topic_1", "Topic_2", "Topic_3"), 
  ncol = 3, 
  cols = c("lightgrey", "purple") 
) 


p21_topics[[1]] <- p21_topics[[1]] + 
  labs(
    title = "Topic 1: Precursor / Progenitor T cells/ TPEX-like",
    color = "Topic 1 Weight"
  )


p21_topics[[2]] <- p21_topics[[2]] + 
  labs(
    title = "Topic 2: Metabolically Active/Activated/Glycolytic",
    color = "Topic 2 Weight"
  )

p21_topics[[3]] <- p21_topics[[3]] + 
  labs(
    title = "Topic 3: Cytotoxic / Effector / (Terminally) Exhausted",
    color = "Topic 3 Weight"
  )

ggsave("plots/main3_titan.png", plot = p21_topics, width = 24, height = 9)

p21_topics



# 4.17 Azimuth Classification with Tonsil-Reference
# ------------------------------------------------------------------------------
# show all available references
AvailableData()
# 
# An attempt at reference based cell type classification was made using Azimuth. The most appropriate references available (PBMC and the Tonsil) were tested. Neither effectively classified cells in the dataset.

options(timeout = 900)
InstallData("tonsilref")

sc_obj_int <- RunAzimuth(sc_obj_int, reference = "tonsilref")

# Azimuth may leave automatic QC metric calculation disabled (False).
getOption("Seurat.object.assay.calcn")
# We reset this option to its default (NULL) to prevent missing QC metrics in subsequent analyses.
options(Seurat.object.assay.calcn = NULL)

# Azimuth tonsilref clusters
p_azimuth <- DimPlot(sc_obj_int, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE, label.size = 3) + 
  ggtitle("Azimuth (Tonsilref L2)") +
  theme(legend.position = "none")

ggsave("plots/4_17_azimuth.png", plot = p_azimuth, width = 8, height = 9, limitsize = FALSE)
p_azimuth

p_azimut_vs_seurat <- p_umap_int | p_azimuth
ggsave("plots/4_17_azimut_vs_seurat.png", plot = p_azimut_vs_seurat, width = 16, height = 9, limitsize = FALSE)
p_azimut_vs_seurat 


