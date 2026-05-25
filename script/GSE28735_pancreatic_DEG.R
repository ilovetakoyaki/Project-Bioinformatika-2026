# =============================================================
# ANALISIS DEG: KANKER PANKREAS (PDAC)
# Dataset  : GSE28735 (45 tumor vs 45 normal, paired samples)
# Platform : GPL6244 - Affymetrix Human Gene 1.0 ST Array
# Tool DEG : limma
# Tambahan : GO & KEGG Enrichment Analysis
# =============================================================

# -------------------------------------------------------
# PART A. INSTALL & LOAD PACKAGES
# -------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install semua package yang dibutuhkan
BiocManager::install(
  c("GEOquery", "limma",
    "hugene10sttranscriptcluster.db",  # annotation GPL6244
    "AnnotationDbi",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"),
  ask = FALSE, update = FALSE
)

install.packages(c("pheatmap", "ggplot2", "dplyr", "umap"),
                 repos = "https://cran.rstudio.com/")

# Load semua library
library(GEOquery)
library(limma)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(umap)

# Buat folder output
if (!dir.exists("results")) dir.create("results")

# -------------------------------------------------------
# PART B. DOWNLOAD DATA DARI GEO
# -------------------------------------------------------

cat("Downloading GSE28735...\n")
gset <- getGEO("GSE28735", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Cek struktur data
cat("Dimensi data:", dim(exprs(gset)), "\n")
cat("Jumlah sampel:", ncol(exprs(gset)), "\n")

# -------------------------------------------------------
# PART C. PRE-PROCESSING DATA EKSPRESI
# -------------------------------------------------------

ex <- exprs(gset)

# Cek apakah perlu log2 transform
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
  cat("Log2 transform diterapkan.\n")
} else {
  cat("Data sudah dalam skala log, tidak perlu transform.\n")
}

# -------------------------------------------------------
# PART D. DEFINISI KELOMPOK SAMPEL
# -------------------------------------------------------

# Cek kolom metadata yang tersedia
cat("\nKolom metadata tersedia:\n")
print(colnames(pData(gset)))

# Cek isi kolom untuk identifikasi grup
cat("\nNilai unik di source_name_ch1:\n")
print(unique(pData(gset)[["source_name_ch1"]]))

# Definisi grup berdasarkan metadata
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)

nama_grup <- levels(gset$group)
cat("\nGrup yang teridentifikasi:\n")
print(nama_grup)

# -------------------------------------------------------
# PART E. VISUALISASI AWAL
# -------------------------------------------------------

# --- E.1 Boxplot distribusi ekspresi ---
group_colors <- as.numeric(gset$group)

png("results/boxplot_distribusi.png", width = 1400, height = 600)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel (GSE28735)",
  ylab = "Expression Value (log2)",
  cex.axis = 0.5
)
legend("topright",
       legend = levels(gset$group),
       fill = unique(group_colors),
       cex = 0.8)
dev.off()
cat("Saved: results/boxplot_distribusi.png\n")

# --- E.2 Density Plot ---
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

p_density <- ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (Kanker Pankreas GSE28735)",
    x = "Expression Value (log2)",
    y = "Density"
  )

ggsave("results/density_plot.png", p_density, width = 8, height = 5, dpi = 300)
cat("Saved: results/density_plot.png\n")

# --- E.3 UMAP ---
umap_input <- t(ex)
umap_input <- umap_input[, apply(umap_input, 2, function(x) !any(is.na(x)))]

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: Tumor vs Normal (GSE28735 Pancreatic Cancer)",
    x = "UMAP 1",
    y = "UMAP 2"
  )

ggsave("results/umap_plot.png", p_umap, width = 8, height = 6, dpi = 300)
cat("Saved: results/umap_plot.png\n")

# -------------------------------------------------------
# PART F. DESIGN MATRIX & LIMMA DEG ANALYSIS
# -------------------------------------------------------

design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# Tentukan perbandingan: tumor vs normal
# Cek nama grup yang muncul dari print(nama_grup) di atas
# Sesuaikan nama di bawah dengan yang muncul
grup_tumor  <- nama_grup[1]   # biasanya tumor
grup_normal <- nama_grup[2]   # biasanya normal

contrast_formula <- paste(grup_tumor, "-", grup_normal)
cat("\nKontras:", contrast_formula, "\n")

# Limma pipeline
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Ambil semua hasil DEG
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)

cat("Total DEG signifikan (padj < 0.05):", nrow(topTableResults), "\n")

# -------------------------------------------------------
# PART G. ANOTASI NAMA GEN
# -------------------------------------------------------

probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

cat("Anotasi selesai. Contoh hasil:\n")
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

# Simpan hasil DEG
write.csv(topTableResults, "results/Hasil_GSE28735_DEG.csv", row.names = FALSE)
cat("Saved: results/Hasil_GSE28735_DEG.csv\n")

# -------------------------------------------------------
# PART H. VISUALISASI DEG
# -------------------------------------------------------

# --- H.1 Volcano Plot ---
volcano_data <- data.frame(
  logFC     = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene      = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1  & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

p_volcano <- ggplot(volcano_data,
                    aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot DEG Kanker Pankreas (GSE28735)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  )

ggsave("results/volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)
cat("Saved: results/volcano_plot.png\n")

# --- H.2 Heatmap Top 50 DEG ---
topTableResults_sorted <- topTableResults[order(topTableResults$adj.P.Val), ]
top50 <- head(topTableResults_sorted, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,
  top50$SYMBOL
)
rownames(mat_heatmap) <- gene_label

# Bersihkan data
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

png("results/heatmap_top50.png", width = 1200, height = 1000)
pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes - Pancreatic Cancer (GSE28735)"
)
dev.off()
cat("Saved: results/heatmap_top50.png\n")

# -------------------------------------------------------
# PART I. GO ENRICHMENT ANALYSIS
# -------------------------------------------------------

# Filter DEG signifikan
deg_filtered <- topTableResults %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

cat("\nJumlah DEG untuk enrichment:", nrow(deg_filtered), "\n")
cat("Upregulated  :", nrow(filter(deg_filtered, logFC > 1)), "\n")
cat("Downregulated:", nrow(filter(deg_filtered, logFC < -1)), "\n")

gene_list <- unique(deg_filtered$SYMBOL)

# Konversi ke Entrez ID
gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
entrez_ids <- gene_entrez$ENTREZID
cat("Entrez ID berhasil dikonversi:", length(entrez_ids), "\n")

# GO Enrichment
ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cat("GO BP terms signifikan:", nrow(as.data.frame(ego)), "\n")

# Visualisasi GO
p_go_dot <- dotplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment\nPancreatic Cancer (GSE28735)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_go_bar <- barplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment\nPancreatic Cancer (GSE28735)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/GO_dotplot.png", p_go_dot, width = 10, height = 7, dpi = 300)
ggsave("results/GO_barplot.png", p_go_bar, width = 10, height = 7, dpi = 300)
write.csv(as.data.frame(ego), "results/GO_BP_enrichment.csv", row.names = FALSE)
cat("Saved: GO plots & CSV\n")

# -------------------------------------------------------
# PART J. KEGG PATHWAY ENRICHMENT
# -------------------------------------------------------

ekegg <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

cat("KEGG pathway signifikan:", nrow(as.data.frame(ekegg)), "\n")

ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p_kegg_dot <- dotplot(ekegg_readable, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment\nPancreatic Cancer (GSE28735)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_kegg_bar <- barplot(ekegg_readable, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment\nPancreatic Cancer (GSE28735)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/KEGG_dotplot.png", p_kegg_dot, width = 10, height = 7, dpi = 300)
ggsave("results/KEGG_barplot.png", p_kegg_bar, width = 10, height = 7, dpi = 300)
write.csv(as.data.frame(ekegg_readable), "results/KEGG_enrichment.csv", row.names = FALSE)
cat("Saved: KEGG plots & CSV\n")

# -------------------------------------------------------
# PART K. RINGKASAN AKHIR
# -------------------------------------------------------

cat("\n========================================\n")
cat("ANALISIS SELESAI!\n")
cat("========================================\n")
cat("Dataset       : GSE28735 (Pancreatic Cancer)\n")
cat("Platform      : GPL6244 Affymetrix HuGene 1.0 ST\n")
cat("Total DEG     :", nrow(topTableResults), "\n")
cat("DEG (|FC|>1)  :", nrow(deg_filtered), "\n")
cat("GO BP terms   :", nrow(as.data.frame(ego)), "\n")
cat("KEGG pathways :", nrow(as.data.frame(ekegg)), "\n")
cat("\nFile tersimpan di folder results/:\n")
print(list.files("results/"))
