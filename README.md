# Project-Bioinformatika-2026
# 🧬 Analisis Differential Expression Genes (DEG) pada Kanker Pankreas
### Menggunakan Dataset GEO Publik (GSE10072) dengan Pendekatan Limma

> **Mini Proyek Bioinformatika** | Bioinformatika | 

---

## 👥 Anggota Kelompok

| Nama | NIM |
|------|-----|
|Ade Riyana P.|499349
|Hanin Izdihar|517672
|Adha Nurkholifah|518426
|Talitha Daris N.|518540
|Devi Korniasari|519029

---

## 📌 Deskripsi Proyek

Proyek ini melakukan re-analisis dataset ekspresi gen publik dari NCBI GEO untuk menganalisis Ekspresi Gen Diferensial dan Enrichment Pathway pada Kanker Pankreas.

Selain analisis DEG standar menggunakan **limma**, proyek ini menambahkan:
- **GO/KEGG Enrichment Analysis** untuk mengetahui jalur biologis yang terlibat
- **Visualisasi komprehensif** (Volcano Plot, Heatmap, UMAP, Dotplot enrichment)

---

## 🗂️ Dataset

| Atribut | Detail |
|---------|--------|
| **GEO Accession** | [GEO GSE28735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) |
| **Platform** | Affymetrix Human Gene 1.0 ST Array  |
| **Sampel** | xxxx |
| **Organisme** | *Homo sapiens* |


---

## 🛠️ Tools & Library yang Digunakan

```r
# Bioconductor
GEOquery       # Download data dari NCBI GEO
limma          # Analisis DEG (Linear Models for Microarray Data)
hgu133a.db     # Anotasi gen untuk platform Affymetrix GPL96
AnnotationDbi  # Interface anotasi gen
clusterProfiler # GO & KEGG Enrichment Analysis  ← TAMBAHAN
org.Hs.eg.db   # Database gen manusia            ← TAMBAHAN

# CRAN
ggplot2        # Volcano plot, visualisasi
pheatmap       # Heatmap ekspresi gen
dplyr          # Manipulasi data
umap           # Dimensionality reduction
```

---

## 🔬 Pipeline Analisis

```
┌─────────────────────────────────────────────────┐
│              PIPELINE ANALISIS DEG              │
├─────────────────────────────────────────────────┤
│                                                 │
│  1. Download data GEO (getGEO)                  │
│         ↓                                       │
│  2. Pre-processing                              │
│     - Cek & log2 transform                      │
│     - Ekstrak matriks ekspresi                  │
│         ↓                                       │
│  3. Definisi grup (Tumor vs Normal)             │
│         ↓                                       │
│  4. Limma DEG Analysis                          │
│     - Design matrix                             │
│     - lmFit → eBayes → topTable                 │
│         ↓                                       │
│  5. Anotasi Gen (hgu133a.db)                    │
│         ↓                                       │
│  6. Visualisasi                                 │
│     - Boxplot distribusi                        │
│     - UMAP plot                                 │
│     - Volcano plot                              │
│     - Heatmap top 50 DEG                        │
│         ↓                                       │
│  7. GO/KEGG Enrichment [TAMBAHAN]               │
│     - Biological Process (BP)                   │
│     - KEGG Pathway                              │
│         ↓                                       │
│  8. Simpan hasil (CSV)                          │
│                                                 │
└─────────────────────────────────────────────────┘
```

---

## 📁 Struktur Repository

```
NomorKelompok_NamaProyek/
├── README.md
├── scripts/
│   └── analysis.R          # Script R lengkap
├── results/
│   ├── boxplot.png
│   ├── umap_plot.png
│   ├── volcano_plot.png
│   ├── heatmap_top50.png
│   ├── GO_enrichment.png   ← TAMBAHAN
│   ├── KEGG_pathway.png    ← TAMBAHAN
│   └── DEG_results.csv
└── paper/
    └── makalah.pdf
```

---

## ▶️ Cara Menjalankan

1. Install R (≥ 4.0) dan RStudio
2. Buka `scripts/analysis.R` di RStudio
3. Jalankan seluruh script (Ctrl+A → Ctrl+Enter)
4. Semua output akan tersimpan di folder `results/`

> **Catatan:** Script akan otomatis menginstall package yang belum ada.

---

## 📊 Hasil Utama

- Total DEG teridentifikasi: **XXX gen** (padj < 0.05, |logFC| > 1)
- Gen upregulated: **XXX**
- Gen downregulated: **XXX**
- Top hub genes: 
- Jalur KEGG paling signifikan: 

---

## 📄 Luaran

| Luaran | 
|--------|
| Script R|
| Makalah (PDF) |
| Poster (A0) |


