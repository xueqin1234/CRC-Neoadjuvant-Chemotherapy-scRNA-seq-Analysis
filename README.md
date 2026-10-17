# Single-Cell RNA Sequencing Analysis of Neoadjuvant Chemotherapy in Colorectal Cancer

This repository contains the bioinformatics analysis code for single-cell RNA sequencing (scRNA-seq) data from a study on the response to neoadjuvant chemotherapy in colorectal cancer (CRC). The pipeline is built on the Seurat framework and covers the complete analytical workflow — from data quality control and cell annotation to downstream multi-omics analyses — and is used to reproduce the figures in the corresponding manuscript.

## Repository Structure



```
├── fig1/                  # Data preprocessing (QC, DoubletFinder) and cell-type annotation
├── fig2/                  # Epithelial cell subclustering (NMF-based malignant cell classification, cell proportion enrichment analysis)
├── fig3/                  # Myeloid cell subclustering
├── fig4/                  # Stromal cell subclustering (e.g., FAP+ fibroblasts) and Monocle trajectory analysis
├── T&B/                   # T cell and B cell subclustering
├── infercnv/              # inferCNV copy-number variation analysis (distinguishing malignant vs. non-malignant epithelial cells)
├── cellphoneDB_example/     # CellPhoneDB cell–cell communication analysis
├── nichenert/             # NicheNet ligand–receptor interaction analysis (neutrophils → FAP+ fibroblasts)
└── README.md
```

## Analysis Modules



* **Data preprocessing and cell annotation**: 10x data loading, quality control (mitochondrial gene percentage, gene counts), doublet removal with DoubletFinder, batch integration with Harmony, and cell-type annotation.

* **Cell subclustering**: Identification and characterization of epithelial, myeloid, stromal, T and B cell subpopulations, including cell proportion and enrichment analyses.

* **Malignant cell identification**: inferCNV-based copy-number analysis to distinguish malignant from non-malignant epithelial cells, followed by NMF-based classification of malignant cells.

* **Cell–cell communication**: CellPhoneDB statistical analysis to identify ligand–receptor interactions, and NicheNet to predict key signaling pathways.

* **Trajectory inference**: Monocle pseudotime analysis of key subpopulations such as FAP+ fibroblasts.

Samples are stratified by treatment response (Sensitive / Resistant / Intermediate) and collected at multiple treatment time points (CRC\_0 / CRC\_2 / CRC\_4), enabling dissection of tumor microenvironment differences among response groups.

## Dependencies



* R (>= 4.0) with key packages: Seurat, harmony, sctransform, DoubletFinder, ComplexHeatmap, ggplot2, NMF, monocle, nichenetr, etc.

* Python / command-line tools: CellPhoneDB (`cellphonedb`), and the inferCNV matrix merging script (`03_merge_infercnv_matrix.py`).

## Usage



1. Each script contains a **path configuration section** at the top (e.g., `root_dir`, `data_dir`); please adjust the paths before running.

2. Raw data (10x count matrices, Seurat RDS objects, etc.) are not included in this repository; prepare them according to the comments in each script.

3. We recommend running the scripts in the following order: `fig1 → fig2/fig3/fig4/T&B → infercnv → cellphone_example / nichenert`.

## License

This code is provided for academic research use only. For licensing inquiries, please contact the corresponding author.