# thyroid_scRNA

Why does papillary thyroid cancer stop responding to radioiodine?

This repo is a small end-to-end attempt at that question using single-cell RNA-seq. It is six Jupyter notebooks, one per phase, each one reading the previous phase's output and writing its own. It starts with ~198k raw cells from 23 patient samples. It ends with a 50-gene "RAI resistance" signature tested against 502 real TCGA patients. There is no framework and no config system, just scanpy, a two-layer GCN and a Cox model.

The short version of the result: the single-cell part works nicely. The clinical validation does not (yet). Both are written up below.

## the pipeline

```
raw 10x matrices (23 samples, 11 patients)
  -> 01 ingest          one AnnData, tagged with patient + tissue stage
  -> 02 QC + Harmony    filter, normalize, 3k HVGs, PCA, batch-correct by patient
  -> 03 cluster         kNN -> UMAP -> Leiden -> marker-score annotation -> pull out thyrocytes
  -> 04 trajectory      diffusion pseudotime, differentiation score, DEG -> 50-gene signature
  -> 05 GCN             predict tumour stage per cell from the cell graph
  -> 06 TCGA            score the signature on bulk RNA-seq, Cox + Kaplan-Meier on PFI
```

The four stages, in the order we expect the disease to move:

| stage | sample suffix | thyrocytes |
|---|---|---|
| Paratumor | `_P` | 6,823 |
| Primary Tumor | `_T` | 20,533 |
| LN Metastasis | `_LeftLN`, `_RightLN` | 16,723 |
| RAI-Refractory Metastasis | `_SC` (PTC4, PTC11) | 5,340 |

## what each notebook does

**01_Phase1_Ingestion.** Reads the 23 `matrix/features/barcodes` triplets in `data/` and concatenates them into `adata_raw_master.h5ad` (197,954 cells x 33,694 genes), with `patient_id` and `tissue_type` in `.obs`. Heads up: this notebook is currently empty in the repo. Grab `adata_raw_master.h5ad` from the release (see below) and start at Phase 2.

**02_Phase2_QC_BatchCorrection.** Keeps cells with 200-7,500 genes and <15% mitochondrial reads, and genes seen in at least 3 cells. That leaves 195,352 cells. Then it normalizes to 10k counts/cell, applies log1p, picks 3,000 HVGs (seurat flavor), scales, runs 50-component PCA and runs Harmony on `patient_id`. Harmony converges in 3 iterations. Raw counts are kept in `adata.raw`.

**03_Phase3_Clustering_Annotation.** Builds a 15-NN graph on the Harmony embedding, then UMAP, then Leiden (res 0.8, 24 clusters). Each cluster gets whichever of 9 lineage marker panels scores highest (thyrocytes, T, NK, B, plasma, myeloid, mast, fibroblasts, endothelial). The panel is deliberately broad so immune clusters don't get forced into "thyrocyte". This gives 49,419 thyrocytes, and those carry forward.

**04_Phase4_Pseudotime_TherapyReversal.** Rebuilds the graph on thyrocytes only and computes a diffusion map. The root is the Paratumor cell at the far end of DC1, on the side away from the RAI-refractory cells. Then it computes DPT pseudotime:

| | Paratumor | Primary | LN met | RAI-refractory |
|---|---|---|---|---|
| mean pseudotime | 0.258 | 0.732 | 0.728 | 0.744 |
| thyroid differentiation score | 1.107 | 0.207 | 0.135 | -0.084 |

The differentiation score (TDS) is 15 iodide-handling genes: `TG`, `TPO`, `SLC5A5` (NIS), `TSHR`, `PAX8`, and so on. This is the whole story in one row. Radioiodine only works if the cell can take up iodide, and these cells gradually forget how. Then a Wilcoxon test of RAI-refractory vs Paratumor: 1,992 genes up, 1,562 down (|logFC| > 1, FDR < 0.05). The down list is basically the thyroid (`TFF3`, `TG`, `TPO`, `SLC26A7`, `GPX3`, metallothioneins). The up list is S100 family, `LGALS3`, `FN1`, `TIMP1`, `KRT19`, `CLDN4`, which reads like dedifferentiation plus an EMT flavor. It drops ribosomal, mito, sex-chromosome and unannotated IDs (130 genes) and keeps the top 50 as `rai_signature_genes.csv`.

**05_Phase5_SingleCell_GraphAI.** Every thyrocyte is a node, the Phase 4 kNN connectivities are edges (1.16M), and the 50 Harmony PCs are features. A 2-layer GCN (`GCNConv` 50 -> 64 -> 4, dropout 0.2) predicts stage. It uses a stratified 70/15/15 split, class-weighted NLL, 200 epochs of Adam, and keeps the checkpoint with the best val balanced accuracy. A logistic regression on the same features, without the graph, is the baseline:

| model | accuracy | balanced accuracy |
|---|---|---|
| GCN | 0.761 | **0.803** |
| logistic regression (no graph) | 0.612 | 0.665 |

The graph is worth ~14 points of balanced accuracy. RAI-refractory recall is 0.815 but precision is 0.521, so the model is eager to call cells refractory.

**06_Phase6_TCGA_Recurrence_Risk.** Downloads TCGA-THCA from the UCSC Xena hub (HiSeqV2 expression, survival, clinical matrix). The notebook skips the download if the files are already in `data/tcga/`, which they are. It keeps only primary tumours (505), z-scores the signature genes (48/50 found) and averages them into an `RAI_Score`. It then joins that with progression-free interval, age and stage (502 patients, 52 events) and fits a multivariable Cox model:

| covariate | HR | 95% CI | p |
|---|---|---|---|
| RAI_Score | 1.003 | 0.65-1.55 | 0.989 |
| age | 1.000 | 0.98-1.02 | 0.971 |
| stage | 1.570 | 1.16-2.12 | 0.003 |

Concordance is 0.617. A median split on RAI_Score gives 32 vs 20 events (high vs low), with log-rank p ≈ 0.12.

## the honest part

The signature does not predict recurrence in TCGA once you adjust for stage. Stage does all the work. There are a few reasons this is not too surprising, and they are the most useful thing in this README:

- **n = 2 patients.** The RAI-refractory group is PTC4 and PTC11. The DEG p-values in Phase 4 are computed per cell, so they look like 0.0, but the real sample size is 2. Treat the signature as candidates, not a result.
- **Wrong target population.** TCGA-THCA is mostly early-stage, mostly curable primary tumours with few events (52/502). A signature learned from refractory *metastases* is being asked to stratify patients who mostly never get there.
- **Bulk vs single cell.** A thyrocyte program is being scored in bulk tissue, where immune and stromal content dilutes it. Some of the up genes (S100A4, FN1, TIMP1) are also expressed by fibroblasts and macrophages, so in bulk they partly measure tumour composition.
- **The GCN number is transductive.** Test cells share kNN edges with training cells from the same specimen, so 0.803 is a within-sample estimate. A patient-held-out split would be the honest version, but you can't hold out a 2-patient class cleanly.

`cox_summary_unfiltered_signature.csv` is the same Cox fit before the technical-gene filter (HR 1.017, p = 0.94). Filtering was not the problem.

## running it

```bash
git clone https://github.com/SatyaKrishna2811/thyroid_scRNA.git
cd thyroid_scRNA
pip install scanpy anndata harmonypy igraph leidenalg torch torch_geometric scikit-learn lifelines pandas numpy matplotlib jupyter
```

The intermediate `.h5ad` files are 0.5-1.4 GB each, so they are not in git. They are attached to the [v1.0.0 release](https://github.com/SatyaKrishna2811/thyroid_scRNA/releases/tag/v1.0.0). Put them in `processed_data/`:

| file | produced by | size |
|---|---|---|
| `adata_raw_master.h5ad` | Phase 1 | 0.64 GB |
| `adata_processed_qc.h5ad` | Phase 2 | 1.35 GB |
| `adata_clustered.h5ad` | Phase 3 | 1.41 GB |
| `adata_thyrocytes.h5ad` | Phase 3 | 0.52 GB |
| `adata_trajectory.h5ad` | Phase 4 | 0.52 GB |

Then run notebooks 02 -> 06 in order from the repo root (paths are `os.getcwd()`-relative). Or download just the one you need and jump straight to the phase you care about. Phase 5 and 6 only need `adata_trajectory.h5ad` and the CSVs already in the repo.

Rough cost on a laptop CPU: Harmony is ~10 min, the full-atlas UMAP ~4 min, the Phase 4 Wilcoxon ~3 min, and the GCN a few minutes. Memory is the real constraint with ~195k cells. The notebook casts `X` to float32 for that reason.

## layout

```
data/
  GSM55851xx_PTC*_{barcodes,features,matrix}.*.gz   raw 10x counts, 23 samples (GEO GSE184362)
  tcga/                                             TCGA-THCA from UCSC Xena
processed_data/
  deg_RAI_vs_Paratumor.csv        full DEG table, Phase 4
  rai_signature_genes.csv         the 50 genes
  cell_gcn.pt                     GCN weights
  gcn_predictions.csv             per-cell prediction + split
  gcn_test_metrics.csv            GCN vs baseline
  final_recurrence_risk.csv       per-patient RAI_Score, PFI, age, stage
  cox_summary*.csv                Cox fits
0[1-6]_Phase*.ipynb               the pipeline
```

## data

- Single-cell: GEO [GSE184362](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184362), papillary thyroid carcinoma, 11 patients, paratumor / primary / lymph node / subcutaneous metastasis (Pu et al., *Nature Communications*, 2021).
- Bulk: TCGA-THCA via the [UCSC Xena](https://xenabrowser.net/) public hub. PFI is the recommended endpoint for THCA, because thyroid cancer rarely kills people, so OS has almost no events.

## todo

- actually commit the Phase 1 ingestion code
- `requirements.txt` with pinned versions (it was run on scanpy 1.12.4)
- pseudobulk DEG per patient instead of per cell, so the p-values mean something
- score the signature on thyrocyte-deconvolved TCGA expression, or on a cohort that actually has RAI-refractory patients
- patient-held-out GCN eval once there are more than 2 refractory patients
- cell-cell communication between dedifferentiated thyrocytes and the fibroblast/myeloid compartments, which are expanding in tumour and LN

## license

No license file yet. Until one is added, all rights are reserved by default.
