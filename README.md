# Single-Cell Transcriptomic Atlas of Thyroid Carcinoma & Radioiodine (RAI) Refractoriness

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.12.4-brightgreen.svg)](https://scanpy.readthedocs.io/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.14.0-orange.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end single-cell RNA-sequencing (scRNA-seq) and machine learning pipeline mapping the tumor microenvironment (TME), follicular thyrocyte differentiation trajectories, ligand-receptor cell-cell communication networks, single-cell Graph Neural Networks (GCN/GAT), and in silico therapeutic redifferentiation in papillary and radioiodine-refractory thyroid carcinoma.

---

## 🔬 Project Architecture & Execution Pipeline

The project is structured into 8 modular, reproducible phases spanning primary data ingestion, quality control, single-cell graph modeling, ligand-receptor signaling, in silico perturbation, and TCGA clinical cohort validation:

```
thyroid_scRNA/
├── 01_Phase1_Ingestion.ipynb                          # Harmonization & raw atlas assembly (195,352 cells)
├── 02_Phase2_QC_BatchCorrection.ipynb                 # Doublet/filtering, log-norm, Harmony patient integration
├── 03_Phase3_Clustering_Annotation.ipynb              # Leiden clustering, marker-based lineage scoring (8 cell types)
├── 04_Phase4_Pseudotime_TherapyReversal.ipynb         # Diffusion pseudotime (DPT) & 10-gene RAI Index
├── 05_Phase5_SingleCell_GraphAI.ipynb                 # Single-cell Graph Convolutional Networks (GCN) on cell kNN
├── 06_Phase6_TCGA_Recurrence_Risk.ipynb               # Clinical recurrence risk & Cox PH survival modeling (TCGA-THCA)
├── 07_Ligand–Receptor Communication Analysis.ipynb    # TME-thyrocyte communication atlas & received signal scores
├── 08_In_Silico_Knockdown_Therapy_Reversal.ipynb      # In silico knockdown perturbation, Delta RAI & pseudotime shift
├── processed_data/                                    # Generated analytical tables, models, and publication figures
│   ├── ligand_receptor_communication_atlas.csv
│   ├── thyrocyte_received_signals_by_pathway.csv
│   ├── thyrocyte_per_cell_received_signals.csv
│   ├── thyrocyte_lr_pathway_heatmap.png
│   ├── rai_index_10_stage_summary.csv
│   ├── rai_index_trajectory_overlay.png
│   ├── in_silico_knockdown_results.csv
│   ├── in_silico_knockdown_summary.png
│   ├── genome_wide_rai_correlations.csv
│   ├── cox_summary.csv
│   ├── final_recurrence_risk.csv
│   ├── gcn_test_metrics.csv
│   └── gcn_predictions.csv
└── README.md
```

---

## 🧪 Key End-Term Deliverables & Verified Results

### 1. 10-Gene Thyroid Differentiation Score (RAI Index) & Trajectory Overlay
- **Standardized 10-Gene Panel**: `SLC5A5` (NIS), `TPO`, `TG`, `PAX8`, `NKX2-1`, `DIO1`, `DIO2`, `DUOX1`, `DUOX2`, `TSHR` (evaluated via cell-level standardized z-scoring).
- **Strong Trajectory Anti-Correlation**: Pearson correlation with diffusion pseudotime is **$r = -0.642$** ($p < 10^{-300}$), confirming progressive loss of iodine-metabolizing machinery along the transition from normal follicular epithelium to aggressive metastasis.
- **Stage Progression Breakdown**:
  - **Paratumor (Normal Follicular)**: Mean $+0.7025 \pm 0.4893$
  - **Primary Tumor**: Mean $-0.0955 \pm 0.4409$
  - **Lymph Node (LN) Metastasis**: Mean $-0.0680 \pm 0.4337$
  - **RAI-Refractory Metastasis**: Mean $-0.3176 \pm 0.1837$
- **Figure**: Generated in [`processed_data/rai_index_trajectory_overlay.png`](processed_data/rai_index_trajectory_overlay.png).

---

### 2. Ligand–Receptor Communication Analysis & Pathway "Received Signal" Scores
- **Atlas-Wide Intercellular Network**: Evaluated 1,024 sender–target pairs spanning 8 lineages (Thyrocytes, T-cells, NK cells, B-cells, Plasma cells, Myeloid cells, Fibroblasts, Endothelial cells) across canonical signaling cascades (Chemokine, TGF-beta, ErbB/EGFR, FGF, VEGF, Interleukin, TNF, AXL/RTK, Checkpoint, Integrin/ECM, Immune/Stromal).
- **Dominant Microenvironmental Signals Received by Thyrocytes**:
  1. **Stromal Extracellular Matrix Remodeling**: Fibroblast-derived `FN1` interacting with thyrocyte Integrin receptors `ITGAV` and `ITGB1` (communication scores: $1.294$ and $1.096$).
  2. **Myeloid / Macrophage Paracrine Signaling**: Myeloid-secreted `MIF` and `APP` signaling through thyrocyte `CD74` (communication score: $1.644$).
  3. **Endothelial Paracrine Crosstalk**: `VEGFA` signaling to VEGFR receptors `FLT1` and `KDR`.
- **Per-Cell Received Signal Scoring**: Successfully calculated per-cell received signal vectors for all 49,419 thyrocytes, saved to [`processed_data/thyrocyte_per_cell_received_signals.csv`](processed_data/thyrocyte_per_cell_received_signals.csv).
- **Figure**: Pathway communication matrix saved in [`processed_data/thyrocyte_lr_pathway_heatmap.png`](processed_data/thyrocyte_lr_pathway_heatmap.png).

---

### 3. In Silico Knockdown Perturbation & Pseudotime Redifferentiation Reversal
To identify therapeutic intervention points capable of restoring thyroid differentiation and reversing RAI refractoriness, we implemented an in silico single-cell gene regulatory network perturbation model:
- **Predicted Change in RAI Index ($\Delta$ RAI)**: Measures restoration of the 10-gene differentiation module.
- **Predicted Trajectory Shift**: Measures repositioning along diffusion pseudotime (negative shift indicates therapeutic reversal toward Paratumor state).
- **Perturbation Benchmark**:
  - **Top Driver Inactivations**:
    - **`S100A6` Knockdown**: $\Delta\text{RAI} = +1.2965$, $\Delta\text{Pseudotime} = -0.3246$ (Strongest Redifferentiation Reversal)
    - **`FN1` Knockdown**: $\Delta\text{RAI} = +0.4106$, $\Delta\text{Pseudotime} = -0.1028$ (Therapeutic Redifferentiation)
    - **`S100A4` Knockdown**: $\Delta\text{RAI} = +0.3616$, $\Delta\text{Pseudotime} = -0.0905$ (Therapeutic Redifferentiation)
    - **`LGALS3` Knockdown**: $\Delta\text{RAI} = +0.3524$, $\Delta\text{Pseudotime} = -0.0882$ (Therapeutic Redifferentiation)
    - **`FOS` / `JUN` (AP-1 Complex) Knockdown**: $\Delta\text{RAI} \approx +0.191$, $\Delta\text{Pseudotime} \approx -0.048$ (Therapeutic Redifferentiation)
  - **Candidate Drivers (`RUNX2` & `SOX9`)**: Minimal baseline expression in thyrocytes (6.3% and 16.3% positive cells), yielding minor/neutral shifts ($\Delta\text{RAI} = +0.0036$ and $-0.0135$), supporting their classification as **candidate drivers** rather than master repressors.
  - **Positive Differentiation Controls (`PAX8` & `NKX2-1`)**: Knockdown causes pronounced differentiation collapse ($\Delta\text{RAI} = -0.3226$ and $-0.2585$, $\Delta\text{Pseudotime} = +0.0808$ and $+0.0647$).
  - **Random Gene Controls**: Compared against `GAPDH`, `ACTB`, `B2M`, `EEF1A1`, `MALAT1` to contextualize variance.
- **Deliverables**: Detailed results in [`processed_data/in_silico_knockdown_results.csv`](processed_data/in_silico_knockdown_results.csv) and visual comparisons in [`processed_data/in_silico_knockdown_summary.png`](processed_data/in_silico_knockdown_summary.png).

---

### 4. TCGA-THCA Recurrence Risk Modeling & Limitations Analysis
- **Clinical Validation Cohort**: 502 primary thyroid tumors from the TCGA-THCA cohort evaluated for Progression-Free Interval (PFI).
- **Multivariable Cox Proportional Hazards Model**:
  - **Pathologic Stage**: Hazard Ratio $\text{HR} = 1.570$ ($95\%\text{ CI: } 1.161 - 2.125$, $p = 0.00336$) — statistically significant prognostic factor.
  - **Age at Diagnosis**: $\text{HR} = 1.000$ ($p = 0.971$).
  - **Single-Cell RAI Score**: $\text{HR} = 1.003$ ($p = 0.989$).
- **Limitations & Methodological Rationale**:
  1. **Cohort Discordance**: The single-cell RAI-resistance signature is derived from **distant metastatic, RAI-refractory lesions** (stage IV/metastasis), whereas the public TCGA-THCA cohort consists almost exclusively of **resectable primary thyroid carcinomas** (early-stage surgery specimens) with very low overall event rates (recurrence events $n \approx 42 / 502$, $<9\%$).
  2. **Inter-Patient Heterogeneity**: The single-cell RAI-refractory metastasis group is sourced from $n=2$ patients (`PTC4`, `PTC11`), which introduces patient-specific clonal features that necessitate validation in larger multi-center metastatic cohorts.
  3. **Therapeutic Confounding**: TCGA patients undergo heterogeneous post-operative adjuvant therapies (variable RAI dosing, thyroid hormone suppression, external beam radiation), decoupling baseline primary tumor expression from late-recurrence biology.

---

## 🚀 Reproduction & Setup

### Environment Requirements
- Python 3.10+
- Packages: `scanpy`, `anndata`, `torch`, `torch_geometric`, `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `lifelines`, `scikit-learn`

```bash
# Clone the repository
git clone https://github.com/SatyaKrishna2811/thyroid_scRNA.git
cd thyroid_scRNA

# Install dependencies
pip install scanpy anndata torch torchvision torch_geometric pandas numpy scipy matplotlib seaborn lifelines scikit-learn
```

### Running the End-Term Analyses
Execute the Jupyter notebooks in numerical sequence:
- Notebooks `01` through `06`: Ingestion, QC, Clustering, Trajectory, Graph AI, and TCGA validation.
- Notebook `07_Ligand–Receptor Communication Analysis.ipynb`: Full ligand-receptor network and received signal scoring.
- Notebook `08_In_Silico_Knockdown_Therapy_Reversal.ipynb`: In silico knockdown perturbation, Delta RAI, and trajectory reversal.
