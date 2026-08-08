# TNBC IIR Score — Integrated Immunotherapy Readiness Biomarker

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![R 4.0+](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

## 📌 Overview

Triple-Negative Breast Cancer (TNBC) is an aggressive subtype with limited treatment options. While immune checkpoint inhibitors (ICIs) have shown promise, only a subset of patients respond. Current biomarkers like PD-L1 and TMB are insufficient for accurate patient selection.

We developed the **Integrated Immunotherapy Readiness (IIR) Score** — a multi-omic biomarker that integrates:

| Component | Description |
|-----------|-------------|
| **Immune infiltration** | xCell deconvolution of immune cell types |
| **Checkpoint activation** | PD1/PD-L1 pathway signature |
| **Spatial architecture** | Immune-Stroma ratio, Exclusion index |
| **Mutational context** | TMB, APOBEC signatures |
| **Drug penetration** | Drug Penetration Index (DPI) |
| **Immune escape** | MHC-I expression, HLA status |

---

## 🎯 Key Findings

| Finding | Result |
|---------|--------|
| **IIR Score validation** | HR = 0.689, p = 0.0028 (METABRIC, n=320) |
| **Survival difference** | 50-month median OS (High vs Poor IIR) |
| **IIR vs PD-L1** | Comparable performance (C-index: 0.582 vs 0.577) |
| **IIR vs TMB** | IIR significantly outperforms TMB |
| **Spatial validation** | r = 0.719 with CD8+ signature |
| **Clinical utility** | Decision curve analysis confirms net benefit |

---

## 🧬 Repository Structure

TNBC-IIR-Score-Immunotherapy/
│
├── Phase_I-XI/ # Original research pipeline
│ ├── Phase_I_Immune_Deconv/ # xCell deconvolution, ImmuneScore
│ ├── Phase_II_Immune_Subtypes/ # BLIS/IM/BLIA subtypes
│ ├── Phase_III_Immunotherapy/ # PD1/PD-L1 signature
│ ├── Phase_IV_Spatial_Metrics/ # Spatial immune-stroma metrics
│ ├── Phase_V_Drug_Penetration/ # Drug Penetration Index (DPI)
│ ├── Phase_VI_Mutations/ # TMB and mutational landscape
│ ├── Phase_VII_Mutational_Signatures/ # APOBEC, Ageing signatures
│ ├── Phase_VIII_IIR_Score/ # Integrated IIR Score
│ ├── Phase_IX_Clinical_Decision/ # Clinical decision engine
│ ├── Phase_X_HLA_ImmuneEscape/ # HLA-I immune escape
│ └── Phase_XI_Escape_Adjusted_ICB/ # Escape-adjusted ICB score
│
├── GAP_FIXING/ # Gap-filling validation studies
│ ├── GAP_01_External_Validation/ # METABRIC validation
│ ├── GAP_02_ICB_Treated_Cohort/ # ICB cohort analysis
│ ├── GAP_03_DPI_Formula/ # DPI formula definition
│ ├── GAP_04_HLA_LOH/ # HLA LOH improvement
│ ├── GAP_05_TMB_Standardization/ # TMB standardization
│ ├── GAP_06_Biomarker_Comparison/ # IIR vs PD-L1 vs TMB
│ ├── GAP_07_Decision_Curves/ # Clinical decision curves
│ ├── GAP_08_Statistical_Power/ # Power analysis
│ ├── GAP_09_Neoantigen_Prediction/ # Neoantigen prediction
│ ├── GAP_10_MultiOmics/ # Multi-omic integration
│ ├── GAP_11_Spatial_Validation/ # Spatial metrics validation
│ └── GAP_12_Drug_Sensitivity/ # Drug sensitivity correlation
│
├── GAP_NOTES/ # Complete scientific reports
│ ├── GAP_01_notes.md
│ ├── GAP_02_notes.md
│ └── ... (all 12 gaps)
│
├── results/ # Key result files
│ ├── METABRIC_validation/
│ ├── TMB_standardization/
│ ├── Biomarker_comparison/
│ └── ...
│
├── README.md # This file
├── LICENSE # MIT License
└── .gitignore # Git ignore rules



---

## 📊 Clinical Implications

| Patient Group | Recommended Strategy |
|---------------|---------------------|
| **High IIR + No escape** | ICB monotherapy |
| **High IIR + Escape** | ICB + priming (IFN-γ, epigenetic, STING) |
| **Poor IIR + No escape** | Chemo + TME priming |
| **Poor IIR + Escape** | Non-ICB + aggressive TME remodeling |

---

## 🛠️ Requirements

### Python

```bash
pip install pandas numpy scipy lifelines matplotlib seaborn scikit-learn statsmodels
R
r
install.packages(c("data.table", "survival", "survminer", "ggplot2", "xCell"))


📈 Key Results

METABRIC Validation
Metric	Value
Samples	320
Events	168
IIR HR	0.689 (95% CI: 0.540-0.879)
p-value	0.0028
C-index	0.582
Biomarker Comparison
Biomarker	C-index	AIC
IIR	0.582	1741.6
PD-L1	0.577	1741.8
ImmuneScore	0.589	1741.2
TMB	0.487	1750.5
Spatial Metrics Validation
Metric	Correlation	p-value
ImmuneScore vs CD8+	r = 0.719	< 0.001
Immune-Stroma Ratio vs CD8+	r = 0.622	< 0.001
Exclusion Index vs CD8+	r = -0.622	< 0.001

📝 Publication Status
Manuscript in preparation

👤 Author
Bhaskararao Ch (Baashi)
Independent Researcher
GitHub: Baashi27-ai

📄 License
MIT License — Free for academic and clinical use with attribution.

🙏 Acknowledgments
TCGA Research Network for providing data

METABRIC Consortium for validation cohort

DepMap/PRISM for drug sensitivity data

📧 Contact
For questions, collaborations, or data access requests, please open an issue on GitHub.

