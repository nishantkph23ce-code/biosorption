# Binary Bio-adsorbent System for Azo Dye Removal: Rice Husk + Spent Tea Leaves

**Authors:** Nishant Kumar et al.

**Repository for:** *Binary Bio-Adsorbent System for Azo Dye Removal Using Rice Husk and Spent Tea Leaves: Kinetics, Four-Isotherm Analysis, and ANN Dose Optimisation*

*(Under review / submitted for publication)*

[![DOI](https://img.shields.io/badge/DOI-pending-blue)](.)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

---

## Repository Structure

```
biosorption-RH-STL/
├── README.md
├── LICENSE
├── requirements.txt
├── data/
│   ├── 00_calibration.csv
│   ├── 01_RH_contact_time.csv
│   ├── 02_STL_contact_time.csv
│   ├── 03_RH_dose_variation.csv
│   ├── 04_STL_dose_variation.csv
│   ├── 05_binary_RH_constant.csv
│   └── 06_binary_STL_constant.csv
├── analysis/
│   └── run_all_analysis.py
└── results/
    └── parameters_summary.csv
```

---

## Experimental Conditions

| Parameter | Value |
|-----------|-------|
| Initial dye concentration (C₀) | 50 mg/L |
| Solution volume (V) | 200 mL |
| Temperature | 29 ± 1°C |
| Solution pH | 9.12 (natural) |
| Stirring speed | 150 rpm |
| Equilibrium time | 120 min |
| Particle size | 300 µm (both adsorbents) |

## Key Results

| Model | Parameter | RH | STL |
|-------|-----------|-----|-----|
| PSO kinetics | R² | 0.955 | 0.985 |
| PSO kinetics | qe,cal (mg/g) | 12.30 | 14.34 |
| Freundlich (best isotherm) | R² | 0.925 | 0.970 |
| D-R mean free energy | E (kJ/mol) | 8.47 | 6.35 |

- **Binary system:** Sf = 1.087 (8.7% synergistic enhancement)
- - **ANN optimum:** RH = 0.66 g + STL = 1.80 g per 200 mL → RE ≥ 95%
 
  - ## Quickstart
 
  - ```bash
    git clone https://github.com/nishantkph23ce-code/biosorption.git
    cd biosorption
    pip install -r requirements.txt
    python analysis/run_all_analysis.py
    ```

    ## Citation

    If you use this work, please cite:

    > Nishant Kumar et al. (2026). *Binary Bio-Adsorbent System for Azo Dye Removal Using Rice Husk and Spent Tea Leaves: Kinetics, Four-Isotherm Analysis, and ANN Dose Optimisation*. (Under review)
    
    ## License

    MIT License. See [LICENSE](LICENSE) for details.

    ## Contact

    Department of Civil Engineering, Katihar Engineering College, Katihar, Bihar, India
