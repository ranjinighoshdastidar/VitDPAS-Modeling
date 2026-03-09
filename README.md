# VitDPAS-Modeling

Pharmacokinetic modeling of serum 25(OH)D3 kinetics, gene regulatory
network (GRN) construction and kernel reduction, ODE-based GRN modeling,
and parametric analysis for the VitDPAS longitudinal vitamin D study.

Part of a multi-repository analytical pipeline:
- [VitDPAS-DGE_RNAseq](https://github.com/ranjinighoshdastidar/VitDPAS-DGE_RNAseq)
- [VitDPAS-Classification](https://github.com/ranjinighoshdastidar/VitDPAS-Classification)
- **VitDPAS-Modeling** ← you are here

## Repository Contents

| File | Step | Language | Description |
|------|------|----------|-------------|
| `VitD_Modeling_script.R` | 4a | R 4.4.2 | Pharmacokinetic modeling (decay + supplementation ODE) |
| `VitaminD_modeling.ipynb` | 4b | Julia 1.11.4 | Numerical verification of PK model |
| `Dorothea_Network_and_Kernel_Reduction.R` | 5 | R 4.4.2 | GRN construction and kernel reduction |
| `GRN_Modeling_Final.ipynb` | 6 | Julia 1.11.4 | ODE-based GRN modeling and parameter optimisation |
| `Parametric_Analysis.R` | 7 | R 4.4.2 | Parametric analysis of optimised parameters |

## Pipeline Overview

Scripts must be run in the following order. Each step depends on
outputs from the previous step.
VitD_Modeling_script.R
↓  produces: ki and di (fixed VitD parameters per replicate)
VitaminD_modeling.ipynb
↓  numerical verification (recommended before proceeding)
Dorothea_Network_and_Kernel_Reduction.R
↓  produces: Model_Equations.txt, interaction matrices
GRN_Modeling_Final.ipynb
↓  produces: optimized_equations_combined.txt
optimized_equations_High1.txt
optimized_equations_Mid1.txt
optimized_equations_Low1.txt
Parametric_Analysis.R
↓  produces: Parametric_Analysis_Results_recent.xlsx

## Input Files Required

Not included in this repository (see Data Availability).

**VitD_Modeling_script.R**
- `VitDPAS_vitD_data.xlsx`   Serum 25(OH)D3 concentrations per participant
  at all 7 time points (d0, d1, d28, d29, d56, d57, d84)

Pharmacokinetic parameters are hardcoded from:
Rybiński M, Ghosh Dastidar R, Zawrotna N et al.
Sci Rep 16, 2997 (2026). https://doi.org/10.1038/s41598-025-32831-z
- DECAY_START:      53.952, 45.694, 55.112 ng/mL
- DECAY_END:        38.103, 41.865, 42.774 ng/mL
- DECAY_DAYS:       26.0, 30.0, 81.0 days
- INCREASE R1:      36.3, 42.2, 54.0 ng/mL at t = 0, 24, 48 h
- INCREASE R2:      38.1, 45.9, 45.7 ng/mL at t = 0, 24, 48 h
- INCREASE R3:      41.9, 48.2, 55.1 ng/mL at t = 0, 24, 48 h

**VitaminD_modeling.ipynb**
- No external files required. All parameters are hardcoded.

**Dorothea_Network_and_Kernel_Reduction.R**
- `normalized_combined.xlsx`   Output of File_Prep_and_Visualisation.R
  (from VitDPAS-DGE_RNAseq repository)
- DoRothEA regulon database accessed via the dorothea R package
- Seed genes (hardcoded): JUNB, ETS2, FOS, NR4A2, KLF10, CSRNP1, HIF1A

**GRN_Modeling_Final.ipynb**
- `Model_Equations.txt`        Hill equations (output of Dorothea script)
- `normalized_combined.xlsx`   CPM data for model fitting
- ki and di parameters from VitD_Modeling_script.R (hardcoded per replicate)

**Parametric_Analysis.R**
- `optimized_equations_combined.txt`   Optimised parameters: full cohort
- `optimized_equations_High1.txt`      Optimised parameters: stable high
- `optimized_equations_Mid1.txt`       Optimised parameters: stable mid
- `optimized_equations_Low1.txt`       Optimised parameters: stable low
  (all four files are outputs of GRN_Modeling_Final.ipynb)

## Output Files

| Script | Key Outputs |
|--------|-------------|
| `VitD_Modeling_script.R` | Decay curve plots, supplementation kinetics plots, combined figure, parameter bar charts |
| `VitaminD_modeling.ipynb` | Verification plots, parameter comparison figures |
| `Dorothea_Network_and_Kernel_Reduction.R` | GRN network figures (4 layouts, PNG + SVG), interaction matrices (.xlsx), Hill equation files |
| `GRN_Modeling_Final.ipynb` | optimized_equations_*.txt files, fitted trajectory plots |
| `Parametric_Analysis.R` | Parametric_Analysis_Results_recent.xlsx |

## Software Requirements

### R (version 4.4.2)
```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("dorothea")
install.packages(c("readxl", "writexl", "openxlsx", "ggplot2",
                   "dplyr", "tidyr", "svglite", "extrafont",
                   "igraph", "stringr"))
```

### Julia (version 1.11.4)
```julia
using Pkg
Pkg.add(["DifferentialEquations", "BlackBoxOptim", "Plots",
         "StatsPlots", "DataFrames", "CSV", "Statistics",
         "Printf", "LinearAlgebra", "IJulia"])
```

### Jupyter (for Julia notebooks)
```bash
pip install jupyterlab
```

Then in Julia:
```julia
using IJulia
IJulia.installkernel("Julia 1.11.4")
```

## How to Run

### Step 4a — Pharmacokinetic Modeling (R)
```r
setwd("/path/to/your/data/folder")
source("VitD_Modeling_script.R")
```

Fits exponential decay C(t) = C0·exp(-d·t) and supplementation ODE
dC/dt = k - d·C to serum 25(OH)D3 data across three replicates.
The ki and di values produced here are used as fixed inputs in
GRN_Modeling_Final.ipynb and are not optimised further.

### Step 4b — Pharmacokinetic Verification (Julia)
```bash
cd /path/to/your/data/folder
jupyter lab
```

Open `VitaminD_modeling.ipynb` → Kernel → Restart Kernel and Run All Cells

Numerically verifies the PK model using DifferentialEquations.jl.
No input files needed — all parameters are hardcoded.

### Step 5 — GRN Construction and Kernel Reduction (R)
```r
setwd("/path/to/your/data/folder")
source("Dorothea_Network_and_Kernel_Reduction.R")
```

Constructs a directed GRN from DoRothEA (confidence levels A and B)
anchored on 7 seed transcription factor genes. Applies iterative kernel
reduction until the network stabilises at 11 genes and 117 interactions
(71 activations, 46 repressions). Outputs network figures, interaction
matrices, and Hill-type equation files.

### Step 6 — ODE-Based GRN Modeling (Julia)
```bash
cd /path/to/your/data/folder
jupyter lab
```

Open `GRN_Modeling_Final.ipynb` → Kernel → Restart Kernel and Run All Cells

Implements the 11-gene Hill-equation ODE system. Fits 373 free parameters
per replicate for 4 models (full cohort, stable high, stable mid, stable low)
using BlackBoxOptim.jl differential evolution. VitD ki and di are fixed
from Step 4a and are not optimised.

⚠ Runtime warning: optimisation is computationally intensive.
Expect several hours per model depending on your hardware.
Do not close the notebook while it is running.

### Step 7 — Parametric Analysis (R)
```r
setwd("/path/to/your/data/folder")
source("Parametric_Analysis.R")
```

Reads the four optimized_equations_*.txt files. Computes composite gene
importance scores (degree centrality + regulatory strength + VitD weighting)
and pairwise parameter ratios between high and low responders.

## 11-Gene Reduced Network

Genes in the ODE model:
ETS2, FOS, HIF1A, EGR1, ETS1, JUN, MYC, NFKB1, STAT1, STAT3, TP53
+ Vitamin D (dynamic variable with parameters fixed from Step 4a)

- Total regulatory interactions: 117 (71 activations, 46 repressions)
- Free parameters per replicate per model: 373
- Models fitted: full cohort, stable high, stable mid, stable low

## Data Availability

Input data files are available from the corresponding author upon
reasonable request. RNA-seq data will be deposited in the NCBI Gene
Expression Omnibus (GEO) upon acceptance of the manuscript.

Pharmacokinetic reference data from:
Rybiński M, Ghosh Dastidar R, Zawrotna N et al.
Sci Rep 16, 2997 (2026). https://doi.org/10.1038/s41598-025-32831-z

## Citation

[Citation to be added upon publication]

## License

MIT License. See LICENSE file for details.
