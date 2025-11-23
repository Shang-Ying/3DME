# 3DME – Meshfree 3D Moment-Equation Solver for Subsurface Flow

This repository contains MATLAB scripts implementing a three-dimensional meshfree generalized finite difference method (GFDM) for steady subsurface flow.  
The code solves both the deterministic groundwater flow equation and the associated first- and second-order moment equations, and includes utilities for Monte Carlo reference solutions via sequential Gaussian simulation.

The scripts were originally developed to generate the numerical results in an accompanying research manuscript on efficient uncertainty quantification for 3-D subsurface flow.

---

## Features

- Meshfree generalized finite difference method (GFDM) in 3-D
- Steady-state deterministic groundwater head solver
- First- and second-order moment-equation solvers for hydraulic head
- Unconditional covariance matrix assembly for Gaussian random fields
- Sequential Gaussian simulation (2-D and 3-D) for Monte Carlo reference solutions
- Node generation / thinning utility for irregular 3-D domains
- Open-source MATLAB implementation under the GPL-3.0 license

---

## Repository contents

All files are plain MATLAB `.m` scripts/functions:

| File | Description |
| ---- | ----------- |
| `GP_MEs.m` | Driver script for solving the steady 3-D moment equations for hydraulic head using the meshfree GFDM formulation. |
| `GP_MCS.m` | Driver script for Monte Carlo simulations with random hydraulic conductivity fields, used as a reference solution. |
| `fct_GFDM_Coef_3D.m` | Builds the 3-D GFDM stencil coefficients for a set of nodes in an irregular domain. |
| `fct_Steady_DeterHead_3D_Acc.m` | Solves the deterministic steady-state head equation in 3-D using the GFDM coefficients. |
| `fct_Steady_1stMEHead_3D.m` | Assembles and solves the first-order moment equation for the hydraulic head. |
| `fct_Steady_2ndMEHead_3D.m` | Assembles and solves the second-order (co)variance / second moment equation for the hydraulic head. |
| `fct_UnCovMat_3D.m` | Constructs the unconditional covariance matrix of the log-conductivity (or conductivity) field in 3-D. |
| `node_drop_3d_GPT.m` | Utility for generating and/or thinning 3-D computational nodes (e.g. advancing-front / node-dropping for irregular geometries). |
| `sgsim.m` | Sequential Gaussian simulation in lower dimensions (e.g. 2-D). |
| `sgsim_3d_fast.m` | Fast 3-D sequential Gaussian simulation used to generate random fields for Monte Carlo tests. |
| `LICENSE` | GPL-3.0 license file for this repository. |

> **Note:** The roles above are indicative. Please see the header comments in each file for the exact inputs, outputs, and options.

---

## Getting started

### Requirements

- MATLAB (recent versions should work; the code is written in standard MATLAB syntax).
- No external libraries beyond MATLAB and its standard toolboxes are required.

### Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/Shang-Ying/3DME.git
2. Add the folder to your MATLAB path (e.g. using addpath or the MATLAB path manager).

## Basic usage

### 1. Moment-equation solver (`GP_MEs.m`)

This script solves:

- Deterministic steady-state hydraulic head,  
- First-order moment (mean) of head, and  
- Second-order moment / variance of head  

for a 3-D irregular domain, using the meshfree GFDM formulation.

Typical workflow:

1. Open `GP_MEs.m` in MATLAB.  
2. Review and, if needed, modify the problem settings in the script, such as:
   - Domain size and node distribution  
   - Boundary conditions  
   - Statistical description of hydraulic conductivity (mean, variance, correlation length, etc.)  
   - Choice of example geometry (`*.dat` files)  
3. Run the script (e.g., press **Run** in the MATLAB editor).

The script will call:

- `fct_GFDM_Coef_3D.m`  
- `fct_Steady_DeterHead_3D_Acc.m`  
- `fct_Steady_1stMEHead_3D.m`  
- `fct_Steady_2ndMEHead_3D.m`  
- `fct_UnCovMat_3D.m`  

and produce the deterministic head and its first- and second-order moments, together with any plots or data saved by the script.

---

### 2. Monte Carlo reference simulations (`GP_MCS.m`)

This script generates Monte Carlo reference solutions using sequential Gaussian simulation.

Typical workflow:

1. Open `GP_MCS.m` in MATLAB.  
2. Adjust the Monte Carlo settings:
   - Number of realizations  
   - Variogram / covariance parameters  
   - Grid / node configuration  
   - Choice of example geometry (`*.dat` files)  
3. Run the script to:
   - Generate multiple random hydraulic conductivity fields (via `sgsim_3d_fast.m` and/or `sgsim.m`),  
   - Solve the deterministic flow equation for each realization, and  
   - Compute sample statistics (e.g., sample mean and variance of hydraulic head).

The resulting Monte Carlo statistics can be compared against the moment-equation results from `GP_MEs.m`.

---

## Example geometries

The four `*.dat` files included in the repository provide example irregular 3-D domains that are used in the test cases. In general:

- The driver scripts read these files to construct the computational domain for both the moment-equation solver and the Monte Carlo solver.

To use your own geometry:

1. Prepare your own `*.dat` files (or other MATLAB-readable format) containing node coordinates and connectivity information compatible with the existing scripts.  
2. Modify the corresponding file paths and loading routines in `GP_MEs.m` and/or `GP_MCS.m` so that they point to your geometry files.

---

## Third-party component: `intriangulation.m`

This repository bundles `intriangulation.m`, a MATLAB function by **Johannes Korsawe** from MATLAB Central File Exchange (Version 1.5.0.0), which determines whether 3-D test points lie inside or outside a watertight triangular surface mesh.

- Original source: MATLAB Central File Exchange –  
  *“intriangulation(vertices,faces,testp,heavytest)” by Johannes Korsawe*  
- The code in `intriangulation.m` is included here for convenience and is used to perform point-in-polyhedron checks when constructing or post-processing irregular 3-D domains.  
- **License:** `intriangulation.m` is distributed under its own license terms, which remain unchanged. Please refer to the license and copyright
  notice inside `intriangulation.m` and on the MATLAB File Exchange page for the exact terms.  
- If you use `intriangulation.m` directly in your own work, please cite it as recommended by the MATLAB Central File Exchange entry.

---

## Reproducibility

The repository is intended to accompany the numerical examples in the associated research manuscript. To reproduce specific figures:

1. Identify which driver script and case were used (typically documented in comments in `GP_MEs.m` and `GP_MCS.m`).  
2. Set the parameters in the script to match the manuscript (geometry, boundary conditions, statistics of hydraulic conductivity, number of realizations, etc.).  
3. Run the script and post-process the outputs (plots, statistics) as needed.

---

## Citation

If you use this code in your research, please cite the accompanying manuscript (once published) and/or this repository. A generic repository citation is:

> Chen, S.-Y. (2025). *3DME – Meshfree 3D Moment-Equation Solver for Subsurface Flow* (MATLAB code). GitHub repository, https://github.com/Shang-Ying/3DME.

If you use `intriangulation.m` directly, please also cite the original MATLAB Central File Exchange entry for `intriangulation(vertices,faces,testp,heavytest)` by Johannes Korsawe.

---

## License

Unless otherwise noted (see below), the code in this repository is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.  
See the `LICENSE` file for the full license text.

- **3DME code (all files except explicitly marked third-party code)**  
  © 2025, Shang-Ying Chen and collaborators.  
  Licensed under GPL-3.0.

- **Third-party code: `intriangulation.m`**  
  © Johannes Korsawe and contributors.  
  Included under its own license terms, which remain in effect and are not overridden by the GPL-3.0 license of this repository.  
  Please refer to the license information inside `intriangulation.m` and on its MATLAB Central File Exchange page.

By using this repository, you agree to comply with **both** the GPL-3.0 license for 3DME and the original license of `intriangulation.m`.

---

## Contact

For questions, comments, or bug reports, please open an issue in this GitHub repository.

If you have questions related to the underlying research manuscript, please contact the corresponding author as indicated in the paper.

