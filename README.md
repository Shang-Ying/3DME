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
