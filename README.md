# HyTT: Hybrid Translation-Transcription Framework for Resource Allocation Modeling

HyTT (**Hy**brid **T**ranslation-**T**ranscription) is a novel, *ab initio* mechanistic-statistical computational platform that structurally integrates machine learning predictions with enzyme-constrained genome-scale metabolic models (ecGEMs). By utilizing Multivariate Adaptive Regression Splines (MARS), HyTT autonomously couples macroscopic spatial proteomic boundaries with microscopic, sequence-derived translational costs to simulate cellular resource allocation and predict complex phenotypes like the "proteome squeeze" and metabolic re-routing without relying on condition-specific omics inputs.

---

## Repository Structure

- `src/`: Core MATLAB functions for constructing and solving the HyTT MILP framework.
- `data/`: Curated biological feature sheets, sequence traits, and ribosomal stoichiometry matrices.
- `tutorial_main.m`: A comprehensive demonstration script illustrating how to run the full pipeline and simulate overflow metabolism (the Crabtree effect).

---

## Prerequisites & Installation

To run the HyTT framework, you need to ensure the following dependencies are installed and properly configured in your MATLAB environment:

1. **MATLAB (R2022a or later recommended)**
2. **COBRA Toolbox:** For constraint-based reconstruction and analysis. [Installation Guide](https://opencobra.github.io/cobratoolbox/stable/installation.html)
3. **GECKO Toolbox (v3.0):** Required for handling enzyme-constrained models (`ecYeastGEM`). [GECKO GitHub](https://github.com/SysBioChalmers/GECKO)
4. **MILP Solver:** Gurobi (v9.0+ recommended) or CPLEX, configured as the default factory solver for the COBRA toolbox.

### Getting Started

1. Clone this repository to your local machine:
   ```bash
   git clone [https://github.com/ehsanmotamedian/HyTT.git](https://github.com/ehsanmotamedian/HyTT.git)
   cd HyTT