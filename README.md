# HyTT: A hybrid machine learning and enzyme-constrained metabolic model for ab initio prediction of proteome reallocation

HyTT (**Hy**brid **T**ranslation-**T**ranscription) is a novel, *ab initio* computational platform that integrates sequence-driven machine learning predictions with enzyme-constrained genome-scale metabolic models (ecGEMs). By embedding Multivariate Adaptive Regression Splines (MARS) penalties as strict constraints within a mixed-integer linear programming (MILP) formulation, HyTT autonomously couples macroscopic spatial capacity limits directly to microscopic, sequence-derived translational costs. This framework enables the accurate prediction of dynamic resource reallocation, metabolic burden, and complex systemic adaptations, such as the Crabtree effect, overflow metabolism, and targeted ribosomal paralog switching, without relying on condition-specific multi-omics inputs.

---

## Repository Structure

- `src/`: Core MATLAB functions for constructing and solving the HyTT MILP framework.
- `data/`: Curated biological feature sheets, sequence traits, and ribosomal stoichiometry matrices.
- `tutorial_main.m`: A comprehensive demonstration script illustrating how to run the full pipeline and simulate overflow metabolism and proteome reallocation.

---

## Prerequisites & Installation

To run the HyTT framework, you need to ensure the following dependencies are installed and properly configured in your MATLAB environment:

1. **MATLAB (R2022a or later recommended)**
2. **COBRA Toolbox:** For constraint-based reconstruction and analysis. [Installation Guide](https://opencobra.github.io/cobratoolbox/stable/installation.html)
3. **GECKO Toolbox (v3.0):** Required for handling enzyme-constrained models (`ecYeastGEM`). [GECKO GitHub](https://github.com/SysBioChalmers/GECKO)
4. **MILP Solver:** Gurobi (v9.0+ recommended) or CPLEX, configured as the default factory solver for the COBRA toolbox.

---

## Getting Started & Usage

### 1. Clone the repository
First, download the framework to your local machine:
```bash
git clone [https://github.com/ehsanmotamedian/HyTT.git](https://github.com/ehsanmotamedian/HyTT.git)
cd HyTT