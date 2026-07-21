# HyTT: A hybrid machine learning and enzyme-constrained metabolic model for ab initio prediction of proteome reallocation

HyTT (**Hy**brid **T**ranscription-**T**ranslation) is a novel, *ab initio* computational platform that integrates sequence-driven machine learning predictions with enzyme-constrained genome-scale metabolic models (ecGEMs). By embedding Multivariate Adaptive Regression Splines (MARS) penalties as strict constraints within a mixed-integer linear programming (MILP) formulation, HyTT autonomously couples macroscopic spatial capacity limits directly to microscopic, sequence-derived translational costs. This framework enables the accurate prediction of dynamic resource reallocation, metabolic burden, and complex systemic adaptations - such as the Crabtree effect, overflow metabolism, and targeted ribosomal paralog switching - without relying on condition-specific multi-omics inputs.

---

## 📁 Repository Structure

- `src/`: Core MATLAB functions for constructing and solving the HyTT MILP framework.
- `data/`: Curated biological feature sheets, sequence traits, and ribosomal stoichiometry matrices.
- `tutorial_main.m`: A comprehensive demonstration script illustrating how to run the full pipeline and simulate overflow metabolism and proteome reallocation.

---

## ⚙️ Prerequisites & Installation

To run the HyTT framework, you need to ensure the following dependencies are installed and properly configured in your MATLAB environment:

1. **MATLAB (R2022a or later recommended)**
2. **COBRA Toolbox:** For constraint-based reconstruction and analysis. [Installation Guide](https://opencobra.github.io/cobratoolbox/stable/installation.html)
3. **GECKO Toolbox (v3.0):** Required for handling enzyme-constrained models (`ecYeastGEM`). [GECKO GitHub](https://github.com/SysBioChalmers/GECKO)
4. **MILP Solver:** Gurobi (v9.0+ recommended) or CPLEX, configured as the default factory solver for the COBRA toolbox.

---

## 🚀 Getting Started & Usage

### 1. Clone the repository
First, download the framework to your local machine:

```bash
git clone https://github.com/ehsanmotamedian/HyTT.git
cd HyTT
```

### 2. Initialization
Open MATLAB and navigate to the `HyTT` directory. Add the repository to your MATLAB path and initialize the required toolboxes:

```matlab
% Add HyTT folders to path
addpath(genpath(pwd));

% Initialize COBRA and GECKO (Ensure they are already installed)
initCobraToolbox(false);
```

### 3. Running the Simulation
To demonstrate the capabilities of the HyTT framework, we have provided a comprehensive tutorial script. Simply run the following command in the MATLAB command window:

```matlab
tutorial_main
```

This script will guide you through:
- Loading the sequence-derived parameters and the base ecGEM.
- Applying the MARS-derived translational penalties.
- Formulating and solving the MILP problem to predict dynamic resource reallocation.
- Visualizing the onset of overflow metabolism (Crabtree effect).

---

## 📖 How to Cite

If you use the HyTT framework in your research, please cite our preprint:

> Motamedian, E., & Nikoloski, Z. (2026). *A hybrid machine learning and enzyme-constrained metabolic model for ab initio prediction of proteome reallocation*. bioRxiv. [DOI will be updated upon publication]

---

## 📄 License
This project is licensed under the MIT License - see the `LICENSE` file for details.

## ✉️ Contact
For questions, bug reports, or collaboration inquiries, please open an issue in this repository or contact the corresponding authors.
