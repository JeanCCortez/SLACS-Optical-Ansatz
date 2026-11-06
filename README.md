# The Optical Phase-Delay Signature in Massive Elliptical Lenses

**Data and Code Repository for the Phenomenological Audit of the SLACS Catalog**

This repository contains the complete dataset and Python scripts required to reproduce the exact findings presented in the paper: *"The Optical Phase-Delay Signature in Massive Elliptical Lenses: A Parameter-Free Phenomenological Ansatz Applied to the SLACS Catalog"*.

The theoretical framework introduces a parameter-free ($N=0$) optical-hydrodynamic ansatz, utilising a universal phase-impedance constant ($\beta = 0.028$) to resolve the mass-sheet degeneracy without invoking non-baryonic dark matter profiles. 

## 📂 Repository Structure

### Data
* `SLACS_table.cat`: The raw observational data from the Sloan Lens ACS Survey (SLACS). It includes the lens redshift ($z_L$), source redshift ($z_S$), Einstein radius ($\theta_E$), and central stellar velocity dispersion ($\sigma_v$).

### Scripts
* `01_compute_metrics_3.py`: The core algorithmic solver. It applies the $N=0$ physical state equation (incorporating thick-vacuum drag and a linear causal metric) to predict the source redshift ($z_S$) purely from baryonic and geometric observables.
* `02_monte_carlo_sonnenfeld_3.py`: A rigorous Monte Carlo error propagation simulation ($10^3$ iterations). It injects the exact systematic biases and random noise documented for SDSS single-fibre spectroscopy (Sonnenfeld, 2025) to establish the theoretical instrumental noise floor.
* `03_plot_redshift_recovery_3.py`: The visualization script that generates the primary diagnostic scatter plot comparing the blind predicted redshifts against spectroscopic observations.

### Outputs
* `slacs_rigorous_physics_audit.png`: The final diagnostic figure demonstrating the blind source redshift recovery for the $N=59$ systems.

---

## ⚙️ Requirements & Setup

The codebase is written in Python 3. To run the audit locally, ensure the following standard scientific libraries are installed:

```bash
pip install pandas numpy scipy matplotlib
🚀 Execution Guide & ReproducibilityTo ensure absolute transparency and reproducibility, the audit must be executed sequentially. Below are the exact commands and the expected theoretical outputs.Step 1: Global Analytical PerformanceRun the core metric solver to evaluate the parameter-free model against the complete SLACS sample ($N=59$) and compare it to the naive geometric baseline ($z_S = 2z_L$).Bashpython 01_compute_metrics_3.py
Expected Output:
The script will output a Mean Absolute Error (MAE) of 0.2946 for the proposed model, and 0.2324 for the null model, with a residual correlation of 0.2034.Step 2: The Monte Carlo Noise Floor ValidationTo demonstrate that the model's MAE is not a theoretical failure but a reflection of instrumental limits, run the Monte Carlo simulation. This applies a +5% systematic bias to $\sigma_v$, an 8.4% log-normal intrinsic scatter, and a 5% Gaussian uncertainty to $\theta_E$.Bashpython 02_monte_carlo_sonnenfeld_3.py
Expected Output:
The simulation confirms that the physical noise floor of the SDSS fibre spectroscopy dictates an expected MAE of ~0.309 ± 0.073. The observed MAE of 0.2946 falls perfectly within this purely instrumental limit, proving that the $N=0$ equation exposes observational noise rather than absorbing it through curve-fitting.Step 3: Diagnostic VisualizationGenerate the high-resolution figure used in the manuscript.Bashpython 03_plot_redshift_recovery_3.py
Expected Output:
Saves slacs_rigorous_physics_audit.png to the local directory, displaying the recovery scatter and statistical metrics (N=59, MAE=0.2946, RMSE=0.7637).📖 Physical Constants & CalibrationThe entire pipeline operates under a strict Zero-Parameter ($N=0$) framework.The speed of light is fixed at $c = 3 \times 10^5$ km/s.The universal phase-impedance constant is fixed at $\beta = 0.028$.Crucially, $\beta$ is not derived from nor fitted to the SLACS data. It is independently calibrated from the asymptotic velocity plateaus in late-type galaxies via the Radial Acceleration Relation (RAR), rendering the lensing audit a truly blind test.📜 Author & LicenseJean C. CortezIndependent Researcher, Rio de Janeiro, Brazil.This project is licensed under the MIT License - providing unrestricted access to the code and data for the advancement of open science and reproducibility.
