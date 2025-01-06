# Combined DCE and DSC Leakage Correction
Matlab code for correcting contrast leakage in DSC MRI using parameters derived from DCE MRI.

# Introduction
This repository contains Matlab code to estimate and compensate for contrast agent (GBCA) leakage in dynamic susceptibility contrast (DSC) MRI. The correction leverages dynamic contrast-enhanced (DCE) MRI–derived vascular parameters (e.g. Ktrans, vc and ve) to accurately model the tissue concentration and remove leakage-induced biases in DSC. By combining DCE and DSC analyses, the method aims to yield a more reliable DSC-based estimate in brain tumors or other pathologies.

# Instruction
The main.m demonstrates a full workflow using example data:
1. Loads or references DCE and DSC signal ratio curves, plus relevant protocol parameters (TR, TE, flip angles, T1 values, etc.).
2. Calls fit_sigR_DCE to estimate vascular permeability parameters (Ktrans, vc and ve) from the DCE data.
3. Calls fit_sigR_DSC to estimate tissue relaxivity (r2*) and post-preload tissue T1, using the DCE-derived parameters.
4. Performs the final leakage correction on the DSC data, generating corrected Δ𝑅2∗ curves.

# Note
The optimized DCE-driven AIF was used for both DCE and DSC signal ratio model fit according to our previous papers. Please refer to (https://doi.org/10.1002/nbm.5038) and (https://doi.org/10.1002/nbm.5225). For large datasets, semi-automatic and automatic algorithm for AIF determination in DSC datasets might be applied. In that case, applying a linear (or nonlinear) conversion for transfering Δ𝑅2∗ measurment in DSC-AIF into concentration unit is neccessary prior to DSC model fit.
In our case, the temporal resolutions in DCE and DSC scan are the same. While using DCE and DSC datasets with different temporal resolution, an interpolation process is encouraged.

# Contact
For questions or additional guidance, please contact:
Chih-Hsien Tseng
Email: <r04548023(at)gmail.com>
If you use or adapt this code for research, we appreciate a citation to the associated paper or this GitHub repository. 
