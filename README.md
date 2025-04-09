# Planetesimal Impact Vapor Plumes and Nebular Shocks form Chondritic Mixtures

This repository presents supplemental materials and replication files for the Impact Vapor and Nebular Shocks (IVANS) model for the formation of chondrules and chondritic mixtures.<p>

Reference:<br>
S. T. Stewart, S. J. Lock, P. J. Carter, E. J. Davies, M. I. Petaev, S. B. Jacobsen, Planetesimal Impact Vapor Plumes and Nebular Shocks form Chondritic Mixtures, <i>[The Planetary Science Journal](https://iopscience.iop.org/journal/2632-3338)</i>, doi:10.3847/PSJ/adbe71, 2025.

The supplementary materials are contained in a Jupyterbook that is archived at <br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14969068.svg)](https://doi.org/10.5281/zenodo.14969068)

This Jupyterbook is hosted online at https://chondrules.net/ivans.

This repository contains:<br>
* Jupyter notebooks (ipynb files) and python scripts (py files) to recreate Figures 1-2, 8-16, and 18-24
    - Data files in text, CSV and python pickle formats
    - The required python environment yaml file
* Table 3 in CSV format
* Example text input files for the CTH code
* Example text input files for the pyKO code
* Equation of state files for ANEOS water and ANEOS fused silica
    - Text input files, text output files, and text SESAME format data tables
* Supplemental animations
    - HTML animation associated with Figure 5 pyKO code calculation
    - QuickTime movie files associated with Figures 6 and 7 CTH code calculations
* Additional files needed to assemble the Supplemental Materials into a Jupyter book

To execute the Jupyter notebooks in this GitHub repository, build the conda environment:<p>

`conda env create -f environment-ivans.yml`<br>
`conda activate ivans`<p>

To rebuild the html Jupyterbook, from the directory above the `ivans` full repository, run:<p>

`jupyter-book build ivans`<p>

The output is in:<p>

`ivans/_build/html`<p>

Move the contents of the `ivans/static` subdirectory to `ivans/_build/html/_static` to complete the Jupyter book. See `ivans/_toc.yml` for the table of contents of the Jupyterbook.<p>

S. T. Stewart<br>
sstewa56@asu.edu<br>
Submitted 11/29/2024<br>
Revised 3/6/2025<br>
Accepted 3/7/2025<br>
Supplemental materials corrections v1.0.3 4/8/2025<br>
