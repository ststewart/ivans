# Planetesimal Impact Vapor Plumes and Nebular Shocks form Chondritic Mixtures

This repository presents supplemental materials and replication files for the Impact Vapor and Nebular Shocks (IVANS) model for the formation of chondrules and chondritic mixtures.<p>

Reference:<br>
S. T. Stewart, S. J. Lock, P. J. Carter, E. J. Davies, M. I. Petaev, S. B. Jacobsen, Planetesimal Impact Vapor Plumes and Nebular Shocks form Chondritic Mixtures, <i>[The Planetary Science Journal](https://iopscience.iop.org/journal/2632-3338)</i>, in press, 2025.

The revised supplementary Jupyter book (v1.0.2) is archived at <br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14969068.svg)](https://doi.org/10.5281/zenodo.14969068)

This Jupyterbook is hosted online at https://chondrules.net/ivans.

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

