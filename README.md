# Emergence and spread of D1.1 reassortant H5N1 viruses in the Americas

Code and analysis workflows associated with:

> [Crespo-Bellido et al. (2026). Emergence of D1.1 reassortant H5N1 avian influenza viruses in North America](https://www.biorxiv.org/content/10.64898/2025.12.19.695329v2.full)

---

## Overview

This repository contains scripts, workflows, and selected data used for:

- phylogenetic subsampling
- continuous phylogeographic reconstruction
- spatial diffusion analyses
- pairwise similarity analyses between groups
- flyway-based phylogenetic visualization
- generation of figures and statistics for the associated manuscript

---

## Data notes

### Sequence data

Some sequence data used in this study were obtained from GISAID and are subject to data-sharing restrictions.

For this reason:
-	sequence alignments and BEAST input files containing raw sequence data are not distributed
-	accession numbers and scripts are provided to allow authorized users to reconstruct the analyses

### Intermediate outputs

Large intermediate outputs (e.g., posterior tree samples, extensive reconstructions) have been excluded to keep the repository manageable.

The repository contains:
-	representative outputs
-	scripts required to regenerate results

### Large files

Shapefiles are stored using **Git LFS**.

After cloning:

  git lfs install
  git lfs pull




[![CC BY 4.0][cc-by-shield]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
