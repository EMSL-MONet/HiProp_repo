# HiProp: Soil Hydraulic Properties and Soil C Cycle Modeling

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository contains the code and analysis for the manuscript **"Simulated Soil Respiration is Sensitive to Soil Hydraulic Properties from Intact versus Repacked Cores"**. The project investigates how different soil hydraulic properties measured from intact versus repacked soil cores affect simulated soil respiration rates using coupled hydrological and biogeochemical models.

## Description

Soil hydraulic properties play a critical role in controlling water movement, oxygen diffusion, and ultimately microbial respiration in soils. This study examines the sensitivity of modeled soil respiration to variations in soil hydraulic properties derived from intact versus repacked soil cores. The repository provides:

- Coupled soil hydrology and biogeochemical modeling framework
- Functions for simulating soil moisture dynamics and respiration
- Analysis tools for comparing model outputs across different hydraulic parameterizations
- Visualization and statistical analysis code

## Getting Started

### Prerequisites

- Python 3.x
- Jupyter Notebook
- Required Python packages (numpy, pandas, matplotlib, scipy, etc.)

### Running the Analysis

Use the [HiProp_manuscript.ipynb](HiProp_manuscript.ipynb) Jupyter notebook to reproduce all figures and analyses from the manuscript. The notebook includes:

- Data loading and preprocessing
- Model simulations with different hydraulic parameters
- Statistical analyses and sensitivity tests
- Figure generation for the manuscript

## Key Files

- `HiProp_manuscript.ipynb` - Main analysis notebook
- `soil_hydrology_module.py` - Soil water dynamics and hydraulic functions
- `millennial_module.py` - Biogeochemical model implementation
- `MEMES_functions.py` - Additional modeling utilities
- `soilpara_in_fit.txt` - Soil hydraulic parameters
- `daily_mean_rainfall.csv` - Input climate data

## Citation

If you use this code in your research, please cite:

```
Citation coming soon
```

## License

This project is licensed under the [MIT License](LICENSE). You are free to use, modify, and distribute this code for academic or commercial purposes with proper attribution.

## Contact

**Arjun Chakrawal**  
Email: arjun.chakrawal@pnnl.gov

Pacific Northwest National Laboratory (PNNL)

## Acknowledgments

This work was supported by the Environmental Molecular Sciences Laboratory (EMSL), a DOE Office of Science User Facility sponsored by the Biological and Environmental Research program. 
