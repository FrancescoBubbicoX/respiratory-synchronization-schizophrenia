# Respiratory Synchronization and Body–Brain Dynamics in Schizophrenia

## Overview

This repository contains the MATLAB and R code used to preprocess physiological recordings (EEG, ECG, and respiration), derive body–brain measures, estimate transfer entropy, and perform the statistical analyses reported in the study of respiratory synchronization and body–brain dynamics in schizophrenia spectrum disorders.

Associated manuscript:
"Predicting Others Through the Body: Respiratory Synchronization and Body–Brain Dynamics in Schizophrenia"

## Related Resources

- **Data repository (Zenodo):** https://doi.org/10.xxxx/zenodo.xxxxxxx
- **Project and preregistration (OSF):** https://doi.org/10.17605/OSF.IO/KQ6YN

## Repository Structure

- **PreprocessingData.m**: Preprocessing of EEG, ECG, and respiration signals.

- **BodyBrain_TimeSeries_Processing.m**: Extraction of RR intervals, respiratory measures, and EEG band-power time series.

- **TransferEntropy_Analysis.m**: Computation of transfer entropy between heart, respiration, and EEG-derived signals, including cluster-based permutation testing and result export.

- **RespiratorySync_Analysis.m**: Computation and analysis of respiratory synchronization measures.

- **BehaviouralAndSurvey_Plotting.R**: Plotting and statistical summaries for behavioral and questionnaire data.

- **BodyBrain_Analysis.R**: Main R script for statistical analyses linking body–brain measures with behavioral outcomes.

- **BodyBrain_Analysis_Functions.R**: Helper functions used by `BodyBrain_Analysis.R`.

- **Functions/**: Custom MATLAB functions used throughout the preprocessing and analysis pipeline.

## Requirements

The analyses were developed and tested using:

### MATLAB

- MATLAB R2024a
- FieldTrip (tested with version 20250106)
- EEGLAB (tested with version 2022.0)
- ITS Toolbox v2.1

Some scripts may require additional MATLAB toolboxes depending on the analysis being performed.

All custom MATLAB functions developed for this project are included in the `Functions/` directory.

### R

The R scripts specify all required packages at the beginning of each file. Users should install the listed packages before running the analyses.

### Data

The datasets required to reproduce the analyses are available in the associated Zenodo repository.

## Analysis Workflow

### Behavioral analyses

1. Respiratory synchronization analyses (`RespiratorySync_Analysis.m`)
2. Behavioral and questionnaire plotting (`BehaviouralAndSurvey_Plotting.R`)

### Body–brain analyses

1. Physiological signal preprocessing (`PreprocessingData.m`)
2. Extraction of body–brain time series (`BodyBrain_TimeSeries_Processing.m`)
3. Transfer entropy estimation and cluster-based permutation analyses (`TransferEntropy_Analysis.m`)
4. Statistical analyses linking body–brain measures with behavioral outcomes (`BodyBrain_Analysis.R`)

## Citation

A formal citation will be added upon publication of the associated manuscript. Until then, please cite the Zenodo dataset and the OSF project when appropriate.

## License

Except where otherwise noted, the code in this repository is Copyright (c) 2026 Francesco Bubbico and is released under the MIT License.

Third-party files retain their original licenses.

## Third-party Code

This repository includes the third-party MATLAB function:

### plot_topography.m

Copyright (c) 2020, Víctor Martínez-Cagigal

Distributed under the BSD 3-Clause License.

Original source:
https://www.mathworks.com/matlabcentral/fileexchange/72729-topographic-eeg-meg-plot

The original copyright notice and license are retained in the source file and in the accompanying license document.
