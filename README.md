
# Fishermathica


A Mathematica code to perform Fisher Matrix forecasts for photometric and spectroscopic galaxy surveys.

This is a fresh commit github to avoid carrying around a lot of commit history and old tracked files. The full commit history (570 commits) can be found on 
https://gitlab.com/santiagocasas/FisherMatrixTools.

## Branch 'alphaAttractors'
This branch was developed during the work behind the publication **Quintessential α-attractor inflation: forecasts for Stage IV galaxy surveys** (doi:[10.1088/1475-7516/2021/04/006](https://dx.doi.org/10.1088/1475-7516/2021/04/006)).

It contains
- New packages: [Quintessence.wl](Quintessence.wl) and [Jacobian.wl](Jacobian.wl),
- Adaptation to [CosmologyFunctions.mCosmologyFunctions.m),
- Self-contained scripts: [SelfContainedRun](SelfContainedRun/).

By running the code, one can reproduce the Fisher information leading to the reported figures in the publication.

### Usage
1. Run ``*generate-script.m`` in [SelfContainedRun/InputFilesGenerator/](SelfContainedRun/InputFilesGenerator/),for the cosmological model of choice.
2. Run ``SURVEYNAME-MODELNAME.m`` in [SelfContainedRun](SelfContainedRun/) to obtain the Fisher matrices for the survey and the model of choice.