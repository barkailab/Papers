# Gera2021
This repository contains all scripts necessary to reproduce the figures and supplementary figures of Gera et al. 2022:
"Evolution of binding preferences among whole-genome duplicated transcription factors".

The matlab version used to run the scripts was 2021 a and b.
In addition to the scripts and functions in this repository, Mat files containing the processed sequencing data are available on dryad (), and should be contained in the main folder of this repository when running the scripts.
GithubDesc file contains a list and short description of all files (scripts from GitHub and data files from dryad) comprising this repository.

The main seuencing data can be divided into 3 types: 
- genome-wide profiles (norm): the normalised read densities along all S cerevisiae bases (~12 mio) - mean for all strains, indidividual repeats for wild-type strains
- promoter binding signal (sumProm): the cumultative read densities over each yeast promoter in the 6701-orf format - mean and individual repeats for all strains.
- 7mer binding signal: the average read density around (20bp) all promoter occurences of each possible 7mer (8192 7mers)- mean and individual repeats for all strains.

The order of the 6701-orf and 8192-7mer format can be found in gene_infoR64, part of group_imp.mat, and nmerRed table, part of nmer7N0.mat, respectively.
Several functions in this repository, e.g. Violin.m, brewermap.m and cbrewer2.m, were taken from Matlab File Exchange.
We thank the respective authors for sharing these files with the wider community - and ask them to contact us if there are any issues with including those files in our repository.
We hope that all functions are sufficiently commented to be accesible for researchers with basic Matlab experience,
but are eager to help via eMail (fxjonas@gmail.com) if any issues arise.

We thank you for your interest in our research and looking forward to your feedback.

Tamar and Felix
