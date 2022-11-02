# Comparison of co-expression networks between species
code for co-expression analysis from the paper "Independent evolution of transcript abundance and regulatory dynamics". 

# Pipeline
The script "compare_PCMs.m" is a pipeline that runs the other functions, coded in matlab.
# Input file
The input is a table of normalized read counts in excel, where each sheet contains gene expression samples collected from a single species.
# Test input
"expression_per_species_test.xlsx"<br/>
This is a subset of the data presented in the article. Full data will be avilable soon.
# Output files
"ggc_(today's_date).mat" is a matlab struct that contains:<br/>
PCM: pairwise correlation matrix, per speices<br/>
expCorr: a bootstrap analysis for each dataset (dataset-control)<br/>
obsT: a table of comparative regulatory similarity score, for all pssible species comparisons<br/>
expT: a table of the median dataset-control (of 10 permutations), per species<br/>
nnT: nearest neighbors, per species. The pair of genes with the most similar co-expression vector<br/>
nnCrossSpT: cross species nearest neighbors. Taking the top n (here n=5) neighbors of gene i in species A, and pick the neighbor that shows maximal similarity to gene i in species B.<br/>
# Pre-requisits
Matlab
# Citing
Independent evolution of transcript abundance and gene regulatory dynamics
Gat Krieger, Offir Lupo, Avraham A. Levy, Naama Barkai
bioRxiv 2020.01.22.915033; doi: https://doi.org/10.1101/2020.01.22.915033

