# SEM_QTL
Documentation of code for SEM-based QTL analyses performed with the lavaan R package.  

The multivariate QTL mapping algorithm treats identification of QTL as updating of a genotype-phenotyping graph with markers affecting phenotypes out of the total population of loci genome-wide. The code is written in R and uses the lavaan sem() function, so is pretty slow. 

ycMapping_dataProcess.R contains code for formatting of genotype and phenotype data on plant growth and yield component phenotypes collected in the LA95135 x SS-MPV57 mapping population

graphConstruction.R contains code for the multivariate QTL mapping algorithm used to identify yield component QTL correcting for plant growth phenotypes, while graphConstructionFunctions.R contain the functions called by graphConstruction.R.

graphConstruction.sh contains code for running the mapping script on an HPC using the Slurm job scheduler.
