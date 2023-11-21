# enriched_OG
This pipeline allows to identify  enriched **Orthologous Groups (OG)** in the internal nodes of a rooted phylogenetic tree.

*enriched_OG* will use a rooted-ultrametric phylogenetic tree in newick format with node names **(see example file)** as well as a presence/absence matrix **(see example file)** as input files. *enriched_OG* will create a stochastic mapping distribution and a null distribution. These distributions can be statistically compared (null distribution is generated without the presence/absence data at the tips) to find OG enriched at specific nodes of the input phylogenetic tree. 

For more information on the stochastic mapping analysis:  

http://www.phytools.org/static.help/make.simmap.html

## STEP 1: CREATING A STOCHASTIC MAPPING AND A NULL DISTRIBUTION ##
R code: 1_stochastic_mapping_null.R

**Libraries required:**  

ape (Paradis and Schliep, 2019)  

phytools (Revell, 2012)  

tidytree (Yu, 2022)  

hash (Brown, 2023)  

This script will perform a stochastic mapping simulation and a null distribution of the orthologous groups gained and lost at each node of the input phylogenetic tree **(replicates to be decided by the user; modify line 11, default value is 100)**. The user can decide what model to use for the stochastic mapping analysis (line 36; see make.simmap documentation https://search.r-project.org/CRAN/refmans/phytools/html/make.simmap.html). The default model is ARD ("all-rates-different)".

**Note: the Q matrix will be modified if one of the lines has 0 as value. In such a case, the replaced value on the Q matrix will be 0.00001**

## STEP 2: PARSING STOCHASTIC MAPPING AND NULL DISTRIBUTION RESULTS ##

Python code: 2_stochastic_mapping_null_parsing.py

How to run:  

python 2_stochastic_mapping_null_parsing.py input_folder output_file

input_folder: location of the output files resulting from 1_stochastic_mapping_null.R  

output_file: file with the stochastic mapping and the null distributions. The column Clusters_gained_or_lost represents the orthologous groups gained or lost at each node on the tree and for each simulation.  

The resulting stochastic and null distributions can be statistically compared (e.g., Wilcox test) to find enriched orhologous groups at nodes of interest.

**This code is part of the data analysis on the following study:**  

**Martinez-Gutierrez, C. A., Uyeda, J. C., & Aylward, F. O. (2022). A timeline of bacterial and archaeal diversification in the ocean. bioRxiv, 2022-10.**  

## REFERENCES ##

ape library: Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526-528. 

phytools library: Revell, L. J. (2012). phytools: an R package for phylogenetic comparative biology (and other things). Methods in ecology and evolution, (2), 217-223. 

tiditree library: https://github.com/YuLab-SMU/tidytree 

hash library: https://cran.r-project.org/web/packages/hash/index.html 


