# Pyrobombus-speciation
This repository contains scripts written for the analysis carried out in Christmas et al. 2020 "Cryptic speciation with gene flow in alpine bumblebees revealed by comparative population genomics"

All custom scripts were written in perl. I wrote these to work specifically on our datasets and so will require some editing (such as input/output filenames and sample names) if you are to apply them to other datasets.<br>

Proportion of GC in non-overlapping slidding windows - <br>
Parsing GenMap output in non-overlapping slidding windows - <br>
Parsing RepeatMasker output in non-overlapping slidding windows - <br>
Parsing LDHat stat output in non-overlapping slidding windows - <br>
Calculating nucleotide diversity (π) - <br>
Calculating absolute divergence (d<sub>XY</sub>) - <br>
Caclulating the average changes in ZF<sub>ST</sub>, π, and d<sub>XY</sub> in 20Kbp steps away from the centres of islands of divergence - <br>
<br>
The majority of the analyses in the paper (identifying islands of divergence (IoDs), comparing metrics inside and outside of IoDs, and all plotting) were performed in R. Theses analyses are detailed in the R script Pyrobombus_speciation_genomics_analysis.R
