# Pyrobombus-speciation
This repository contains scripts written for the analysis carried out in Christmas et al. 2020 "Cryptic speciation with gene flow in alpine bumblebees revealed by comparative population genomics"

All custom scripts were written in perl. I wrote these to work specifically on our datasets and so will require some editing (such as input/output filenames and sample names) if you are to apply them to other datasets.<br>

Proportion of GC in non-overlapping slidding windows - [count_GC_in_windows.pl](count_GC_in_windows.pl) <br>
Parsing GenMap output in non-overlapping slidding windows - [mappability_in_windows.pl](mappability_in_windows.pl) <br>
Parsing RepeatMasker output in non-overlapping slidding windows - <br>
[parse_repeatmasker_out.pl](parse_repeatmasker_out.pl) <br>
[get_window_measure_of_repeat_content.pl](get_window_measure_of_repeat_content.pl) <br>
Parsing LDHat stat output in non-overlapping slidding windows - [get_average_rho_over_20kb_windows.pl](get_average_rho_over_20kb_windows.pl)<br>
Calculating nucleotide diversity (π) - [window_pi_calcs.pl](window_pi_calcs.pl)<br>
Calculating absolute divergence (d<sub>XY</sub>) - [window_dxy_calcs.pl](window_dxy_calcs.pl)<br>
Caclulating the average changes in 20Kbp steps away from the centres of islands of divergence:  <br>
ZF<sub>ST</sub> - [calc_average_zfst_from_peaks.pl](calc_average_zfst_from_peaks.pl) <br>
π  - [calc_average_pi_from_peaks.pl](calc_average_pi_from_peaks.pl) <br>
d<sub>XY</sub> - [calc_average_dxy_from_peaks.pl](calc_average_dxy_from_peaks.pl) <br>
<br>
The majority of the analyses in the paper (identifying islands of divergence (IoDs), comparing metrics inside and outside of IoDs, and all plotting) were performed in R. Theses analyses are detailed in the R script [Pyrobombus_speciation_genomics_analysis.R](Pyrobombus_speciation_genomics_analysis.R)
