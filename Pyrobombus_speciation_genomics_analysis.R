library(ggpubr)
library(cowplot)
library(reshape2)
library(ggsignif)
library(readr)
library(ggplot2)
library(dplyr)
library(boot)
library(ggpointdensity)
library(viridis)

# Write a function to calculate mean from the data

samplemean <- function(x, d) {  # This is the function to calculate mean from the data (x), using a bootstrap sample, d.
  return(mean(x[d]))
}

#Import all 20 kbp window datasets

# Depth of coverage
DoC_Bsyl <- read_delim("DoC_Bsyl_windows.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
DoC_Bsyl$BIN_START= DoC_Bsyl$BIN_START+1

# Contig genome positions
contig_lengths <- read_delim("Bsyl_contig_cum_lengths_ordered_by_Bt_LG.txt","\t", escape_double = FALSE, trim_ws = TRUE)

# Chromosome end coordinates
chrom_lengths <- read_delim("Bt_LG_end_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
# GC content
gc_20kb <- read_delim("Bsyl_gc_content.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(gc_20kb)[names(gc_20kb) == "genome_mid"] <- "GEN_MID"
gc_20kb$BIN_START= gc_20kb$BIN_START+1

LDhat_results <- read_delim("Bsyl_whole_genome_LDhat_stat_out_bpen_1_20kb_windows.txt","\t", escape_double = FALSE, trim_ws = TRUE)
#LDhat_results_bpen5 <- read_delim("Bsyl_whole_genome_LDhat_stat_out_bpen_5_20kb_windows.txt","\t", escape_double = FALSE, trim_ws = TRUE)

genmap_k150_e2_20kb_windows <- read_delim("genmap_k150_e2_20kb_windows.txt","\t", escape_double = FALSE, trim_ws = TRUE)
genmap_k150_e2_20kb_windows$BIN_START= genmap_k150_e2_20kb_windows$BIN_START+1

## RR carried over from Bt. linkage map
RR <- read_delim("RR_20Kb_windows.txt","\t", escape_double = FALSE, trim_ws = TRUE)


# Exon content
exons <- read_delim("Bsyl_exon_content_20kb.txt","\t", escape_double=FALSE, trim_ws=TRUE)

# Gene counts in windows
gene_counts <- read_delim("Bsyl_gene_counts_20kb.txt","\t", escape_double=FALSE, trim_ws=TRUE)


repeat_content <- read_delim("repeat_content_20kb_windows.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # THIS HAS DUPLICATED ENTRIES - NEEDS
#FIXING. FOR NOW, REMOVE DUPLICATES IN EXCEL BEFORE IMPORTING
repeat_content$BIN_START= repeat_content$BIN_START+1

## Import positions of 15bp abundant repeat - is it associated with centromeres?
STR_15bp_positions <- read_delim("abundant_15bp_repeat_genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# and regions where the repeat appears in clusters from Vale's analysis
STR_15bp_cluster_positions <- read_delim("sat15_genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)



# number of SNPs per window
Bs_Bi_snps_per_20kb_window <- read_delim("Bs_Bi_snps_per_20kb_window.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Bs_Bi_snps_per_20kb_window$BIN_START= Bs_Bi_snps_per_20kb_window$BIN_START+1
names(Bs_Bi_snps_per_20kb_window)[names(Bs_Bi_snps_per_20kb_window) == "number_of_snps"] <- "Bs_Bi_number_of_snps"

Bb_Bv_snps_per_20kb_window <- read_delim("Bb_Bv_snps_per_20kb_window.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Bb_Bv_snps_per_20kb_window$BIN_START= Bb_Bv_snps_per_20kb_window$BIN_START+1
names(Bb_Bv_snps_per_20kb_window)[names(Bb_Bv_snps_per_20kb_window) == "number_of_snps"] <- "Bb_Bv_number_of_snps"

# Pi
Bs_niwot_windowed_pi <- read_delim("Bs_niw_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(Bs_niwot_windowed_pi)[names(Bs_niwot_windowed_pi) == "PI"] <- "Bs_niw_PI"
Bs_niwot_windowed_pi$GEN_MID = Bs_niwot_windowed_pi$GEN_MID - 0.5

Bs_quail_windowed_pi <- read_delim("Bs_qua_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(Bs_quail_windowed_pi)[names(Bs_quail_windowed_pi) == "PI"] <- "Bs_qua_PI"
Bs_quail_windowed_pi$GEN_MID = Bs_quail_windowed_pi$GEN_MID - 0.5

sylv_windowed_pi <- read_delim("sylv_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(sylv_windowed_pi)[names(sylv_windowed_pi) == "PI"] <- "sylv_PI"
sylv_windowed_pi$GEN_MID = sylv_windowed_pi$GEN_MID - 0.5
inco_windowed_pi <- read_delim("inco_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(inco_windowed_pi)[names(inco_windowed_pi) == "PI"] <- "inco_PI"
inco_windowed_pi$GEN_MID = inco_windowed_pi$GEN_MID - 0.5
mela_windowed_pi <- read_delim("mela_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(mela_windowed_pi)[names(mela_windowed_pi) == "PI"] <- "mela_PI"
mela_windowed_pi$GEN_MID = mela_windowed_pi$GEN_MID - 0.5

bifa_windowed_pi <- read_delim("bifa_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(bifa_windowed_pi)[names(bifa_windowed_pi) == "PI"] <- "bifa_PI"
bifa_windowed_pi$GEN_MID = bifa_windowed_pi$GEN_MID - 0.5

vanco_windowed_pi <- read_delim("vanco_20kb.windowed.pi.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(vanco_windowed_pi)[names(vanco_windowed_pi) == "PI"] <- "vanco_PI"
vanco_windowed_pi$GEN_MID = vanco_windowed_pi$GEN_MID - 0.5

# Fst
# Within species
Bs_niwot_quail_fst_window <- read_delim("bsyl_niw_qua_phased.windowed.weir.fst", "\t", escape_double = FALSE, trim_ws = TRUE)
names(Bs_niwot_quail_fst_window)[names(Bs_niwot_quail_fst_window) == "WEIGHTED_FST"] <- "Bs.N.Q.fst"

# Between species
sylv_inco_fst_window <- read_delim("sylv_inco_20kb.windowed.weir.fst.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(sylv_inco_fst_window)[names(sylv_inco_fst_window) == "WEIGHTED_FST"] <- "Bs.Bi.fst"
sylv_inco_fst_window$GEN_MID = sylv_inco_fst_window$GEN_MID - 0.5
sylv_mela_fst_window <- read_delim("sylv_mela_20kb.windowed.weir.fst.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(sylv_mela_fst_window)[names(sylv_mela_fst_window) == "WEIGHTED_FST"] <- "Bs.Bm.fst"
sylv_mela_fst_window$GEN_MID = sylv_mela_fst_window$GEN_MID - 0.5
inco_mela_fst_window <- read_delim("inco_mela_20kb.windowed.weir.fst.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(inco_mela_fst_window)[names(inco_mela_fst_window) == "WEIGHTED_FST"] <- "Bi.Bm.fst"
inco_mela_fst_window$GEN_MID = inco_mela_fst_window$GEN_MID - 0.5

bifa_vanco_fst_window <- read_delim("bifa_vanco_20kb.windowed.weir.fst.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(bifa_vanco_fst_window)[names(bifa_vanco_fst_window) == "WEIGHTED_FST"] <- "Bb.Bv.fst"
bifa_vanco_fst_window$GEN_MID = bifa_vanco_fst_window$GEN_MID - 0.5

bifa_mela_fst_window <- read_delim("bifa_mela_20kb.windowed.weir.fst.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(bifa_mela_fst_window)[names(bifa_mela_fst_window) == "WEIGHTED_FST"] <- "Bb.Bm.fst"
bifa_mela_fst_window$GEN_MID = bifa_mela_fst_window$GEN_MID - 0.5

vanco_mela_fst_window <- read_delim("vanco_mela_20kb.windowed.weir.fst.genome_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)
names(vanco_mela_fst_window)[names(vanco_mela_fst_window) == "WEIGHTED_FST"] <- "Bv.Bm.fst"
vanco_mela_fst_window$GEN_MID = vanco_mela_fst_window$GEN_MID - 0.5


#TajD
#sylv_TD <- read_delim("../window_analysis/tajimas_D/sylv_20kb.Tajima.D.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#names(sylv_TD)[names(sylv_TD) == "TD"] <- "sylv_TD"
#inco_TD <- read_delim("../window_analysis/tajimas_D/inco_20kb.Tajima.D.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#names(inco_TD)[names(inco_TD) == "TD"] <- "inco_TD"
#mela_TD <- read_delim("../window_analysis/tajimas_D/mela_20kb.Tajima.D.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#names(mela_TD)[names(mela_TD) == "TD"] <- "mela_TD"
#bifa_TD <- read_delim("../window_analysis/tajimas_D/bifa_20kb.Tajima.D.genome_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#names(bifa_TD)[names(bifa_TD) == "TD"] <- "bifa_TD"

#Dxy
Bs_NQ_dxy <- read_delim("sylv_niwot_quail_dxy_20kb_windows_COMPLETE.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(Bs_NQ_dxy)[names(Bs_NQ_dxy) == "dxy"] <- "Bs_NQ_dxy"
Bs_NQ_dxy$BIN_START= Bs_NQ_dxy$BIN_START+1

sylv_inco_dxy <- read_delim("sylv_inco_dxy_20kb_windows_complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(sylv_inco_dxy)[names(sylv_inco_dxy) == "dxy"] <- "sylv_inco_dxy"
sylv_inco_dxy$BIN_START= sylv_inco_dxy$BIN_START+1


bifa_vanco_dxy <- read_delim("bifa_vanco_dxy_20kb_windows.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(bifa_vanco_dxy)[names(bifa_vanco_dxy) == "dxy"] <- "bifa_vanco_dxy"
bifa_vanco_dxy$BIN_START= bifa_vanco_dxy$BIN_START+1


## Merge all the data together so we have one data frame with coordinates, GC content, FST, Pi, Tajima's D, Dxy, 
# Mapping depth, and RR, and then calculate PBS

all_data_20kb <- merge(x = gc_20kb, y = DoC_Bsyl[ , c("CHROM","BIN_START", "DoC_Bsyl", "DoC_C2")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = exons[ , c("CHROM","BIN_START", "prop_exon")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = gene_counts[ , c("CHROM","BIN_START", "number_of_overlapping_genes")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = RR[ , c("Bt_LG","GEN_START", "RR")], by = c("Bt_LG","GEN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = Bs_niwot_windowed_pi[ , c("CHROM","BIN_START", "Bs_niw_PI")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = Bs_quail_windowed_pi[ , c("CHROM","BIN_START", "Bs_qua_PI")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = sylv_windowed_pi[ , c("CHROM","BIN_START", "sylv_PI")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = inco_windowed_pi[ , c("CHROM","BIN_START", "inco_PI")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = bifa_windowed_pi[ , c("CHROM","BIN_START", "bifa_PI")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = mela_windowed_pi[ , c("CHROM","BIN_START", "mela_PI")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = vanco_windowed_pi[ , c("CHROM","BIN_START", "vanco_PI")], by = c("CHROM","BIN_START"))


all_data_20kb <- merge(x = all_data_20kb, y = Bs_niwot_quail_fst_window[ , c("CHROM","BIN_START", "Bs.N.Q.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = sylv_inco_fst_window[ , c("Bt_LG","CHROM","BIN_START","GEN_START","GEN_END","GEN_MID", "Bs.Bi.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = sylv_mela_fst_window[ , c("CHROM","BIN_START", "Bs.Bm.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = inco_mela_fst_window[ , c("CHROM","BIN_START", "Bi.Bm.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = bifa_vanco_fst_window[ , c("CHROM","BIN_START", "Bb.Bv.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = bifa_mela_fst_window[ , c("CHROM","BIN_START", "Bb.Bm.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = vanco_mela_fst_window[ , c("CHROM","BIN_START", "Bv.Bm.fst")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = Bs_NQ_dxy[ , c("CHROM","BIN_START", "Bs_NQ_dxy")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = sylv_inco_dxy[ , c("CHROM","BIN_START", "sylv_inco_dxy")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = bifa_vanco_dxy[ , c("CHROM","BIN_START", "bifa_vanco_dxy")], by = c("CHROM","BIN_START"))
all_data_20kb$sylv_inco_dxy = all_data_20kb$sylv_inco_dxy/100 # I had *100 to make this a percentage in the perl script but it shouldn't be
all_data_20kb$Bs_NQ_dxy = all_data_20kb$Bs_NQ_dxy/100 # I had *100 to make this a percentage in the perl script but it shouldn't be
# LDhat results
all_data_20kb <- merge(x = all_data_20kb, y = LDhat_results[ , c("CHROM","BIN_START", "rho_per_kb")], by = c("CHROM","BIN_START"))

#drops <- c("mappability")
#all_data_20kb <- all_data_20kb[ , !(names(all_data_20kb) %in% drops)] # To drop the mappability column when I want to update with a new one
all_data_20kb <- merge(x = all_data_20kb, y = genmap_k150_e2_20kb_windows[ , c("CHROM","BIN_START", "mappability")], by = c("CHROM","BIN_START"))
all_data_20kb <- merge(x = all_data_20kb, y = repeat_content[ , c("CHROM","BIN_START", "prop_repeat")], by = c("CHROM","BIN_START")) # For some reason this is resulting in
# double entries for some windows....

# Number of SNPs per window for each comparison
#all_data_20kb <- merge(x = all_data_20kb, y = Bs_Bi_snps_per_20kb_window[ , c("CHROM","BIN_START", "Bs_Bi_number_of_snps")], by = c("CHROM","BIN_START")) 
#all_data_20kb <- merge(x = all_data_20kb, y = Bb_Bv_snps_per_20kb_window[ , c("CHROM","BIN_START", "Bb_Bv_number_of_snps")], by = c("CHROM","BIN_START")) 

## Calculate PBS from FST values using the equations from Yi et al 2010, Science

# Transform FST value:
# T = -log(1-FST)

all_data_20kb$Bs.Bi.T <- -log(1-(all_data_20kb$Bs.Bi.fst-0.01))
all_data_20kb$Bs.Bm.T <- -log(1-(all_data_20kb$Bs.Bm.fst-0.01))
all_data_20kb$Bi.Bm.T <- -log(1-(all_data_20kb$Bi.Bm.fst-0.01))

all_data_20kb$Bb.Bv.T <- -log(1-(all_data_20kb$Bb.Bv.fst-0.01))
all_data_20kb$Bb.Bm.T <- -log(1-(all_data_20kb$Bb.Bm.fst-0.01))
all_data_20kb$Bv.Bm.T <- -log(1-(all_data_20kb$Bv.Bm.fst-0.01))

# Calculate PBS for each species
all_data_20kb$B.s.PBS <- (all_data_20kb$Bs.Bi.T + all_data_20kb$Bs.Bm.T - all_data_20kb$Bi.Bm.T)/2
all_data_20kb$B.i.PBS <- (all_data_20kb$Bs.Bi.T + all_data_20kb$Bi.Bm.T - all_data_20kb$Bs.Bm.T)/2
all_data_20kb$B.b.PBS <- (all_data_20kb$Bb.Bv.T + all_data_20kb$Bb.Bm.T - all_data_20kb$Bv.Bm.T)/2
all_data_20kb$B.v.PBS <- (all_data_20kb$Bb.Bv.T + all_data_20kb$Bv.Bm.T - all_data_20kb$Bb.Bm.T)/2

# Remove data points for DoC < 20

all_data_20kb <-all_data_20kb[all_data_20kb$DoC_Bsyl < 20,]
all_data_20kb <-all_data_20kb[all_data_20kb$DoC_C2 < 20,]


# all windows mean values
all_data_gc_mean <- mean(all_data_20kb$GC_prop)
all_data_cov_mean <- mean(all_data_20kb$DoC_Bsyl)

quantile(all_data_20kb$B.s.PBS, c(.9))
quantile(all_data_20kb$B.i.PBS, c(.9))
quantile(all_data_20kb$B.b.PBS, c(.9))
quantile(all_data_20kb$B.v.PBS, c(.9))

## Convert FSTs to Z-scores
# Bs Niwot-Quail

Bs_Niwot_Quail_mean_fst = mean(all_data_20kb$Bs.N.Q.fst)
Bs_Niwot_Quail_median_fst =median(all_data_20kb$Bs.N.Q.fst)
Bs_Niwot_Quail_sd_fst = sd(all_data_20kb$Bs.N.Q.fst)
all_data_20kb$Bs.N.Q.zfst = (all_data_20kb$Bs.N.Q.fst-Bs_Niwot_Quail_median_fst)/Bs_Niwot_Quail_sd_fst
Bs_Niwot_Quail_sd_Zfst = sd(all_data_20kb$Bs.N.Q.zfst)
sum(all_data_20kb$Bs.N.Q.zfst >= 2)
median(all_data_20kb$Bs.N.Q.zfst)
mean(all_data_20kb$Bs.N.Q.zfst)

# % of genome that is 'highly divergent'
Bs_N_Q_genome_divergent_perc <- (486*20000)/252081862*100
# 3.86%

#Bs-Bi
Bs_Bi_mean_fst = mean(all_data_20kb$Bs.Bi.fst)
Bs_Bi_median_fst =median(all_data_20kb$Bs.Bi.fst)
Bs_Bi_sd_fst = sd(all_data_20kb$Bs.Bi.fst)
#all_data_20kb$Bs.Bi.zfst = (all_data_20kb$Bs.Bi.fst-Bs_Bi_mean_fst)/Bs_Bi_sd_fst
all_data_20kb$Bs.Bi.zfst = (all_data_20kb$Bs.Bi.fst-Bs_Bi_median_fst)/Bs_Bi_sd_fst # Use median instead of mean due to bimodal distribution
Bs_Bi_sd_Zfst = sd(all_data_20kb$Bs.Bi.zfst)
sum(all_data_20kb$Bs.Bi.zfst >= 2)
median(all_data_20kb$Bs.Bi.zfst)
mean(all_data_20kb$Bs.Bi.zfst)

# % of genome that is 'highly divergent'
BsBi_genome_divergent_perc <- (1758*20000)/252081862*100
# 13.95%

#Bb-Bv
Bb_Bv_mean_fst = mean(all_data_20kb$Bb.Bv.fst)
Bb_Bv_median_fst = median(all_data_20kb$Bb.Bv.fst)
Bb_Bv_sd_fst = sd(all_data_20kb$Bb.Bv.fst)
#all_data_20kb$Bb.Bv.zfst = (all_data_20kb$Bb.Bv.fst-Bb_Bv_mean_fst)/Bb_Bv_sd_fst
all_data_20kb$Bb.Bv.zfst = (all_data_20kb$Bb.Bv.fst-Bb_Bv_median_fst)/Bb_Bv_sd_fst
Bb_Bv_sd_Zfst = sd(all_data_20kb$Bb.Bv.zfst)
sum(all_data_20kb$Bb.Bv.zfst >= 2)
median(all_data_20kb$Bb.Bv.zfst)
mean(all_data_20kb$Bb.Bv.zfst)

# % of genome that is 'highly divergent'
BbBv_genome_divergent_perc <- (1060*20000)/252081862*100
#8.41%

# Plot histograms of Fst

# set up x-axis plot
zfst_x_axis_plot <- ggplot(data=all_data_20kb, aes(x=Bs.Bi.fst)) +
  labs(x=expression(bolditalic(F["ST"]))) +
  scale_x_continuous(limits=c(-0.0075,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.001)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=14,face="bold"),axis.title.y=element_blank(), 
        axis.text.y=element_blank(),legend.position = "none", axis.ticks.y = element_blank())

zfst_x_axis_plot


## Set values <0 to be =0
NQ_below_zero = all_data_20kb$Bs.N.Q.fst<0
all_data_20kb$Bs.N.Q.fst[NQ_below_zero] = 0



Bs_Niwot_Quail_fst_hist <- ggplot(data=all_data_20kb, aes(x=Bs.N.Q.fst)) +
  geom_histogram(binwidth = 0.005,aes(y=..count../sum(..count..)),position = 'identity', fill="green4") +
  #geom_density(aes(y=..count../sum(..count..)), alpha=0.3) +
  annotate("text", x = 0.7, y = 0.16, size = 4, label = expression(paste(bolditalic("B. sylvicola "), bold("Niwot ridge - Quail Mountain")))) +
  scale_x_continuous(limits=c(-0.0075,1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_manual(values=c("red4","goldenrod3"),labels = c("B. sylvicola - B. incognita (S)", "B. bifarius - B. vancouverensis (A)")) +
  labs(x=expression(F["ST"]), y=expression("Frequency")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = c(0.7,0.8),
        axis.title.x=element_blank(),axis.text.x=element_blank())
Bs_Niwot_Quail_fst_hist


BsBi_fst_hist <- ggplot(data=all_data_20kb, aes(x=Bs.Bi.fst)) +
  geom_histogram(binwidth = 0.005, aes(y=..count../sum(..count..)),position = 'identity', fill="magenta") +
  annotate("text", x = 0.7, y = 0.025, size = 4, label = expression(paste(bolditalic("B. sylvicola - B. incognita ")))) +
  #geom_density(aes(y=..count../sum(..count..)), alpha=0.3) +
  #geom_vline(xintercept=0.93) +
  scale_x_continuous(limits=c(-0.0075,1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_manual(values=c("red4","goldenrod3"),labels = c("B. sylvicola - B. incognita (S)", "B. bifarius - B. vancouverensis (A)")) +
  labs(x=expression(F["ST"]), y=expression("Frequency")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = c(0.7,0.8),
        axis.title.x=element_blank(),axis.text.x=element_blank())
BsBi_fst_hist

BbBv_fst_hist <- ggplot(data=all_data_20kb, aes(x=Bb.Bv.fst)) +
  geom_histogram(binwidth = 0.005, aes(y=..count../sum(..count..)),position = 'identity', fill="blue") +
  annotate("text", x = 0.7, y = 0.05, size = 4, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis")))) +
  #geom_vline(xintercept=0.17) +
  #geom_density(aes(y=..count../sum(..count..)), alpha=0.3) +
  scale_x_continuous(limits=c(-0.0075,1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_manual(values=c("red4","goldenrod3"),labels = c("B. sylvicola - B. incognita (S)", "B. bifarius - B. vancouverensis (A)")) +
  labs(x=expression(F["ST"]), y=expression("Frequency")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = c(0.7,0.8),
        axis.title.x=element_blank(),axis.text.x=element_blank())
        
BbBv_fst_hist


ggarrange(Bs_Niwot_Quail_fst_hist,BsBi_fst_hist,BbBv_fst_hist,zfst_x_axis_plot, heights = c(1,1,1,0.2),ncol = 1, nrow = 4, align = "v")


# Define highly divergent windows based on FST, where windows with ZFST >= 2 are 'highly divergent'

Bs_N_Q_high_FST <- all_data_20kb$Bs.N.Q.zfst >= 2
Bs_Bi_high_FST <- all_data_20kb$Bs.Bi.zfst >= 2
Bb_Bv_high_FST <- all_data_20kb$Bb.Bv.zfst >= 2

all_data_20kb$Bs_N_Q_divergent_windows <- "Background"
all_data_20kb$Bs_N_Q_divergent_windows[Bs_N_Q_high_FST] <- "Highly divergent"

all_data_20kb$Bs_Bi_divergent_windows <- "Background"
all_data_20kb$Bs_Bi_divergent_windows[Bs_Bi_high_FST] <- "Highly divergent"

all_data_20kb$Bb_Bv_divergent_windows <- "Background"
all_data_20kb$Bb_Bv_divergent_windows[Bb_Bv_high_FST] <- "Highly divergent"

## Need to merge neighbouring divergent windows into single divergent blocks using perl scripts (see "merging_divergent_windows" directory)

# Output coordinates for all windows in order to assign IoD classification to each window using perl scripts
window_coords <- data.frame(all_data_20kb$Bt_LG, all_data_20kb$GEN_START)
write.table(window_coords, file="window_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Output Bs-Bi outlier windows
BsBi_fst_outliers <- subset(all_data_20kb, Bs.Bi.zfst >= 2)
BsBi_fst_outliers_coords <- data.frame(BsBi_fst_outliers$Bt_LG, BsBi_fst_outliers$GEN_START, BsBi_fst_outliers$CHROM,BsBi_fst_outliers$BIN_START)
write.table(BsBi_fst_outliers_coords, file="BsBi_fst_outliers_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Then import in IoD assignments
BsBi_window_coords_with_IoD_assignment <- read_delim("BsBi_window_coords_with_IoD_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
drops <- c("BsBi_IoD")
all_data_20kb <- all_data_20kb[ , !(names(all_data_20kb) %in% drops)] # To drop the column when I want to update with a new one
all_data_20kb <- merge(x = all_data_20kb, y = BsBi_window_coords_with_IoD_assignment[ , c("Bt_LG","GEN_START", "BsBi_IoD")], by = c("Bt_LG","GEN_START"))

## Do the same for Bb - Bv comparison

# Output Bb-Bv outlier windows
BbBv_fst_outliers <- subset(all_data_20kb, Bb.Bv.zfst >= 2)
BbBv_fst_outliers_coords <- data.frame(BbBv_fst_outliers$Bt_LG, BbBv_fst_outliers$GEN_START, BbBv_fst_outliers$CHROM,BbBv_fst_outliers$BIN_START)
write.table(BbBv_fst_outliers_coords, file="BbBv_fst_outliers_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Then import in IoD assignments
BbBv_window_coords_with_IoD_assignment <- read_delim("BbBv_window_coords_with_IoD_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
drops <- c("BbBv_IoD")
all_data_20kb <- all_data_20kb[ , !(names(all_data_20kb) %in% drops)] # To drop the column when I want to update with a new one
all_data_20kb <- merge(x = all_data_20kb, y = BbBv_window_coords_with_IoD_assignment[ , c("Bt_LG","GEN_START", "BbBv_IoD")], by = c("Bt_LG","GEN_START"))

## And for Bs Niwot-Quail comparison

# Output Bs Niwot-Quail outlier windows
Bs.N.Q_fst_outliers <- subset(all_data_20kb, Bs.N.Q.zfst >= 2)
Bs.N.Q_fst_outliers_coords <- data.frame(Bs.N.Q_fst_outliers$Bt_LG, Bs.N.Q_fst_outliers$GEN_START, Bs.N.Q_fst_outliers$CHROM,Bs.N.Q_fst_outliers$BIN_START)
write.table(Bs.N.Q_fst_outliers_coords, file="Bs.N.Q_fst_outliers_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Then import in IoD assignments
Bs.N.Q.window_coords_with_IoD_assignment <- read_delim("Bs.N.Q.window_coords_with_IoD_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
drops <- c("Bs.N.Q_IoD")
all_data_20kb <- all_data_20kb[ , !(names(all_data_20kb) %in% drops)] # To drop the column when I want to update with a new one
all_data_20kb <- merge(x = all_data_20kb, y = Bs.N.Q.window_coords_with_IoD_assignment[ , c("Bt_LG","GEN_START", "Bs.N.Q_IoD")], by = c("Bt_LG","GEN_START"))

## Then look at overlap
all_data_20kb$shared_divergent_window = "FALSE"

for (row in 1:nrow(all_data_20kb)) {
  
  BsNQ  <- all_data_20kb[row, "Bs.N.Q_IoD"]
  BsBi <- all_data_20kb[row, "BsBi_merged_IoD"]
  BbBv  <- all_data_20kb[row, "BbBv_merged_IoD"]
  
  if((BsBi == "IoD") && (BbBv == "IoD")) {
    all_data_20kb[row, "shared_divergent_window"] = "TRUE"
  }
}

sum(all_data_20kb$shared_divergent_window == "TRUE")
sum(all_data_20kb$BbBv_merged_IoD=="IoD")



# 522
# Number of 20Kb windows in IoDs in BbBv comparison:
sum(all_data_20kb$BbBv_IoD == "IoD") # 1,333
perc_BbBv_IoDs_in_BsBi_IoDs <- 522/1333*100
# 39.16%

perc_BsBi_IoDs_overlapping_with_BbBv_IoDs <- 522/1894*100
#27.56%

sum(all_data_20kb$BsBi_IoD == "IoD")
# Check significance of overlap using answer from here: https://stats.stackexchange.com/questions/267/how-do-i-calculate-if-the-degree-of-overlap-between-two-lists-is-significant
# and here https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
n_A = 1894;n_B = 1333; n_C = 11923; n_A_B = 522
phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE) 

## Bs.N.Q overlap
all_data_20kb$shared_divergent_window = "FALSE"

for (row in 1:nrow(all_data_20kb)) {
  
  BsBi <- all_data_20kb[row, "BsBi_merged_IoD"]
  BsNQ  <- all_data_20kb[row, "Bs.N.Q_merged_IoD_over_100kb"]
  
  if((BsBi == "IoD") && (BsNQ == "IoD")) {
    all_data_20kb[row, "shared_divergent_window"] = "TRUE"
  }
}

sum(all_data_20kb$shared_divergent_window == "TRUE")
sum(all_data_20kb$Bs.N.Q_merged_IoD_over_100kb == "IoD")

# 388
# Number of 20Kb windows in IoDs in BsNQ comparison:
sum(all_data_20kb$Bs.N.Q_IoD == "IoD") # 590
perc_BsNQ_IoDs_in_BsBi_IoDs <- 388/590*100
# 65.76%

# Number of 20Kb windows in IoDs in BsBi comparison:
sum(all_data_20kb$BsBi_IoD == "IoD") # 1894
perc_BsBi_IoDs_overlapping_with_BsNQ_IoDs <- 388/1894*100
# 20.49%

# Check significance of overlap using answer from here: https://stats.stackexchange.com/questions/267/how-do-i-calculate-if-the-degree-of-overlap-between-two-lists-is-significant
n_A = 1894;n_B = 590; n_C = 11923; n_A_B = 388
phyper(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE) 

# Calculate overlap between all comparisons for Venn Diagram

all_data_20kb$shared_divergent_window = "FALSE"

for (row in 1:nrow(all_data_20kb)) {
  
  BsBi <- all_data_20kb[row, "BsBi_IoD"]
  BbBv <- all_data_20kb[row, "BbBv_IoD"]
  BsNQ  <- all_data_20kb[row, "Bs.N.Q_IoD"]
  
  if((BsBi == "background") && (BbBv == "background") && (BsNQ == "IoD")) {
    all_data_20kb[row, "shared_divergent_window"] = "TRUE"
  }
}

sum(all_data_20kb$shared_divergent_window == "TRUE")

## Make Venn Diagram of number of outlier windows

library(VennDiagram)
overrideTriple=T
grid.newpage();
venn.plot <- draw.triple.venn(
  area1 = 1894,
  area2 = 1333,
  area3 = 590,
  n12 = 522,
  n23 = 137,
  n13 = 388,
  n123 = 87,
  euler.d = TRUE, scaled = TRUE,
  category = c("Sympatric", "Allopatric", "Within-species"),
  fill = c("green4","magenta", "orange"),
  lty = "blank",
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cex = 1,
  cat.cex = 1.2,
  cat.col = c("green4","magenta", "orange"),print.mode = "raw"
);
grid.draw(venn.plot);


## As percentages of the genome
grid.newpage();
venn.plot <- draw.triple.venn(
  area1 = 15.9,
  area2 = 11.2,
  area3 = 4.9,
  n12 = 4.4,
  n23 = 1.1,
  n13 = 3.3,
  n123 = 0.7,
  euler.d = TRUE, scaled = TRUE,
  category = c("Sympatric", "Allopatric", "Within-species"),
  fill = c("green4","magenta", "orange"),
  lty = "blank",
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cex = 1,
  cat.cex = 1.2,
  cat.col = c("green4","magenta", "orange"),print.mode = "raw"
);
grid.draw(venn.plot);

## As sequence length of IoDs
grid.newpage();
venn.plot <- draw.triple.venn(
  area1 = 37.9,
  area2 = 26.7,
  area3 = 11.8,
  n12 = 10.4,
  n23 = 2.7,
  n13 = 7.8,
  n123 = 1.7,
  euler.d = TRUE, scaled = TRUE,
  category = c("Sympatric", "Allopatric", "Within-species"),
  fill = c("green4","magenta", "orange"),
  lty = "blank",
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cex = 1,
  cat.cex = 1.2,
  cat.col = c("green4","magenta", "orange"),print.mode = "raw"
);
grid.draw(venn.plot);


## IoDs of at least 100 Kb in length merged if within 1Mb on a chromosome
# Bs Niw-Qua
#First without merging windows within 1Mb
BsNQ_window_coords_with_100kb_IoD_assignment <- read_delim("Bs.N.Q.window_coords_with_IoD_over_100kb_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
all_data_20kb <- merge(x = all_data_20kb, y = BsNQ_window_coords_with_100kb_IoD_assignment[ , c("Bt_LG","GEN_START", "Bs.N.Q_IoD_over_100kb")], by = c("Bt_LG","GEN_START"))
#After merging windows within 1Mb
BsNQ_window_coords_with_100kb_IoD_merged_assignment <- read_delim("Bs.N.Q.window_coords_with_IoD_over_100kb_merged_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
all_data_20kb <- merge(x = all_data_20kb, y = BsNQ_window_coords_with_100kb_IoD_merged_assignment[ , c("Bt_LG","GEN_START", "Bs.N.Q_merged_IoD_over_100kb")], by = c("Bt_LG","GEN_START"))

# Bs - Bi
BsBi_window_coords_with_100kb_IoD_assignment <- read_delim("window_coords_with_100kb_IoDs_merged_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
drops <- c("BsBi_100kb_IoD")
all_data_20kb <- all_data_20kb[ , !(names(all_data_20kb) %in% drops)] # To drop the column when I want to update with a new one
all_data_20kb <- merge(x = all_data_20kb, y = BsBi_window_coords_with_100kb_IoD_assignment[ , c("Bt_LG","GEN_START", "BsBi_merged_IoD")], by = c("Bt_LG","GEN_START"))

## And for Bb-Bv
#First without merging windows within 1Mb
BbBv_window_coords_with_100kb_IoD_assignment <- read_delim("BbBv_window_coords_with_IoD_over_100kb_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
all_data_20kb <- merge(x = all_data_20kb, y = BbBv_window_coords_with_100kb_IoD_assignment[ , c("Bt_LG","GEN_START", "BbBv_IoD_over_100kb")], by = c("Bt_LG","GEN_START"))

#After merging windows within 1Mb
BbBv_window_coords_with_100kb_IoD_merged_assignment <- read_delim("BbBv_window_coords_with_merged_IoD_assignment.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
all_data_20kb <- merge(x = all_data_20kb, y = BbBv_window_coords_with_100kb_IoD_merged_assignment[ , c("Bt_LG","GEN_START", "BbBv_merged_IoD")], by = c("Bt_LG","GEN_START"))


# Test significance of overlap using cooccur package
library(cooccur)

# Make data frame for each comparison

# B.s NQ - v BsBi
# Check the overlap between the sympatric and within-species largest IoDs
BsNQ_BsBi_cooccur_df <- data.frame(all_data_20kb$Bs.N.Q_merged_IoD_over_100kb, all_data_20kb$BsBi_merged_IoD,stringsAsFactors = FALSE)
names(BsNQ_BsBi_cooccur_df)[names(BsNQ_BsBi_cooccur_df) == "all_data_20kb.Bs.N.Q_merged_IoD_over_100kb"] <- "BsNQ"
names(BsNQ_BsBi_cooccur_df)[names(BsNQ_BsBi_cooccur_df) == "all_data_20kb.BsBi_merged_IoD"] <- "BsBi"

BsNQ_BsBi_cooccur_df <- na.omit(BsNQ_BsBi_cooccur_df)

# Set windows in and out of IoDs to be 1s and 0s resepctively (equivalent to 'presence' and 'absence')
BsNQ_BsBi_cooccur_df[BsNQ_BsBi_cooccur_df=="background"] <- 0
BsNQ_BsBi_cooccur_df[BsNQ_BsBi_cooccur_df=="IoD"] <- 1

#write.table(BsNQ_BsBi_cooccur_df , file="BsNQ_BsBi_cooccur_df.txt", sep="\t", row.names=TRUE, quote=FALSE)

sapply(BsNQ_BsBi_cooccur_df, mode) # To check what columns are set at (need to be numeric)
BsNQ_BsBi_cooccur_df <- transform(BsNQ_BsBi_cooccur_df, BsNQ = as.numeric(BsNQ), BsBi = as.numeric(BsBi))
BsNQ_BsBi_cooccur_df_t <- t(BsNQ_BsBi_cooccur_df) # To transpose the data 

BsNQ_BsBi_cooccur <- cooccur(mat=BsNQ_BsBi_cooccur_df_t, type = "spp_site", thresh = FALSE, spp_names = TRUE,
        true_rand_classifier = 0.1, prob = "hyper",
        site_mask = NULL, only_effects = FALSE,
        eff_standard = TRUE, eff_matrix = FALSE)
BsNQ_BsBi_cooccur

# BsBi - v BbBv
# Check the overlap between the sympatric largest IoDs and allopatric divergent windows (ZFST > 2)
BsBi_BbBv_cooccur_df <- data.frame(all_data_20kb$BsBi_merged_IoD,all_data_20kb$BbBv_IoD_over_100kb,stringsAsFactors = FALSE)
names(BsBi_BbBv_cooccur_df)[names(BsBi_BbBv_cooccur_df) == "all_data_20kb.BbBv_IoD_over_100kb"] <- "BbBv"
names(BsBi_BbBv_cooccur_df)[names(BsBi_BbBv_cooccur_df) == "all_data_20kb.BsBi_merged_IoD"] <- "BsBi"

# Set windows in and out of IoDs to be 1s and 0s resepctively (equivalent to 'presence' and 'absence')
BsBi_BbBv_cooccur_df[BsBi_BbBv_cooccur_df=="background"] <- 0
BsBi_BbBv_cooccur_df[BsBi_BbBv_cooccur_df=="IoD"] <- 1

sapply(BsBi_BbBv_cooccur_df, mode)
BsBi_BbBv_cooccur_df <- transform(BsBi_BbBv_cooccur_df, BbBv = as.numeric(BbBv), BsBi = as.numeric(BsBi))
BsBi_BbBv_cooccur_df_t <- t(BsBi_BbBv_cooccur_df)

BsBi_BbBv_cooccur <- cooccur(mat=BsBi_BbBv_cooccur_df_t, type = "spp_site", thresh = FALSE, spp_names = TRUE,
                             true_rand_classifier = 0.1, prob = "hyper",
                             site_mask = NULL, only_effects = FALSE,
                             eff_standard = TRUE, eff_matrix = FALSE)
BsBi_BbBv_cooccur

# B.s NQ - v BbBv
# Check the overlap between the sympatric and within-species largest IoDs
BsNQ_BbBv_cooccur_df <- data.frame(all_data_20kb$Bs.N.Q_merged_IoD_over_100kb, all_data_20kb$BbBv_IoD_over_100kb,stringsAsFactors = FALSE)
names(BsNQ_BbBv_cooccur_df)[names(BsNQ_BbBv_cooccur_df) == "all_data_20kb.Bs.N.Q_merged_IoD_over_100kb"] <- "BsNQ"
names(BsNQ_BbBv_cooccur_df)[names(BsNQ_BbBv_cooccur_df) == "all_data_20kb.BbBv_IoD_over_100kb"] <- "BbBv"

BsNQ_BbBv_cooccur_df <- na.omit(BsNQ_BbBv_cooccur_df)

# Set windows in and out of IoDs to be 1s and 0s resepctively (equivalent to 'presence' and 'absence')
BsNQ_BbBv_cooccur_df[BsNQ_BbBv_cooccur_df=="background"] <- 0
BsNQ_BbBv_cooccur_df[BsNQ_BbBv_cooccur_df=="IoD"] <- 1

#write.table(BsNQ_BbBv_cooccur_df , file="BsNQ_BbBv_cooccur_df.txt", sep="\t", row.names=TRUE, quote=FALSE)

sapply(BsNQ_BbBv_cooccur_df, mode) # To check what columns are set at (need to be numeric)
BsNQ_BbBv_cooccur_df <- transform(BsNQ_BbBv_cooccur_df, BsNQ = as.numeric(BsNQ), BbBv = as.numeric(BbBv))
BsNQ_BbBv_cooccur_df_t <- t(BsNQ_BbBv_cooccur_df) # To transpose the data 

BsNQ_BbBv_cooccur <- cooccur(mat=BsNQ_BbBv_cooccur_df_t, type = "spp_site", thresh = FALSE, spp_names = TRUE,
                             true_rand_classifier = 0.1, prob = "hyper",
                             site_mask = NULL, only_effects = FALSE,
                             eff_standard = TRUE, eff_matrix = FALSE)
BsNQ_BbBv_cooccur



mean_map<-mean(all_data_20kb$mappability)
#mean_map_IoD_both <- mean(all_data_20kb$mappability[IoD_Both_windows])

### WHOLE GENOME PLOTS ######
# Import coordinates for B terrestris chromsome borders
Bt_LG_coords <- read_delim("Bt_LG_coords.txt","\t", escape_double = FALSE, trim_ws = TRUE)

# Set up x-axis plot
x_axis_plot <- ggplot(data=Bt_LG_coords, aes(x=x,y=y)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(data=STR_15bp_cluster_positions, aes(x=start, y=0.05, colour=type), shape=15, alpha=0.5, size=4) +
  labs(x="Genome position (Mbp)") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0.03,0.07)) +
  scale_colour_manual(values=c("green", "purple")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"),axis.title.y=element_blank(), 
        axis.text.y=element_blank(),legend.position = "none", axis.ticks.y = element_blank())

x_axis_plot

## RR plot from LDHat
RR_genome_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = rho_per_kb)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  geom_line(size=0.2) +
  labs(x="Genome position (Mbp)", y="rho/kb") +
  scale_x_continuous(expand=c(0,0),breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),axis.title.x=element_blank(), legend.position = "none",axis.text.x=element_blank())
RR_genome_plot


## FST plots

Bs_N_Q_fst_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bs.N.Q.zfst)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  geom_point(size=0.2,aes(colour=Bs.N.Q_merged_IoD_over_100kb)) +
  geom_hline(yintercept=Bs_Niwot_Quail_mean_fst, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 18, size = 3, label = expression(paste(bolditalic("B. sylvicola "), bold("Niwot ridge - Quail Mountain")))) +
  labs(x="Genome position (Mbp)", y=expression(paste(bolditalic(ZF["ST"]))),colour="Islands of Divergence") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  #scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  #scale_colour_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                             "red","orange","red","orange","red","orange", "darkgrey")) +
  scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
                             "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),axis.title.x=element_blank(), legend.position = "none",axis.text.x=element_blank())
Bs_N_Q_fst_plot


Bs_Bi_fst_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bs.Bi.zfst)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=1, lty =2) +
  geom_point(size=0.2,aes(colour=BsBi_merged_IoD)) +
  geom_hline(yintercept=Bs_Bi_mean_fst, lty=2, colour="black") +
  #geom_vline(data = contig_lengths, aes(xintercept = End_pos), colour = "grey", alpha=0.5, lty =2) +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  #geom_point(data=STR_15bp_positions, aes(x=genome_pos, y=2.8), size=1) +
  annotate("text", x = 190000000, y = 3, size = 3, label = expression(paste(bolditalic("B. sylvicola - B. incognita ")))) +
  labs(x="Genome position (Mbp)", y=expression(paste(bolditalic(ZF["ST"]))),colour="Islands of Divergence") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  #scale_y_continuous(limits=c(0,1.2)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  #scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                          "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"),axis.title.x=element_blank(), legend.position = "none",axis.text.x=element_blank())
Bs_Bi_fst_plot


Bb_Bv_fst_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bb.Bv.zfst)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2,aes(colour=BbBv_IoD_over_100kb)) +
  geom_hline(yintercept=Bb_Bv_mean_fst, lty=2, colour="black") +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  annotate("text", x = 190000000, y = 5.5, size = 3, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis ")))) +
  labs(x="Genome position (Mbp)", y=expression(paste(bolditalic(ZF["ST"]))),colour="Islands of Divergence") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  #scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  #scale_colour_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                             "red","orange","red","orange","red","orange", "darkgrey")) +
  #scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                          "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"), axis.title.x=element_blank(),axis.text.x=element_blank(), legend.position = "none")
Bb_Bv_fst_plot

ggarrange(Bs_N_Q_fst_plot,Bs_Bi_fst_plot,Bb_Bv_fst_plot,RR_genome_plot,x_axis_plot, heights = c(1,1,1,0.8,0.5),ncol = 1, nrow = 5, align = "v")

ggarrange(Bs_N_Q_fst_plot,Bb_Bv_fst_plot,Bs_Bi_fst_plot,RR_genome_plot,x_axis_plot, heights = c(1,1,1,1,0.2),ncol = 1, nrow = 5, align = "v")
# PBS

#Bs

Bs_mean_pbs <- mean(all_data_20kb$B.s.PBS)

Bs_pbs_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = B.s.PBS)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  geom_hline(yintercept=Bs_mean_pbs, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 2.4, size = 3, label = expression(paste(bolditalic("B. sylvicola ")))) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  labs(x="Genome position (Mbp)", y=expression(paste(bold("PBS")))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.08,2.6)) +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  #scale_colour_manual(values=c("blue","red","red","red","darkgrey")) +
  #scale_colour_manual(values=c("orange","blue","purple","green","red","darkgrey")) +
  #scale_colour_manual(values=c("red","khaki","orange","blue","green", "pink", "darkgrey", "darkred")) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bs_pbs_plot

Bi_mean_pbs <- mean(all_data_20kb$B.i.PBS)

Bi_pbs_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = B.i.PBS)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  geom_hline(yintercept=Bi_mean_pbs, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 2.4, size = 3, label = expression(paste(bolditalic("B. incognita ")))) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  labs(x="Genome position (Mbp)", y=expression(paste(bold("PBS")))) +
  #geom_hline(yintercept=0.4651099, lty=2, colour="red") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.08,2.6)) +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  #scale_colour_manual(values=c("blue","red","red","red","darkgrey")) +
  #scale_colour_manual(values=c("orange","blue","purple","green","red","darkgrey")) +
  #scale_colour_manual(values=c("red","khaki","orange","blue","green", "pink", "darkgrey", "darkred")) +
  scale_colour_manual(values=c("darkgrey", "red")) +
  scale_fill_manual(values=c("goldenrod2","brown3","goldenrod2","brown3","goldenrod2","brown3","goldenrod2","brown3","goldenrod2","brown3","goldenrod2","brown3",
                               "goldenrod2","brown3","goldenrod2","brown3","goldenrod2","brown3", "darkgrey")) +
  
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bi_pbs_plot

Bb_mean_pbs <- mean(all_data_20kb$B.b.PBS)
Bb_pbs_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = B.b.PBS)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bb_Bv_divergent_windows)) +
  geom_hline(yintercept=Bb_mean_pbs, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 1.4, size = 3, label = expression(paste(bolditalic("B. bifarius ")))) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  labs(x="Genome position (Mbp)", y=expression(paste(bold("PBS")))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.08,1.5)) +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  #scale_colour_manual(values=c("blue","red","red","red","darkgrey")) +
  #scale_colour_manual(values=c("red","khaki","orange","blue","green", "pink", "darkgrey", "darkred")) +
  #scale_colour_manual(values=c("orange","blue","purple","green","red","darkgrey")) +
  scale_colour_manual(values=c("darkgrey", "red")) +
  scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
                             "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bb_pbs_plot


Bv_mean_pbs <- mean(all_data_20kb$B.v.PBS)
Bv_pbs_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = B.v.PBS)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_vline(data = contig_lengths, aes(xintercept = End_pos), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bb_Bv_divergent_windows)) +
  geom_hline(yintercept=Bv_mean_pbs, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 1.4, size = 3, label = expression(paste(bolditalic("B. vancouverensis ")))) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  labs(x="Genome position (Mbp)", y=expression(paste(bold("PBS")))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.08,1.5)) +
  #scale_colour_manual(values=c("red","khaki","orange","blue","green", "pink", "darkgrey", "darkred")) +
  #scale_colour_manual(values=c("orange","blue","purple","green","red","darkgrey")) +
  #scale_colour_manual(values=c("blue","red","red","red","darkgrey")) +
  scale_colour_manual(values=c("darkgrey","red")) +
  scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
                             "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")

Bv_pbs_plot



ggarrange(Bs_N_Q_fst_plot,Bs_pbs_plot,Bi_pbs_plot,Bb_pbs_plot,Bv_pbs_plot,x_axis_plot, heights = c(1,1,1,1,1,0.5),ncol = 1, nrow = 6, align = "v")



# Pi

Bs_niw_PI_mean <- mean(all_data_20kb$Bs_niw_PI)

Bs_Niwot_pi_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bs_niw_PI)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  geom_point(size=0.2,aes(colour=Bs_N_Q_divergent_windows)) +
  geom_hline(yintercept=Bs_niw_PI_mean, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 0.007, size = 3, label = expression(paste(bold("Niwot Ridge")))) +
  labs(x="Genome position (Mbp)", y="π",colour="Islands of Divergence") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  #scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  #scale_colour_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                             "red","orange","red","orange","red","orange", "darkgrey")) +
  scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
                             "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),axis.title.x=element_blank(), legend.position = "none",axis.text.x=element_blank())
Bs_Niwot_pi_plot


Bs_qua_PI_mean <- mean(all_data_20kb$Bs_qua_PI)

Bs_Quail_pi_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bs_qua_PI)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  geom_point(size=0.2,aes(colour=Bs_N_Q_divergent_windows)) +
  geom_hline(yintercept=Bs_qua_PI_mean, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 0.007, size = 3, label = expression(paste(bold("Quail Mountain")))) +
  labs(x="Genome position (Mbp)", y="π",colour="Islands of Divergence") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  #scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  #scale_colour_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                             "red","orange","red","orange","red","orange", "darkgrey")) +
  scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
                             "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),axis.title.x=element_blank(), legend.position = "none",axis.text.x=element_blank())
Bs_Quail_pi_plot


Bs_PI_mean <- mean(all_data_20kb$sylv_PI)

Bs_pi_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = sylv_PI)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  geom_hline(yintercept=Bs_PI_mean, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 0.0095, size = 3, label = expression(paste(bolditalic("B. sylvicola")))) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  labs(x="Genome position (Mbp)", y="π") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.01)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=14, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bs_pi_plot

Bi_PI_mean <- mean(all_data_20kb$inco_PI)

Bi_pi_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = inco_PI)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  geom_hline(yintercept=Bi_PI_mean, lty=2, colour="black") +
  annotate("text", x = 190000000, y = 0.0095, size = 3, label = expression(paste(bolditalic("B. incognita")))) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  labs(x="Genome position (Mbp)", y="π") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.01)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=14, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bi_pi_plot

Bb_PI_mean <- mean(all_data_20kb$bifa_PI)

Bb_pi_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = bifa_PI)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bb_Bv_divergent_windows)) +
  geom_hline(yintercept=Bb_PI_mean, lty=2, colour="black") +
  labs(x="Genome position (Mbp)", y="π") +
  annotate("text", x = 190000000, y = 0.0095, size = 3, label = expression(paste(bolditalic("B. bifarius")))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.01)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=14, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bb_pi_plot

Bv_PI_mean <- mean(all_data_20kb$vanco_PI)

Bv_pi_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = vanco_PI)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bb_Bv_divergent_windows)) +
  geom_hline(yintercept=Bv_PI_mean, lty=2, colour="black") +
  labs(x="Genome position (Mbp)", y="π") +
  annotate("text", x = 190000000, y = 0.0095, size = 3, label = expression(paste(bolditalic("B. vancouverensis")))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.01)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=14, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none")
Bv_pi_plot

ggarrange(Bs_pi_plot,Bi_pi_plot,Bb_pi_plot,Bv_pi_plot,x_axis_plot, heights = c(1,1,1,1,0.5),ncol = 1, nrow = 5, align = "v")

# dxy
Bs.N.Q_mean_dxy <- mean(all_data_20kb$ Bs_NQ_dxy)

Bs_N_Q_dxy_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bs_NQ_dxy)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line(aes(y=mappability), colour="blue", alpha=0.5) +
  #geom_polygon(data=Bt_LG_coords, aes(x=x,y=y,fill=LG_Bt)) +
  geom_point(size=0.2,aes(colour=Bs_N_Q_divergent_windows)) +
  geom_hline(yintercept=Bs.N.Q_mean_dxy, lty=2, colour="black") +
  #annotate("text", x = 190000000, y = 0.014, size = 3, label = expression(paste(bolditalic("B. sylvicola "), bold("Niwot ridge - Mt. Quail")))) +
  labs(x="Genome position (Mbp)", y=expression(paste(bolditalic(d["xy"]))),colour="Islands of Divergence") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  #scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  #scale_colour_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
  #                             "red","orange","red","orange","red","orange", "darkgrey")) +
  scale_fill_manual(values=c("red","orange","red","orange","red","orange","red","orange","red","orange","red","orange",
                             "red","orange","red","orange","red","orange")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),axis.title.x=element_blank(), legend.position = "none",axis.text.x=element_blank())
Bs_N_Q_dxy_plot


BsBi_mean_dxy <- mean(all_data_20kb$sylv_inco_dxy)

Bs_Bi_dxy_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = sylv_inco_dxy)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  #annotate("text", x = 190000000, y = 0.028, size = 3, label = expression(paste(bolditalic("B. sylvicola - B. incognita")))) +
  geom_hline(yintercept=0.0088116, lty=2, colour="black") +
  labs(x="Genome position (Mbp)", y=expression(paste(bold(d["xy "])))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.03)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank())
Bs_Bi_dxy_plot

BbBv_mean_dxy <- mean(all_data_20kb$bifa_vanco_dxy)

Bb_Bv_dxy_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = bifa_vanco_dxy)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_point(size=0.2, aes (colour=Bb_Bv_divergent_windows)) +
  #annotate("text", x = 190000000, y = 0.028, size = 3, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis")))) +
  geom_hline(yintercept=BbBv_mean_dxy, lty=2, colour="black") +
  labs(x="Genome position (Mbp)", y=expression(paste(bold(d["xy "])))) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.03)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"), axis.title.x=element_blank(), axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank())
Bb_Bv_dxy_plot


ggarrange(Bs_Bi_dxy_plot,Bb_Bv_dxy_plot, x_axis_plot, heights = c(1,1,0.2),ncol = 1, nrow = 3, align = "v")

ggarrange(Bs_N_Q_fst_plot,Bs_Niwot_pi_plot,Bs_Quail_pi_plot,Bs_N_Q_dxy_plot,x_axis_plot, heights = c(1,1,1,1,0.5),
          ncol = 1, nrow = 5, align = "v")

ggarrange(Bs_Bi_fst_plot, Bs_pbs_plot,Bi_pbs_plot,Bs_pi_plot,Bi_pi_plot,Bs_Bi_dxy_plot,x_axis_plot, heights = c(1,1,1,1,1,1,0.5),
          ncol = 1, nrow = 7, align = "v")

ggarrange(Bb_Bv_fst_plot,Bb_pbs_plot,Bv_pbs_plot,Bb_pi_plot,Bv_pi_plot,Bb_Bv_dxy_plot,x_axis_plot, heights = c(1,1,1,1,1,1,0.5),
          ncol = 1, nrow = 7, align = "v")


av.gc <- mean(all_data_20kb$GC_prop)

## GC plot
GC_genome_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = GC_prop)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes(colour=Bs_Bi_divergent_windows)) +
  geom_hline(aes(yintercept=av.gc), colour="blue",lty=2) +
  labs(x="Genome position (Mbp)", y="Proportion GC") +
  #scale_x_continuous(limits=c(186027913,187959241)) + # contig 54
  #scale_x_continuous(limits=c(123655592,126676346)) + # contig 28
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"), axis.title.x =element_blank(), axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank())
GC_genome_plot

## Mappability plot
mappability_genome_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = mappability)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  labs(x="Genome position (Mbp)", y="Mappability") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"), axis.title.x =element_blank(), axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank())
mappability_genome_plot

## Gene content plot
gene_content_genome_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = prop_gene)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes (colour=IoD)) +
  labs(x="Genome position (Mbp)", y="Proportion genic") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  #scale_colour_manual(values=c("purple","green","red","darkgrey")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.position = "none")
gene_content_genome_plot

## Depth of coverage plot
DoC_genome_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = DoC_Bsyl)) +
  geom_vline(data = contig_lengths, aes(xintercept = End_pos), colour = "grey", alpha=0.5, lty =2) +
  geom_line() +
  labs(x="Genome position (Mbp)", y="Depth of coverage") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  #scale_colour_manual(values=c("purple","green","red","darkgrey")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.position = "none")
DoC_genome_plot

## Repeat content plot
repeats_genome_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = prop_repeat)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  labs(x="Genome position (Mbp)", y="Proportion repeat content") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"), axis.title.x =element_blank(), axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank())
repeats_genome_plot

## Recombination rates plot
recom_rate_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = rho_per_kb)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes (colour=Bs_Bi_divergent_windows)) +
  labs(x="Genome position (Mbp)", y="Recombination rate (rho per Kbp)") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12, face="bold"),axis.title=element_text(size=12,face="bold"), axis.title.x =element_blank(), axis.text.x = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank())
recom_rate_plot 

ggarrange(GC_genome_plot,mappability_genome_plot,repeats_genome_plot,x_axis_plot, heights = c(1,1,1,0.2), ncol = 1, nrow = 4, align = "v")

## Bs-Bi number of SNPs per window
BsBi_snps_per_window_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bs_Bi_number_of_snps)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes (colour=IoD)) +
  labs(x="Genome position (Mbp)", y="Number of SNPs per window (B. sylvicola - B. incognita)") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.position = "none")
BsBi_snps_per_window_plot

## Bb-Bv number of SNPs per window
BbBv_snps_per_window_plot <- ggplot(data=all_data_20kb, aes(x = GEN_MID, y = Bb_Bv_number_of_snps)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  #geom_line() +
  geom_point(size=0.2, aes (colour=IoD)) +
  labs(x="Genome position (Mbp)", y="Number of SNPs per window (B. sylvicola - B. incognita)") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.position = "none")
BbBv_snps_per_window_plot

## BOXPLOTS #####

# First, convert to long format for easy plotting, separately for each comparison
# Bs Niwot-Quail
BsNQ_all_data_20kbp_long <- melt(all_data_20kb,idvar="Bs.N.Q_IoD")

# Define different metrics
Niwot.pi.lines <- BsNQ_all_data_20kbp_long$variable == "Bs_niw_PI"
BsNQ_all_data_20kbp_long$metric[Niwot.pi.lines] <- "pi"

Quail.pi.lines <- BsNQ_all_data_20kbp_long$variable == "Bs_qua_PI"
BsNQ_all_data_20kbp_long$metric[Quail.pi.lines] <- "pi"

NQ_dxy.lines <- BsNQ_all_data_20kbp_long$variable == "Bs_NQ_dxy"
BsNQ_all_data_20kbp_long$metric[NQ_dxy.lines] <- "dxy"


# Bs - Bi
BsBi_all_data_20kbp_long <- melt(all_data_20kb,idvar="BsBi_IoD")

# Define different metrics
Bs.pi.lines <- BsBi_all_data_20kbp_long$variable == "sylv_PI"
BsBi_all_data_20kbp_long$metric[Bs.pi.lines] <- "pi"

Bi.pi.lines <- BsBi_all_data_20kbp_long$variable == "inco_PI"
BsBi_all_data_20kbp_long$metric[Bi.pi.lines] <- "pi"

Bs.pbs.lines <- BsBi_all_data_20kbp_long$variable == "B.s.PBS"
BsBi_all_data_20kbp_long$metric[Bs.pbs.lines] <- "PBS"

Bi.pbs.lines <- BsBi_all_data_20kbp_long$variable == "B.i.PBS"
BsBi_all_data_20kbp_long$metric[Bi.pbs.lines] <- "PBS"

si_dxy.lines <- BsBi_all_data_20kbp_long$variable == "sylv_inco_dxy"
BsBi_all_data_20kbp_long$metric[si_dxy.lines] <- "dxy"

si_fst.lines <- BsBi_all_data_20kbp_long$variable == "Bs.Bi.fst"
BsBi_all_data_20kbp_long$metric[si_fst.lines] <- "fst"


# Bb - Bv
BbBv_all_data_20kbp_long <- melt(all_data_20kb,idvar="BbBv_IoD")

Bb.pi.lines <- BbBv_all_data_20kbp_long$variable == "bifa_PI"
BbBv_all_data_20kbp_long$metric[Bb.pi.lines] <- "pi"

Bv.pi.lines <- BbBv_all_data_20kbp_long$variable == "vanco_PI"
BbBv_all_data_20kbp_long$metric[Bv.pi.lines] <- "pi"

Bb.pbs.lines <- BbBv_all_data_20kbp_long$variable == "B.b.PBS"
BbBv_all_data_20kbp_long$metric[Bb.pbs.lines] <- "PBS"

Bv.pbs.lines <- BbBv_all_data_20kbp_long$variable == "B.v.PBS"
BbBv_all_data_20kbp_long$metric[Bv.pbs.lines] <- "PBS"

bv_dxy.lines <- BbBv_all_data_20kbp_long$variable == "bifa_vanco_dxy"
BbBv_all_data_20kbp_long$metric[bv_dxy.lines] <- "dxy"

bv_fst.lines <- BbBv_all_data_20kbp_long$variable == "Bb.Bv.fst"
BbBv_all_data_20kbp_long$metric[bv_fst.lines] <- "fst"


## boxplots of π

# Define "background" and "IoD" in each case first for stats
bsNQ_background <- all_data_20kb$Bs.N.Q_merged_IoD_over_100kb == "background"
bsNQ_IoD <- all_data_20kb$Bs.N.Q_merged_IoD_over_100kb == "IoD"

bsbi_background <- all_data_20kb$BsBi_merged_IoD == "background"
bsbi_IoD <- all_data_20kb$BsBi_merged_IoD == "IoD"

bbbv_background <- all_data_20kb$BbBv_merged_IoD == "background"
bbbv_IoD <- all_data_20kb$BbBv_merged_IoD == "IoD"

# Check significance of differences
wilcox.test(all_data_20kb$Bs_niw_PI[bsNQ_background], all_data_20kb$Bs_niw_PI[bsNQ_IoD])
wilcox.test(all_data_20kb$Bs_qua_PI[bsNQ_background], all_data_20kb$Bs_qua_PI[bsNQ_IoD])

# Calculate means
Bs_niw_PI_background_mean <- mean(all_data_20kb$Bs_niw_PI[bsNQ_background])
Bs_niw_PI_IoD_mean <- mean(all_data_20kb$Bs_niw_PI[bsNQ_IoD])
# Calculate % decrease
(Bs_niw_PI_background_mean - Bs_niw_PI_IoD_mean)/Bs_niw_PI_background_mean * 100

# Calculate means
Bs_qua_PI_background_mean <- mean(all_data_20kb$Bs_qua_PI[bsNQ_background])
Bs_qua_PI_IoD_mean <- mean(all_data_20kb$Bs_qua_PI[bsNQ_IoD])
# Calculate % decrease
(Bs_qua_PI_background_mean - Bs_qua_PI_IoD_mean)/Bs_qua_PI_background_mean * 100

BsNQ_pi_data <- subset(BsNQ_all_data_20kbp_long, metric=="pi")

BsNQ_pi_box_plot <- ggplot(data=BsNQ_pi_data, aes(y=value, x=variable, fill=Bs.N.Q_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #geom_signif(comparisons = list(c("background", "IoD")), y_position=0.007, map_signif_level=TRUE, textsize=6) +
  labs(x="Species", y="π", fill="") +
  scale_x_discrete(labels=c("Niwot Ridge", "Quail Mountain")) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"),axis.title.x=element_blank(), legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.008))

BsNQ_pi_box_plot


# Bs Bi Pi
# Check significance
wilcox.test(all_data_20kb$sylv_PI[bsbi_background], all_data_20kb$sylv_PI[bsbi_IoD])
wilcox.test(all_data_20kb$inco_PI[bsbi_background], all_data_20kb$inco_PI[bsbi_IoD])

BsBi_pi_data <- subset(BsBi_all_data_20kbp_long, metric=="pi")

# Calculate means
Bs_Pi_background_mean <- mean(all_data_20kb$sylv_PI[bsbi_background])
Bs_Pi_IoD_mean <- mean(all_data_20kb$sylv_PI[bsbi_IoD])
# Calculate % decrease
(Bs_Pi_background_mean - Bs_Pi_IoD_mean)/Bs_Pi_background_mean * 100

# Calculate means
Bi_Pi_background_mean <- mean(all_data_20kb$inco_PI[bsbi_background])
Bi_Pi_IoD_mean <- mean(all_data_20kb$inco_PI[bsbi_IoD])
# Calculate % decrease
(Bi_Pi_background_mean - Bi_Pi_IoD_mean)/Bi_Pi_background_mean * 100


BsBi_pi_box_plot <- ggplot(data=BsBi_pi_data, aes(y=value, x=variable, fill=BsBi_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #geom_signif(comparisons = list(c("background", "IoD")), y_position=0.007, map_signif_level=TRUE, textsize=6) +
  labs(x="Species", y="π") +
  scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text=element_text(size=12), axis.title=element_blank(),axis.text.y=element_blank(),legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.008))

BsBi_pi_box_plot

# Bb - Bv Pi
# Check significance
wilcox.test(all_data_20kb$bifa_PI[bbbv_background], all_data_20kb$bifa_PI[bbbv_IoD])
wilcox.test(all_data_20kb$vanco_PI[bbbv_background], all_data_20kb$vanco_PI[bbbv_IoD])

# Calculate means
Bb_Pi_background_mean <- mean(all_data_20kb$bifa_PI[bbbv_background])
Bb_Pi_IoD_mean <- mean(all_data_20kb$bifa_PI[bbbv_IoD])
# Calculate % decrease
(Bb_Pi_background_mean - Bb_Pi_IoD_mean)/Bb_Pi_background_mean * 100

# Calculate means
Bv_Pi_background_mean <- mean(all_data_20kb$vanco_PI[bbbv_background])
Bv_Pi_IoD_mean <- mean(all_data_20kb$vanco_PI[bbbv_IoD])
# Calculate % decrease
(Bv_Pi_background_mean - Bv_Pi_IoD_mean)/Bv_Pi_background_mean * 100


BbBv_pi_data <- subset(BbBv_all_data_20kbp_long, metric=="pi")

BbBv_pi_box_plot <- ggplot(data=BbBv_pi_data, aes(y=value, x=variable, fill=BbBv_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #geom_signif(comparisons = list(c("background", "IoD")), y_position=0.007, map_signif_level=TRUE, textsize=6) +
  labs(x="Species", y="π",fill="") +
  scale_x_discrete(labels=c(expression(italic("B.bifarius")), expression(italic("B.vancouverensis")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text=element_text(size=12), axis.title=element_blank(),axis.text.y=element_blank(),legend.position = c(0.8,0.9)) +
  coord_cartesian(ylim = c(0, 0.008))

BbBv_pi_box_plot

## boxplot of PBS
# Check significance
wilcox.test(all_data_20kb$B.s.PBS[bsbi_background], all_data_20kb$B.s.PBS[bsbi_IoD])
wilcox.test(all_data_20kb$B.i.PBS[bsbi_background], all_data_20kb$B.i.PBS[bsbi_IoD])
wilcox.test(all_data_20kb$B.b.PBS[bbbv_background], all_data_20kb$B.b.PBS[bbbv_IoD])
wilcox.test(all_data_20kb$B.v.PBS[bbbv_background], all_data_20kb$B.v.PBS[bbbv_IoD])
# Calculate means

#Bs
Bs_pbs_background_mean <- mean(all_data_20kb$B.s.PBS[bsbi_background])
Bs_pbs_IoD_mean <- mean(all_data_20kb$B.s.PBS[bsbi_IoD])
# Calculate % increase
(Bs_pbs_IoD_mean - Bs_pbs_background_mean)/Bs_pbs_background_mean * 100
# 474% increase

#Bi
Bi_pbs_background_mean <- mean(all_data_20kb$B.i.PBS[bsbi_background])
Bi_pbs_IoD_mean <- mean(all_data_20kb$B.i.PBS[bsbi_IoD])
# Calculate % increase
(Bi_pbs_IoD_mean - Bi_pbs_background_mean)/Bi_pbs_background_mean * 100
# 393% increase


BsBi_pbs_data <- subset(BsBi_all_data_20kbp_long, metric=="PBS")

BsBi_pbs_box_plot <- ggplot(data=BsBi_pbs_data, aes(y=value, x=variable, fill=BsBi_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #geom_signif(comparisons = list(c("background", "IoD")), y_position=0.007, map_signif_level=TRUE, textsize=6) +
  labs(x="Species", y="PBS") +
  scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),axis.title.x=element_blank(),legend.position = "none") +
coord_cartesian(ylim = c(0, 2.5))

BsBi_pbs_box_plot

#BbBv
# Calculate means 

#Bb
Bb_pbs_background_mean <- mean(all_data_20kb$B.b.PBS[bbbv_background])
Bb_pbs_IoD_mean <- mean(all_data_20kb$B.b.PBS[bbbv_IoD])
# Calculate % increase
(Bb_pbs_IoD_mean - Bb_pbs_background_mean)/Bb_pbs_background_mean * 100
# 326% increase

#Bv
Bv_pbs_background_mean <- mean(all_data_20kb$B.v.PBS[bbbv_background])
Bv_pbs_IoD_mean <- mean(all_data_20kb$B.v.PBS[bbbv_IoD])
# Calculate % increase
(Bv_pbs_IoD_mean - Bv_pbs_background_mean)/Bv_pbs_background_mean * 100
# 373% increase

BbBv_pbs_data <- subset(BbBv_all_data_20kbp_long, metric=="PBS")

BbBv_pbs_box_plot <- ggplot(data=BbBv_pbs_data, aes(y=value, x=variable, fill=BbBv_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #geom_signif(comparisons = list(c("background", "IoD")), y_position=0.007, map_signif_level=TRUE, textsize=6) +
  labs(x="Species", y="PBS") +
  scale_x_discrete(labels=c(expression(italic("B.bifarius")), expression(italic("B.vancouverensis")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12),axis.text.y=element_blank(), axis.title=element_blank(),
        axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position = c(0.8,0.8), legend.title=element_blank()) +
coord_cartesian(ylim = c(0, 2.5))

BbBv_pbs_box_plot

ggarrange(BsBi_pbs_box_plot, BbBv_pbs_box_plot,heights = c(1,1),ncol = 2, nrow = 1, align = "v")


## boxplot of dxy

wilcox.test(all_data_20kb$Bs_NQ_dxy[bsNQ_background], all_data_20kb$Bs_NQ_dxy[bsNQ_IoD])

# Calculate means
Bs_NQ_dxy_background_mean <- mean(all_data_20kb$Bs_NQ_dxy[bsNQ_background])
Bs_NQ_dxy_IoD_mean <- mean(all_data_20kb$Bs_NQ_dxy[bsNQ_IoD])
# Calculate % decrease
(Bs_NQ_dxy_background_mean - Bs_NQ_dxy_IoD_mean)/Bs_NQ_dxy_background_mean * 100


BsNQ_dxy_box_plot <- ggplot(data=all_data_20kb, aes(y=Bs_NQ_dxy, x=Bs.N.Q_IoD, fill=Bs.N.Q_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=0.013, map_signif_level=TRUE, textsize=6) +
  labs(x="Within-species", y=expression(paste(bold(d["xy"])))) +
  scale_x_discrete(labels=c("Background", "Islands of Divergence")) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"),axis.text.x=element_blank(), legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.03))

BsNQ_dxy_box_plot


wilcox.test(all_data_20kb$sylv_inco_dxy[bsbi_background], all_data_20kb$sylv_inco_dxy[bsbi_IoD])
# Calculate means
BsBi_dxy_background_mean <- mean(all_data_20kb$sylv_inco_dxy[bsbi_background])
BsBi_dxy_IoD_mean <- mean(all_data_20kb$sylv_inco_dxy[bsbi_IoD])
# Calculate % decrease
(BsBi_dxy_IoD_mean - BsBi_dxy_background_mean)/BsBi_dxy_background_mean * 100

BsBi_dxy_box_plot <- ggplot(data=all_data_20kb, aes(y=sylv_inco_dxy, x=BsBi_IoD, fill=BsBi_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=0.025, map_signif_level=TRUE, textsize=6) +
  labs(x="Sympatric", y=expression(paste(bold(d["xy"])))) +
  scale_x_discrete(labels=c("Background", "Islands of Divergence")) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y=element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"),axis.title.y=element_blank(),
        axis.text.x=element_blank(),legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.03))

BsBi_dxy_box_plot


wilcox.test(all_data_20kb$bifa_vanco_dxy[bbbv_background], all_data_20kb$bifa_vanco_dxy[bbbv_IoD])
# Calculate means
BbBv_dxy_background_mean <- mean(all_data_20kb$bifa_vanco_dxy[bbbv_background])
BbBv_dxy_IoD_mean <- mean(all_data_20kb$bifa_vanco_dxy[bbbv_IoD])
# Calculate % decrease
(BbBv_dxy_background_mean - BbBv_dxy_IoD_mean)/BbBv_dxy_background_mean * 100

BbBv_dxy_box_plot <- ggplot(data=all_data_20kb, aes(y=bifa_vanco_dxy, x=BbBv_IoD, fill=BbBv_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=0.013, map_signif_level=TRUE, textsize=6) +
  labs(x="Allopatric", y=expression(paste(bold(d["xy"])))) +
  scale_x_discrete(labels=c("Background", "Islands of Divergence")) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y=element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.line.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"),axis.title.y=element_blank(),
        axis.text.x=element_blank(),legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.03))

BbBv_dxy_box_plot

ggarrange(BsBi_pi_box_plot,BsBi_pbs_box_plot,BsBi_dxy_box_plot, heights = c(1,1,1),ncol = 1, nrow = 3, align = "v")
ggarrange(BbBv_pi_box_plot,BbBv_pbs_box_plot,BbBv_dxy_box_plot, heights = c(1,1,1),ncol = 1, nrow = 3, align = "v")

# Combined plot of Pi and DXY
ggarrange(BsNQ_pi_box_plot,BsBi_pi_box_plot,BbBv_pi_box_plot,BsNQ_dxy_box_plot,BsBi_dxy_box_plot,BbBv_dxy_box_plot,
          heights = c(1,1,1,1,1,1),ncol = 3, nrow = 2, align = "v")


## boxplot of dxy per chromosome

BsBi_dxy_per_chrom_box_plot <- ggplot(data=all_data_20kb, aes(y=sylv_inco_dxy, x=Bt_LG, fill=BsBi_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  annotate("text", x = 3, y = 0.028, size = 3, label = expression(paste(bolditalic("B. sylvicola - B. incognita")))) +
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y=expression(paste(bold(d["xy"])))) +
  #scale_x_discrete(labels=c("Background", "Islands of Divergence")) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),axis.title.x=element_blank(),legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.03))

BsBi_dxy_per_chrom_box_plot

BbBv_dxy_per_chrom_box_plot <- ggplot(data=all_data_20kb, aes(y=bifa_vanco_dxy, x=Bt_LG, fill=BbBv_IoD)) +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  annotate("text", x = 3, y = 0.014, size = 3, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis")))) +
  labs(x=(expression(italic("B. bifarius - B. vancouverensis"))), y=expression(paste(bold(d["xy"])))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),axis.title.x=element_blank(),legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.015))

BbBv_dxy_per_chrom_box_plot


# Box plots of Fst per chromosome

Bs_Bi_fst_per_chrom_box_plot <- ggplot(data=all_data_20kb, aes(x = Bt_LG, y = Bs.Bi.zfst, fill=BsBi_IoD)) +
  geom_boxplot(outlier.shape=NA) + 
  annotate("text", x = 3, y = 3, size = 3, label = expression(paste(bolditalic("B. sylvicola - B. incognita")))) +
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y=expression(paste(bolditalic(ZF["ST"])))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),axis.title.x=element_blank(),legend.position = "none")
  #coord_cartesian(ylim = c(0, 0.015))
  
Bs_Bi_fst_per_chrom_box_plot

Bb_Bv_fst_per_chrom_box_plot <- ggplot(data=all_data_20kb, aes(x = Bt_LG, y = Bb.Bv.zfst, fill=BbBv_IoD)) +
  geom_boxplot(outlier.shape=NA) + 
  annotate("text", x = 3, y = 5, size = 3, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis")))) +
  labs(x=(expression(italic("B. bifarius - B. vancouverensis"))), y=expression(paste(bolditalic(ZF["ST"])))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),axis.title.x=element_blank(),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.015))

Bb_Bv_fst_per_chrom_box_plot

ggarrange(Bs_Bi_fst_per_chrom_box_plot,BsBi_dxy_per_chrom_box_plot,Bb_Bv_fst_per_chrom_box_plot,BbBv_dxy_per_chrom_box_plot, heights = c(1,1,1,1),ncol = 1, nrow = 4, align = "v")


## boxplot of exon content

# Average exon content
av.exon.content <- mean(all_data_20kb$prop_exon) # 14.0%

BsNQ_IoD_av.exon.content <- mean(all_data_20kb$prop_exon[bsNQ_IoD]) # 11.6 %
BsNQ_background_av.exon.content <- mean(all_data_20kb$prop_exon[bsNQ_background]) # 14.1 %

BsBi_IoD_av.exon.content <- mean(all_data_20kb$prop_exon[bsbi_IoD]) # 14.8 %
BsBi_background_av.exon.content <- mean(all_data_20kb$prop_exon[bsbi_background]) # 13.9 %

BbBv_IoD_av.exon.content <- mean(all_data_20kb$prop_exon[bbbv_IoD]) # 23.8%
BbBv_background_av.exon.content <- mean(all_data_20kb$prop_exon[bbbv_background]) # 12.5%

# Test for significance
wilcox.test(all_data_20kb$prop_exon[bsNQ_background], all_data_20kb$prop_exon[bsNQ_IoD])
wilcox.test(all_data_20kb$prop_exon[bsbi_background], all_data_20kb$prop_exon[bsbi_IoD])
wilcox.test(all_data_20kb$prop_exon[bbbv_background], all_data_20kb$prop_exon[bbbv_IoD])

Bs_NQ_IoD_exon_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=Bs.N.Q_merged_IoD_over_100kb, fill=Bs.N.Q_merged_IoD_over_100kb)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x=(expression(italic("B. sylvicola: Niwot Ridge - Quail Mountain"))), y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

Bs_NQ_IoD_exon_box


Bs_NQ_IoD_exon_per_chrom_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=Bt_LG, fill=Bs.N.Q_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x=(expression(italic("B. sylvicola: Niwot Ridge - Quail Mountain"))), y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

Bs_NQ_IoD_exon_per_chrom_box

BsBi_IoD_exon_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=BsBi_merged_IoD, fill=BsBi_merged_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BsBi_IoD_exon_box

BsBi_IoD_exon_per_chrom_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=Bt_LG, fill=BsBi_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BsBi_IoD_exon_per_chrom_box

BbBv_IoD_exon_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=BbBv_merged_IoD, fill=BbBv_merged_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x=(expression(italic("B. bifarius - B. vancouverensis"))), y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BbBv_IoD_exon_box

BbBv_IoD_exon_per_chrom_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=Bt_LG, fill=BbBv_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x=(expression(italic("B. bifarius - B. vancouverensis"))), y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BbBv_IoD_exon_per_chrom_box

ggarrange(Bs_NQ_IoD_exon_box,BsBi_IoD_exon_box,BbBv_IoD_exon_box, heights = c(1,1,1),ncol = 1, nrow = 3, align = "v")

ggarrange(Bs_NQ_IoD_exon_per_chrom_box, BsBi_IoD_exon_per_chrom_box,BbBv_IoD_exon_per_chrom_box, heights = c(1,1,1),ncol = 1, nrow = 3, align = "v")

## Look at exon content of shared divergent windows
BsBi_BbBv_shared_divergent_windows <- subset(all_data_20kb,shared_divergent_window == "TRUE")
mean(BsBi_BbBv_shared_divergent_windows$prop_exon) # 26.6%

BsBi_BbBv_shared_IoD_exon_box <- ggplot(data=all_data_20kb, aes(y=prop_exon, x=shared_divergent_window, fill=shared_divergent_window)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  labs(x="Common IoD", y="Exon content") +
  #scale_x_discrete(labels=c(expression(italic("B.sylvicola")), expression(italic("B.incognita")))) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BsBi_BbBv_shared_IoD_exon_box

## Number of genes inside and outside of IoDs
mean(all_data_20kb$number_of_overlapping_genes) # 1.26 genes per 20 kbp
mean(BsNQ_IoD$number_of_overlapping_genes) # 1.26 genes per 20 kbp
mean(BsBi_IoD$number_of_overlapping_genes) # 1.38 genes per 20 Kbp
mean(BbBv_IoD$number_of_overlapping_genes) # 2.04 genes per 20 Kbp
mean(BsBi_BbBv_shared_divergent_windows$number_of_overlapping_genes) # 1.79

## Check for difference in mappability in and out of IoDs
wilcox.test(all_data_20kb$mappability[bsNQ_background], all_data_20kb$mappability[bsNQ_IoD])
wilcox.test(all_data_20kb$mappability[bsbi_background], all_data_20kb$mappability[bsbi_IoD])
wilcox.test(all_data_20kb$mappability[bbbv_background], all_data_20kb$mappability[bbbv_IoD])

## Mappability
# Calculate means
within_species_map_background_mean <- mean(all_data_20kb$mappability[bsNQ_background])
within_species_map_IoD_mean <- mean(all_data_20kb$mappability[bsNQ_IoD])
# Calculate % decrease
(within_species_map_background_mean - within_species_map_IoD_mean)/within_species_map_background_mean * 100
# 1.7% reduction

# Calculate means
BsBi_map_background_mean <- mean(all_data_20kb$mappability[bsbi_background])
BsBi_map_IoD_mean <- mean(all_data_20kb$mappability[bsbi_IoD])
# Calculate % decrease
(BsBi_map_background_mean - BsBi_map_IoD_mean)/BsBi_map_background_mean * 100
#1.5% reduction

# Calculate means
BbBv_map_background_mean <- mean(all_data_20kb$mappability[bbbv_background])
BbBv_map_IoD_mean <- mean(all_data_20kb$mappability[bbbv_IoD])
# Calculate % decrease
(BbBv_map_background_mean - BbBv_map_IoD_mean)/BbBv_map_background_mean * 100
#2.6% reduction


## boxplot of mappability

Bs_NQ_IoD_mappability_box <- ggplot(data=all_data_20kb, aes(y=mappability, x=Bs.N.Q_IoD, fill=Bs.N.Q_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=1.1, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. sylvicola: Niwot Ridge-Quail Mountain"))), y="Mappability") +
  scale_y_continuous(expand=c(0,0), limits=c(0.25,1.2)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.text.x=element_blank(),axis.title=element_text(size=12, face="bold"), legend.position = "none", axis.title.x=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

Bs_NQ_IoD_mappability_box


BsBi_IoD_mappability_box <- ggplot(data=all_data_20kb, aes(y=mappability, x=BsBi_IoD, fill=BsBi_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=1.1, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y="Mappability") +
  scale_y_continuous(expand=c(0,0), limits=c(0.25,1.2)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_blank(), legend.position = "none",axis.title=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

BsBi_IoD_mappability_box

BbBv_IoD_mappability_box <- ggplot(data=all_data_20kb, aes(y=mappability, x=BbBv_IoD, fill=BbBv_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=1.1, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. bifarius - B. vancouverensis"))), y="Mappability") +
  scale_y_continuous(expand=c(0,0), limits=c(0.25,1.2)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_blank(), legend.position = "none",axis.title=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

BbBv_IoD_mappability_box

## GC content
# Calculate means
within_species_GC_background_mean <- mean(all_data_20kb$GC_prop[bsNQ_background])
within_species_GC_IoD_mean <- mean(all_data_20kb$GC_prop[bsNQ_IoD])
# Calculate % decrease
(within_species_GC_background_mean - within_species_GC_IoD_mean)/within_species_GC_background_mean * 100
# 4.4% reduction

# Calculate means
BsBi_GC_background_mean <- mean(all_data_20kb$GC_prop[bsbi_background])
BsBi_GC_IoD_mean <- mean(all_data_20kb$GC_prop[bsbi_IoD])
# Calculate % decrease
(BsBi_GC_background_mean - BsBi_GC_IoD_mean)/BsBi_GC_background_mean * 100
#6.3% reduction

# Calculate means
BbBv_GC_background_mean <- mean(all_data_20kb$GC_prop[bbbv_background])
BbBv_GC_IoD_mean <- mean(all_data_20kb$GC_prop[bbbv_IoD])
# Calculate % decrease
(BbBv_GC_background_mean - BbBv_GC_IoD_mean)/BbBv_GC_background_mean * 100
#11.9% reduction

## Check for difference in GC in and out of IoDs
wilcox.test(all_data_20kb$GC_prop[bsNQ_background], all_data_20kb$GC_prop[bsNQ_IoD])
wilcox.test(all_data_20kb$GC_prop[bsbi_background], all_data_20kb$GC_prop[bsbi_IoD])
wilcox.test(all_data_20kb$GC_prop[bbbv_background], all_data_20kb$GC_prop[bbbv_IoD])

## GC content
Bs_NQ_IoD_GC_box <- ggplot(data=all_data_20kb, aes(y=GC_prop, x=Bs.N.Q_IoD, fill=Bs.N.Q_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=0.6, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. sylvicola: Niwot Ridge-Quail Mountain"))), y="Proportion GC") +
  scale_y_continuous(expand=c(0,0), limits=c(0.2,0.7)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.text.x=element_blank(),axis.title=element_text(size=12, face="bold"), legend.position = "none", axis.title.x=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

Bs_NQ_IoD_GC_box

BsBi_IoD_GC_box <- ggplot(data=all_data_20kb, aes(y=GC_prop, x=BsBi_IoD, fill=BsBi_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=0.6, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y="GC content") +
  scale_y_continuous(expand=c(0,0), limits=c(0.2,0.7)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_blank(), legend.position = "none",axis.title=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

BsBi_IoD_GC_box

BbBv_IoD_GC_box <- ggplot(data=all_data_20kb, aes(y=GC_prop, x=BbBv_IoD, fill=BbBv_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=0.6, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. bifarius - B. vancouverensis"))), y="GC content") +
  scale_y_continuous(expand=c(0,0), limits=c(0.2,0.7)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_blank(), legend.position = "none",axis.title=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

BbBv_IoD_GC_box

## Boxplots of RR from LDHat

## Check for difference in RR in and out of IoDs
wilcox.test(all_data_20kb$rho_per_kb[bsNQ_background], all_data_20kb$rho_per_kb[bsNQ_IoD])
wilcox.test(all_data_20kb$rho_per_kb[bsbi_background], all_data_20kb$rho_per_kb[bsbi_IoD])
wilcox.test(all_data_20kb$rho_per_kb[bbbv_background], all_data_20kb$rho_per_kb[bbbv_IoD])

# Calculate means
BsNQ_rho_background_mean <- mean(all_data_20kb$rho_per_kb[bsNQ_background])
BsNQ_rho_IoD_mean <- mean(all_data_20kb$rho_per_kb[bsNQ_IoD])
# Calculate % decrease
(BsNQ_rho_background_mean - BsNQ_rho_IoD_mean)/BsNQ_rho_background_mean * 100
# 83% decrease inside IoDs

# Calculate means
BsBi_rho_background_mean <- mean(all_data_20kb$rho_per_kb[bsbi_background])
BsBi_rho_IoD_mean <- mean(all_data_20kb$rho_per_kb[bsbi_IoD])
# Calculate % decrease
(BsBi_rho_background_mean - BsBi_rho_IoD_mean)/BsBi_rho_background_mean * 100
# 90% decrease inside IoDs

# Calculate means
BbBv_rho_background_mean <- mean(all_data_20kb$rho_per_kb[bbbv_background])
BbBv_rho_IoD_mean <- mean(all_data_20kb$rho_per_kb[bbbv_IoD])
# Calculate % decrease
(BbBv_rho_background_mean - BbBv_rho_IoD_mean)/BbBv_rho_background_mean * 100
# 56% decrease inside IoDs

Bs_NQ_IoD_RR_box <- ggplot(data=all_data_20kb, aes(y=rho_per_kb, x=Bs.N.Q_IoD, fill=Bs.N.Q_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #annotate("text", x = 2.2, y = 128, size = 3, label = expression(bold("Within-species"))) +
  geom_signif(comparisons = list(c("background", "IoD")), y_position=120, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. sylvicola: Niwot Ridge-Quail"))), y="Recombination rate\n(rho/Kbp)") +
  scale_y_continuous(limits=c(0,130)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.text.x=element_blank(),axis.title=element_text(size=12, face="bold"), legend.position = "none", axis.title.x=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

Bs_NQ_IoD_RR_box

BsBi_IoD_RR_box <- ggplot(data=all_data_20kb, aes(y=rho_per_kb, x=BsBi_IoD, fill=BsBi_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #annotate("text", x = 2.2, y = 128, size = 3, label = expression(bold("Sympatric"))) +
  geom_signif(comparisons = list(c("background", "IoD")), y_position=120, map_signif_level=TRUE, textsize=6) +
  labs(x=(expression(italic("B. sylvicola - B. incognita"))), y="Recombination rate (rho/Kbp)") +
  scale_y_continuous(limits=c(0,130)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_blank(), legend.position = "none",axis.title=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

BsBi_IoD_RR_box

BbBv_IoD_RR_box <- ggplot(data=all_data_20kb, aes(y=rho_per_kb, x=BbBv_IoD, fill=BbBv_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  #annotate("text", x = 2.2, y = 128, size = 3, label = expression(bold("Allopatric"))) +
  geom_signif(comparisons = list(c("background", "IoD")), y_position=120, map_signif_level=TRUE, textsize=6) +
  labs(x="Position in genome", y="Recombination rate (rho/Kbp)") +
  scale_y_continuous(limits=c(0,130)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_blank(), legend.position = "none",axis.title=element_blank())
#coord_cartesian(ylim = c(0, 0.008))

BbBv_IoD_RR_box

ggarrange(Bs_NQ_IoD_RR_box,BsBi_IoD_RR_box,BbBv_IoD_RR_box, heights = c(1,1,1),
          ncol = 1, nrow = 3, align = "v")

## Repeat content
# Calculate means
within_species_repeat_background_mean <- mean(all_data_20kb$prop_repeat[bsNQ_background])
within_species_repeat_IoD_mean <- mean(all_data_20kb$prop_repeat[bsNQ_IoD])
# Calculate % decrease
(within_species_repeat_IoD_mean - within_species_repeat_background_mean)/within_species_repeat_background_mean * 100
# 45.6% increase

# Calculate means
BsBi_repeat_background_mean <- mean(all_data_20kb$prop_repeat[bsbi_background])
BsBi_repeat_IoD_mean <- mean(all_data_20kb$prop_repeat[bsbi_IoD])
# Calculate % decrease
(BsBi_repeat_IoD_mean - BsBi_repeat_background_mean)/BsBi_repeat_background_mean * 100
#38.6% increase

# Calculate means
BbBv_repeat_background_mean <- mean(all_data_20kb$prop_repeat[bbbv_background])
BbBv_repeat_IoD_mean <- mean(all_data_20kb$prop_repeat[bbbv_IoD])
# Calculate % decrease
(BbBv_repeat_IoD_mean - BbBv_repeat_background_mean)/BbBv_repeat_background_mean * 100
#55.0% reduction

## Check for difference in repeat content in and out of IoDs
wilcox.test(all_data_20kb$prop_repeat[bsNQ_background], all_data_20kb$prop_repeat[bsNQ_IoD])
wilcox.test(all_data_20kb$prop_repeat[bsbi_background], all_data_20kb$prop_repeat[bsbi_IoD])
wilcox.test(all_data_20kb$prop_repeat[bbbv_background], all_data_20kb$prop_repeat[bbbv_IoD])

## Box plots
Bs_NQ_IoD_repeats_box <- ggplot(data=all_data_20kb, aes(y=prop_repeat, x=Bs.N.Q_IoD, fill=Bs.N.Q_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=1.05, map_signif_level=TRUE, textsize=6) +
  labs(x="Within-species", y="Proportion repeat sequence") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.2)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12),axis.title=element_text(size=12, face="bold"), legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

Bs_NQ_IoD_repeats_box

BsBi_IoD_repeats_box <- ggplot(data=all_data_20kb, aes(y=prop_repeat, x=BsBi_IoD, fill=BsBi_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=1.05, map_signif_level=TRUE, textsize=6) +
  labs(x="Sympatric", y="Repeat content") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.2)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12),axis.text.y=element_blank(),axis.title=element_text(size=12, face="bold"),axis.title.y=element_blank(), legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BsBi_IoD_repeats_box

BbBv_IoD_repeats_box <- ggplot(data=all_data_20kb, aes(y=prop_repeat, x=BbBv_IoD, fill=BbBv_IoD)) +
  geom_boxplot() +
  geom_boxplot(outlier.shape=NA) + # This hides the outliers
  geom_signif(comparisons = list(c("background", "IoD")), y_position=1.05, map_signif_level=TRUE, textsize=6) +
  labs(x="Allopatric", y="Repeat content") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1.2)) +
  scale_fill_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12),axis.text.y=element_blank(),axis.title=element_text(size=12, face="bold"),axis.title.y=element_blank(), legend.position = "none")
#coord_cartesian(ylim = c(0, 0.008))

BbBv_IoD_repeats_box

ggarrange(Bs_NQ_IoD_RR_box,BsBi_IoD_RR_box,BbBv_IoD_RR_box,Bs_NQ_IoD_GC_box,BsBi_IoD_GC_box,BbBv_IoD_GC_box,
          Bs_NQ_IoD_mappability_box,BsBi_IoD_mappability_box,BbBv_IoD_mappability_box,
          Bs_NQ_IoD_repeats_box, BsBi_IoD_repeats_box, BbBv_IoD_repeats_box, heights = c(1,1,1,1,1,1,1,1,1,1,1,1),
          ncol = 3, nrow = 4, align = "v")


### Scatter plots

BsBi_fst_dxy_plot <- ggplot(data=all_data_20kb , aes(x=sylv_inco_dxy, y=Bs.Bi.fst, colour=BsBi_IoD)) +
  geom_point(alpha=0.5) +
  labs(x="B. sylvicola - B. incognita dxy", y="B. sylvicola - B. incognita Fst") +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

BsBi_fst_dxy_plot

BbBv_fst_dxy_plot <- ggplot(data=all_data_20kb , aes(x=bifa_vanco_dxy, y=Bb.Bv.fst, colour=BbBv_IoD)) +
  geom_point(alpha=0.5) +
  labs(x="B. bifarius - B. vancouverensis dxy", y="B. bifarius - B. vancouverensis Fst") +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

BbBv_fst_dxy_plot

Bs.PBS_dxy_plot <- ggplot(data=all_data_20kb , aes(x=sylv_inco_dxy, y=B.s.PBS, colour=BsBi_IoD)) +
  geom_point() +
  labs(x="B. sylvicola - B. incognita dxy", y="B. sylvicola PBS") +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

Bs.PBS_dxy_plot

Bi.PBS_dxy_plot <- ggplot(data=all_data_20kb , aes(x=sylv_inco_dxy, y=B.i.PBS, colour=BsBi_IoD)) +
  geom_point() +
  labs(x="B. sylvicola - B. incognita dxy", y="B. incognita PBS") +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

Bi.PBS_dxy_plot

ggarrange(Bs.PBS_dxy_plot,Bi.PBS_dxy_plot,fst_dxy_plot,Bi.PBS_Bs.PBS_plot, heights = c(1,1,1,1),
          ncol = 2, nrow = 2, align = "v")

Bi.PBS_Bs.PBS_plot <- ggplot(data=all_data_20kb , aes(x=B.s.PBS, y=B.i.PBS)) +
  geom_pointdensity() +
  scale_color_viridis() +
  scale_y_continuous(limits=c(-0.5,1.8)) +
  scale_x_continuous(limits=c(-1.0,2.5)) +
  labs(x="B. sylvicola PBS", y="B. incognita PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position = "none")

Bi.PBS_Bs.PBS_plot

Bb.PBS_Bv.PBS_plot <- ggplot(data=all_data_20kb , aes(x=B.b.PBS, y=B.v.PBS)) +
  geom_pointdensity() +
  scale_color_viridis() +
  scale_y_continuous(limits=c(-0.5,1.8)) +
  scale_x_continuous(limits=c(-1.0,2.5)) +
  labs(x="B. bifarius PBS", y="B. vancouverensis PBS") +
  #scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"),legend.position = "none")

Bb.PBS_Bv.PBS_plot

cor.test(all_data_20kb$B.s.PBS, all_data_20kb$B.i.PBS, method=c("spearman"))
cor.test(all_data_20kb$B.b.PBS, all_data_20kb$B.v.PBS, method=c("spearman"))


ggarrange(Bi.PBS_Bs.PBS_plot,Bb.PBS_Bv.PBS_plot, heights = c(1,1),
          ncol = 1, nrow = 2, align = "v")

cor.test(all_data_20kb$B.i.PBS, all_data_20kb$Bs.Bi.fst, method=c("spearman"))
cor.test(all_data_20kb$B.v.PBS, all_data_20kb$Bb.Bv.fst, method=c("spearman"))

## Compare GC and RR

## From linkage map

cor.test(all_data_20kb$RR, all_data_20kb$GC_prop, method=c("spearman"))


GC_RR_plot <- ggplot(data=all_data_20kb , aes(x=RR, y=GC_prop)) +
  geom_point() +
  labs(x="Recombination rate", y="Proportion GC") +
  geom_smooth() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position = "none")

GC_RR_plot

## From LDHat
# Bpen-1
cor.test(all_data_20kb$rho_per_kb, all_data_20kb$GC_prop, method=c("spearman"))
cor.test(all_data_20kb$rho_per_kb, all_data_20kb$Bs.N.Q.zfst, method=c("spearman"))
cor.test(all_data_20kb$rho_per_kb, all_data_20kb$Bs.Bi.zfst, method=c("spearman"))
cor.test(all_data_20kb$rho_per_kb, all_data_20kb$Bb.Bv.zfst, method=c("spearman"))
cor.test(all_data_20kb$rho_per_kb, all_data_20kb$sylv_inco_dxy, method=c("spearman"))

GC_RR_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=GC_prop)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="Recombination rate (ρ/kbp)", y="Proportion GC") +
  geom_smooth() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="none")

GC_RR_plot

RR_Bs_Pi_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=sylv_PI, colour=BsBi_IoD)) +
  geom_point(alpha=0.5) +
  labs(x="Recombination rate (ρ/kbp)", y=(expression(bolditalic("B. sylvicola π")))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_Bs_Pi_plot

RR_Bi_Pi_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=inco_PI, colour=BsBi_IoD)) +
  geom_point(alpha=0.5) +
  labs(x="Recombination rate (ρ/kbp)", y=(expression(bolditalic("B. incognita π")))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_Bi_Pi_plot

RR_Bb_Pi_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=bifa_PI, colour=BbBv_IoD)) +
  geom_point(alpha=0.5) +
  labs(x="Recombination rate (ρ/kbp)", y=(expression(bolditalic("B. bifarius π")))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_Bb_Pi_plot

RR_Bv_Pi_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=vanco_PI, colour=BbBv_IoD)) +
  geom_point(alpha=0.5) +
  labs(x="Recombination rate (ρ/kbp)", y=(expression(bolditalic("B. vancouverensis π")))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_Bv_Pi_plot

RR_Bs.N.Q.fst_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=Bs.N.Q.zfst, colour=Bs.N.Q_IoD)) +
  geom_point(alpha=0.5) +
  annotate("text", x = 70, y = 15, size = 3, label = expression(paste(bolditalic("B. sylvicola: Niwot Ridge - Quail Mountain")))) +
  labs(x="Recombination rate (ρ/kbp)", y=expression(paste(bolditalic(ZF["ST"])))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_Bs.N.Q.fst_plot

RR_BsBi_zfst_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=Bs.Bi.zfst, colour=BsBi_IoD)) +
  geom_point(alpha=0.5) +
  annotate("text", x = 70, y = 3, size = 3, label = expression(paste(bolditalic("B. sylvicola - B. incognita")))) +
  labs(x="Recombination rate (ρ/kbp)", y=expression(paste(bolditalic(ZF["ST"])))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_BsBi_zfst_plot

RR_BbBv_zfst_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=Bb.Bv.zfst, colour=BbBv_IoD)) +
  geom_point(alpha=0.5) +
  annotate("text", x = 70, y = 6, size = 3, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis")))) +
  labs(x="Recombination rate (ρ/kbp)", y=expression(paste(bolditalic(ZF["ST"])))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_BbBv_zfst_plot


RR_BsBi_dxy_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=sylv_inco_dxy, colour=BsBi_IoD)) +
  geom_point(alpha=0.5) +
  annotate("text", x = 70, y = 0.025, size = 3, label = expression(paste(bolditalic("B. sylvicola - B. incognita")))) +
  labs(x="Recombination rate (ρ/kbp)", y=expression(paste(bolditalic(d["xy"])))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_BsBi_dxy_plot

RR_BbBv_dxy_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb, y=bifa_vanco_dxy, colour=BbBv_IoD)) +
  geom_point(alpha=0.5) +
  annotate("text", x = 70, y = 0.015, size = 3, label = expression(paste(bolditalic("B. bifarius - B. vancouverensis")))) +
  labs(x="Recombination rate (ρ/kbp)", y=expression(paste(bolditalic(d["xy"])))) +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), 
        axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position = "none")

RR_BbBv_dxy_plot

ggarrange(RR_Bs_Pi_plot,RR_Bb_Pi_plot,RR_Bi_Pi_plot,RR_Bv_Pi_plot,
          RR_BsBi_zfst_plot,RR_BbBv_zfst_plot,RR_BsBi_dxy_plot,RR_BbBv_dxy_plot,
          GC_RR_plot, heights = c(1,1,1,1,1,1,1,1,1),ncol = 2, nrow = 5, align = "v")


## RR acros the genome
RR_Bpen_1_genome_plot <- ggplot(data=all_data_20kb , aes(y=rho_per_kb, x=GEN_START)) +
  geom_vline(data = chrom_lengths, aes(xintercept = chrom_end), colour = "grey", alpha=0.5, lty =2) +
  geom_line() +
  labs(x="Genome position (Mbp)", y="Recombination rate (ρ/kbp)") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,50000000,100000000,150000000,200000000,250000000), labels=c(0,50, 100, 150, 200, 250), limits=c(0,252081862)) +
  scale_y_continuous(expand=c(0,0)) +
  #scale_colour_manual(values=c("purple","green","red","darkgrey")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold"),legend.position = "none")
RR_Bpen_1_genome_plot

ggarrange(Bs_Bi_fst_plot,RR_Bpen_1_genome_plot, heights = c(2,1),
          ncol = 1, nrow = 2, align = "v")

cor.test(all_data_20kb$rho_per_kb, all_data_20kb$Bs.Bi.fst, method=c("pearson"))

fst_RR_Bpen_1_plot <- ggplot(data=all_data_20kb , aes(x=rho_per_kb.x, y=Bs.Bi.fst)) +
  geom_point() +
  labs(x="Recombination rate (ρ/kbp)", y="FST") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

fst_RR_Bpen_1_plot

cor.test(all_data_20kb$rho_per_kb, all_data_20kb$sylv_PI, method=c("pearson"))
cor.test(all_data_20kb$rho_per_kb, all_data_20kb$sylv_inco_dxy, method=c("pearson"))
cor.test(all_data_20kb$Bs.Bi.fst, all_data_20kb$sylv_inco_dxy, method=c("pearson"))

Bs_NQ_fst_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=Bs.N.Q.fst)) +
  geom_point() +
  labs(x="GC prop", y="FST") +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

Bs_NQ_fst_GC_plot
cor.test(all_data_20kb$Bs.N.Q.fst, all_data_20kb$GC_prop, method=c("pearson"))

BsBi_fst_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=Bs.Bi.fst)) +
  geom_point() +
  labs(x="GC prop", y="FST") +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

BsBi_fst_GC_plot
cor.test(all_data_20kb$Bs.Bi.fst, all_data_20kb$GC_prop, method=c("pearson"))

BbBv_fst_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=Bb.Bv.fst)) +
  geom_point() +
  labs(x="GC prop", y="FST") +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

BbBv_fst_GC_plot
cor.test(all_data_20kb$Bb.Bv.fst, all_data_20kb$GC_prop, method=c("pearson"))

## Fst - repeats correlations
cor.test(all_data_20kb$Bs.N.Q.fst, all_data_20kb$prop_repeat, method=c("pearson"))
cor.test(all_data_20kb$Bs.Bi.fst, all_data_20kb$prop_repeat, method=c("pearson"))
cor.test(all_data_20kb$Bb.Bv.fst, all_data_20kb$prop_repeat, method=c("pearson"))


fst_mappability_plot <- ggplot(data=all_data_20kb , aes(x=mappability, y=Bs.N.Q.fst)) +
  geom_point() +
  labs(x="mappability", y="FST") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

fst_mappability_plot

Bb_Bv_fst_mappability_plot <- ggplot(data=all_data_20kb , aes(x=mappability, y=Bb.Bv.fst,colour=BbBv_IoD)) +
  geom_point() +
  labs(x="mappability", y="FST") +
  scale_colour_manual(values=c("darkgrey","red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

Bb_Bv_fst_mappability_plot

DoC_mappability_plot <- ggplot(data=all_data_20kb , aes(x=DoC_Bsyl, y=mappability, colour=IoD)) +
  geom_point() +
  labs(x="Depth of coverage", y="mappability") +
  scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

DoC_mappability_plot

repeats_mappability_plot <- ggplot(data=all_data_20kb , aes(x=prop_repeat, y=mappability)) +
  geom_point() +
  labs(x="Repeat content", y="mappability") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

repeats_mappability_plot


Bs_PBS_mappability_plot <- ggplot(data=all_data_20kb , aes(x=mappability, y=sylv_PI, colour=IoD)) +
  geom_point() +
  labs(x="mappability", y="B.sylvicola PBS") +
  scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

Bs_PBS_mappability_plot

Bv_PBS_mappability_plot <- ggplot(data=all_data_20kb , aes(x=mappability, y=B.v.PBS, colour=IoD)) +
  geom_point() +
  labs(x="mappability", y="B.vancouverensis PBS") +
  #scale_colour_manual(values=c("blue","purple","green","red","darkgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

Bv_PBS_mappability_plot


## GC - PBS correlations

Bs_PBS_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=B.s.PBS)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="GC content", y="B.sylvicola PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="blank")

Bs_PBS_GC_plot

Bi_PBS_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=B.i.PBS)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="GC content", y="B.incognita PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="blank")

Bi_PBS_GC_plot

Bb_PBS_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=B.b.PBS)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="GC content", y="B.bifarius PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="blank")

Bb_PBS_GC_plot


Bv_PBS_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=B.v.PBS)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="GC content", y="B.vancouverensis PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="blank")

Bv_PBS_GC_plot

ggarrange(Bs_PBS_GC_plot,Bi_PBS_GC_plot,Bb_PBS_GC_plot,Bv_PBS_GC_plot, heights = c(1,1,1,1),
          ncol = 2, nrow = 2, align = "v")

# GC - DXY scatter plots

cor.test(all_data_20kb$GC_prop, all_data_20kb$sylv_inco_dxy, method=c("pearson"))
BsBi_dxy_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=sylv_inco_dxy)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="GC content", y="B. sylvicola - B. incognita dxy") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="blank")

BsBi_dxy_GC_plot

cor.test(all_data_20kb$GC_prop, all_data_20kb$bifa_vanco_dxy, method=c("pearson"))
BbBv_dxy_GC_plot <- ggplot(data=all_data_20kb , aes(x=GC_prop, y=bifa_vanco_dxy)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(x="GC content", y="B. bifarius - B. vancouverensis dxy") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"), legend.position="blank")

BbBv_dxy_GC_plot

ggarrange(BsBi_dxy_GC_plot,BbBv_dxy_GC_plot, heights = c(1,1), ncol = 1, nrow = 2, align = "v")

## Compare number of snps per window between comparisons

BsBi_BbBv_snps_plot <- ggplot(data=all_data_20kb , aes(x=Bs_Bi_number_of_snps, y=Bb_Bv_number_of_snps)) +
  geom_point() +
  #labs(x="GC content", y="B.bifarius PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

BsBi_BbBv_snps_plot
cor.test(all_data_20kb$Bs_Bi_number_of_snps, all_data_20kb$Bb_Bv_number_of_snps, method=c("pearson"))
# Corr = 0.6568656 


# Number of SNPs and mappability
snps_v_mappability_plot <- ggplot(data=all_data_20kb , aes(x=mappability, y=Bb_Bv_number_of_snps)) +
  geom_point() +
  #labs(x="GC content", y="B.bifarius PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

snps_v_mappability_plot
cor.test(all_data_20kb$mappability, all_data_20kb$Bb_Bv_number_of_snps, method=c("pearson"))
# cor =0.3100051

# Depth v number of SNPs
snps_v_depth_plot <- ggplot(data=all_data_20kb , aes(x=DoC_Bsyl, y=Bs_Bi_number_of_snps)) +
  geom_point() +
  #labs(x="GC content", y="B.bifarius PBS") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12), axis.title=element_text(size=16, face="bold"))

snps_v_depth_plot
cor.test(all_data_20kb$DoC_Bsyl, all_data_20kb$Bs_Bi_number_of_snps, method=c("pearson"))


## Output window coords for assigning RR

window_coords <- data.frame(all_data_20kb$Bt_LG, all_data_20kb$GEN_START)
write.table(window_coords, file="window_genome_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

## Output DXY, PI, and ZFST values for looking at change from peaks

# dxy
dxy_coords <- data.frame(all_data_20kb$Bt_LG, all_data_20kb$GEN_START, all_data_20kb$Bs_NQ_dxy, all_data_20kb$sylv_inco_dxy, all_data_20kb$bifa_vanco_dxy)
write.table(dxy_coords, file="dxy_genome_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

# ZFST
zfst_coords <- data.frame(all_data_20kb$Bt_LG, all_data_20kb$GEN_START, all_data_20kb$Bs.N.Q.zfst, all_data_20kb$Bs.Bi.zfst, all_data_20kb$Bb.Bv.zfst)
write.table(zfst_coords, file="zfst_genome_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Pi
pi_coords <- data.frame(all_data_20kb$Bt_LG, all_data_20kb$GEN_START, all_data_20kb$Bs_niw_PI,all_data_20kb$Bs_qua_PI,all_data_20kb$sylv_PI, all_data_20kb$inco_PI, all_data_20kb$bifa_PI, all_data_20kb$vanco_PI)
write.table(pi_coords, file="pi_genome_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)



# RR
RR_coords <- data.frame(all_data_20kb$CHROM, all_data_20kb$BIN_START, all_data_20kb$rho_per_kb)
RR_coords$all_data_20kb.BIN_START = RR_coords$all_data_20kb.BIN_START-1
write.table(RR_coords, file="RR_coords.txt", sep="\t", row.names=FALSE, quote=FALSE)

