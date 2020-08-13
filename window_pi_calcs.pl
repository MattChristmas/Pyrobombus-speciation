#!/usr/bin/perl -w

use strict; #enforce some good programming rules
use warnings;
use diagnostics;
use Data::Dumper; 
use IO::Handle;


## First get cumulative contig lengths for genome positions

open IN1, "sylv_genome_contig_cum_lengths.txt" or die; #This file contains two columns, listing the contig names and cumulative lengths for the assembly

my %lengths;

while (my $line = <IN1>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);	
	
	my $contig = $data[0];
	my $length = $data[1];
	
	#print "$contig\t$length\n";
	
	$lengths{$contig} = $length;
									
}

close IN1;

#print "Done with setting contig lengths\n";

open IN2, "maker_gene_coords_no_overlap.txt" or die;

my %gene_coords;
my $gene=1;

while (my $line = <IN2>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);	
	
	my $contig = $data[0];
	my $start = $data[3];
	my $stop = $data[4];
	my $coords = "$start-$stop";
	
	#print "$contig\t$length\n";
	
	$gene_coords{$contig}->{$gene}=$coords; #This will create a hash where the (arbitrary) gene number is the key and the coordinates are the value.
											# The reference to this hash is stored as the value of the scaffold key in the gene_coords hash
											
	$gene++;
										
}

#print "Done with reading in gff file\n";

# my %contig_001_hash = %{ $gene_coords{"contig_001"} };  # FOR TESTING
# 
# foreach (sort keys %contig_001_hash) {
#     print "$_ : $contig_001_hash{$_}\n";
#   }



open IN3, "Bs.Bi.Bm.Bb_whole_genome_snps_filtered_20pc_missingness.recode.vcf" or die;

open my ($output_fh), ">>bifa_pi_20kb_windows_intergenic.txt" or die;
$output_fh->autoflush(1); # turn off buffering

print $output_fh "contig\tcontig_start\tcontig_mid\tcontig_end\tGEN_START\tGEN_MID\tGEN_END\tbifa_intergen_pi\n";


my %samples = ( # Set up a sample array where the key is the position of the sample in the vcf file
	9 => '108_S_Bor_4_W',
	10 => '120_S_Bor_6_W',
	11 => '121_S_Bor_6_W',
	12 => '132_S_Bor_6_W',
	13 => '134_S_Bor_6_W',
	14 => '139_S_Bor_6_W',
	15 => '142_S_Bor_3_W',
	16 => '146_S_Bor_3_W',
	17 => '149_S_Bor_3_W',
	18 => '150_S_Bor_3_W',
	19 => '151_S_Bor_3_W',
	20 => '181_S_Bor_8_W',
	21 => '191_S_Bor_8_W',
	22 => '193_S_Bor_8_W',
	23 => '216_S_Hor_9_W',
	24 => '217_S_Hor_9_W',
	25 => '218_S_Hor_9_W',
	26 => '219_S_Hor_9_W',
	27 => '220_S_Hor_9_W',
	28 => '221_S_Hor_9_W',
	29 => '223_S_Hor_9_W',
	30 => '226_S_Hor_9_W',
	31 => '228_S_Hor_9_W',
	32 => '230_S_Hor_9_W',
	33 => '232_S_Hor_9_W',
	34 => '234_S_Hor_9_W',
	35 => '235_S_Hor_9_W',
	36 => '236_S_Hor_9_W',
	37 => '237_S_Hor_9_W',
	38 => '238_S_Hor_9_W',
	39 => '240_S_Hor_9_W',
	40 => '241_S_Hor_9_W',
	41 => '242_S_Hor_9_W',
	42 => '243_S_Hor_9_W',
	43 => '244_S_Hor_9_W',
	44 => '245_S_Hor_9_W',
	45 => '246_S_Hor_9_W',
	46 => '249_S_Hor_9_W',
	47 => '250_S_Hor_9_W',
	48 => '260_S_Hor_10_W',
	49 => '262_S_Hor_10_W',
	50 => '263_S_Hor_10_W',
	51 => '265_S_Hor_10_W',
	52 => '266_S_Hor_10_W',
	53 => '267_S_Hor_10_W',
	54 => '268_S_Hor_10_W',
	55 => '269_S_Hor_10_W',
	56 => '270_S_Hor_10_W',
	57 => '274_S_Hor_10_W',
	58 => '275_S_Hor_10_W',
	59 => '277_S_Hor_10_W',
	60 => '285_S_Hor_10_W',
	61 => '288_S_Hor_10_W',
	62 => '290_S_Hor_10_W',
	63 => '292_S_Hor_10_W',
	64 => '293_S_Hor_10_W',
	65 => '294_S_Hor_10_W',
	66 => '297_S_Hor_10_W',
	67 => '298_S_Hor_10_W',
	68 => '300_S_Hor_10_W',
	69 => '304_S_Hor_10_W',
	70 => '306_S_Hor_10_W',
	71 => '307_S_Hor_10_W',
	72 => '309_S_Hor_10_W',
	73 => '311_S_Hor_10_W',
	74 => '313_S_Hor_10_W',
	75 => '314_S_Hor_11_W',
	76 => '315_S_Hor_11_W',
	77 => '320_S_Hor_11_W',
	78 => '323_S_Hor_11_W',
	79 => '329_S_Hor_11_W',
	80 => '330_S_Hor_11_W',
	81 => '334_S_Hor_11_W',
	82 => '342_S_Hor_11_W',
	83 => '344_S_Hor_11_W',
	84 => '350_S_Hor_11_W',
	85 => '352_S_Dem_12_W',
	86 => '358_S_Dem_12_W',
	87 => '370_S_Dem_12_W',
	88 => '388_S_Dem_12_W',
	89 => '38_S_Dem_1_W',
	90 => '397_S_Dem_13_W',
	91 => '401_S_Dem_13_W',
	92 => '412_S_Dem_13_W',
	93 => '419_S_Dem_13_W',
	94 => '423_S_Dem_13_W',
	95 => '425_S_Dem_13_W',
	96 => '431_S_Dem_13_W',
	97 => '441_S_Eva_14_W',
	98 => '442_S_Eva_14_W',
	99 => '443_S_Eva_14_W',
	100 => '444_S_Eva_14_W',
	101 => '445_S_Eva_14_W',
	102 => '446_S_Eva_14_W',
	103 => '447_S_Eva_14_W',
	104 => '450_S_Eva_14_W',
	105 => '451_S_Eva_14_W',
	106 => '453_S_Eva_14_W',
	107 => '454_S_Eva_14_W',
	108 => '456_S_Eva_14_W',
	109 => '458_S_Eva_14_W',
	110 => '45_S_Dem_1_W',
	111 => '463_S_Eva_14_W',
	112 => '464_S_Eva_14_W',
	113 => '465_S_Eva_14_W',
	114 => '466_S_Eva_14_W',
	115 => '467_S_Eva_14_W',
	116 => '468_S_Eva_14_W',
	117 => '470_S_Eva_14_W',
	118 => '471_S_Eva_14_W',
	119 => '473_S_Eva_14_W',
	120 => '474_S_Eva_14_W',
	121 => '476_S_Eva_14_W',
	122 => '47_S_Dem_1_W',
	123 => '480_S_Eva_14_W',
	124 => '484_S_Eva_15_W',
	125 => '48_S_Dem_1_W',
	126 => '492_S_Eva_15_W',
	127 => '493_S_Eva_15_W',
	128 => '499_S_Eva_15_W',
	129 => '49_S_Dem_1_W',
	130 => '502_S_Eva_15_W',
	131 => '504_S_Eva_15_W',
	132 => '509_S_Eva_15_W',
	133 => '50_S_Dem_1_Q',
	134 => '510_S_Eva_15_W',
	135 => '512_S_Eva_15_W',
	136 => '518_S_Eva_15_W',
	137 => '519_S_Eva_15_W',
	138 => '51_S_Dem_1_W',
	139 => '524_S_Eva_16_W',
	140 => '527_S_Eva_16_W',
	141 => '528_S_Eva_16_W',
	142 => '529_S_Eva_16_W',
	143 => '530_S_Eva_16_W',
	144 => '534_S_Eva_16_W',
	145 => '535_S_Eva_16_W',
	146 => '536_S_Eva_16_W',
	147 => '538_S_Eva_16_W',
	148 => '541_S_Eva_16_W',
	149 => '546_S_Eva_16_W',
	150 => '54_S_Dem_1_Q',
	151 => '552_S_Eva_16_W',
	152 => '555_S_Eva_16_W',
	153 => '559_S_Eva_16_W',
	154 => '55_S_Dem_1_W',
	155 => '560_S_Eva_16_W',
	156 => '561_S_Eva_16_W',
	157 => '563_S_Eva_16_W',
	158 => '56_S_Dem_1_W',
	159 => '571_S_Pen_17_W',
	160 => '574_S_Pen_17_W',
	161 => '576_S_Pen_17_W',
	162 => '577_S_Pen_17_W',
	163 => '586_S_Pen_17_W',
	164 => '58_S_Dem_1_W',
	165 => '592_S_Pen_17_W',
	166 => '594_S_Pen_17_W',
	167 => '598_S_Pen_17_W',
	168 => '599_S_Pen_17_W',
	169 => '601_S_Pen_17_W',
	170 => '60_S_Dem_1_W',
	171 => '612_S_Pen_18_W',
	172 => '618_S_Pen_18_W',
	173 => '61_S_Dem_1_W',
	174 => '621_S_Pen_18_W',
	175 => '623_S_Pen_18_W',
	176 => '624_S_Pen_18_W',
	177 => '625_S_Pen_18_W',
	178 => '626_S_Pen_18_W',
	179 => '629_S_Pen_18_W',
	180 => '632_S_Pen_18_W',
	181 => '633_S_Pen_18_W',
	182 => '634_S_Pen_18_W',
	183 => '637_S_Pen_18_W',
	184 => '639_S_Pen_18_W',
	185 => '641_S_Pen_18_W',
	186 => '64_S_Dem_1_W',
	187 => '678_S_Pen_19_W',
	188 => '679_S_Pen_19_W',
	189 => '67_S_Dem_1_W',
	190 => '682_S_Pen_19_W',
	191 => '683_S_Pen_19_W',
	192 => '684_S_Pen_19_W',
	193 => '687_S_Pen_19_W',
	194 => '688_S_Pen_19_W',
	195 => '689_S_Pen_19_W',
	196 => '690_S_Pen_19_W',
	197 => '691_S_Pen_19_W',
	198 => '699_S_Pen_20_W',
	199 => '69_S_Dem_1_W',
	200 => '702_S_Pen_20_W',
	201 => '709_S_Pen_20_W',
	202 => '70_S_Dem_1_W',
	203 => '719_S_Qua_21_W',
	204 => '721_S_Qua_21_W',
	205 => '723_S_Qua_21_W',
	206 => '727_S_Qua_21_W',
	207 => '72_S_Bor_2_W',
	208 => '730_S_Qua_21_W',
	209 => '733_S_Qua_21_W',
	210 => '739_S_Qua_21_W',
	211 => '73_S_Bor_2_W',
	212 => '740_S_Qua_21_W',
	213 => '743_S_Qua_21_W',
	214 => '74_S_Bor_2_W',
	215 => '751_S_Qua_22_W',
	216 => '756_S_Qua_22_W',
	217 => '757_S_Qua_22_W',
	218 => '758_S_Qua_22_W',
	219 => '75_S_Bor_2_W',
	220 => '765_S_Qua_22_W',
	221 => '766_S_Qua_22_W',
	222 => '769_S_Qua_22_W',
	223 => '76_S_Bor_2_W',
	224 => '770_S_Qua_22_W',
	225 => '773_S_Qua_22_W',
	226 => '774_S_Qua_22_W',
	227 => '775_S_Qua_22_W',
	228 => '778_S_Qua_22_W',
	229 => '77_S_Bor_2_W',
	230 => '780_S_Qua_22_W',
	231 => '781_S_Qua_22_W',
	232 => '782_S_Qua_22_W',
	233 => '783_S_Qua_22_W',
	234 => '784_S_Qua_22_W',
	235 => '785_S_Qua_22_W',
	236 => '787_S_Qua_22_W',
	237 => '788_S_Qua_22_W',
	238 => '794_S_Qua_23_W',
	239 => '79_S_Bor_2_W',
	240 => '804_S_Qua_23_W',
	241 => '806_S_Qua_23_W',
	242 => '807_S_Qua_23_W',
	243 => '808_S_Qua_23_W',
	244 => '80_S_Bor_2_W',
	245 => '812_S_Qua_23_W',
	246 => '815_S_Qua_23_W',
	247 => '816_S_Qua_23_W',
	248 => '817_S_Qua_23_W',
	249 => '81_S_Bor_2_W',
	250 => '824_S_Qua_23_W',
	251 => '83_S_Bor_2_W',
	252 => '852_S_Niw_25_W',
	253 => '854_S_Niw_25_W',
	254 => '855_S_Niw_25_W',
	255 => '859_S_Niw_25_W',
	256 => '85_S_Bor_2_W',
	257 => '860_S_Niw_25_W',
	258 => '863_S_Niw_25_W',
	259 => '864_S_Niw_25_W',
	260 => '865_S_Niw_25_W',
	261 => '869_S_Niw_25_W',
	262 => '86_S_Bor_2_W',
	263 => '870_S_Niw_25_W',
	264 => '878_S_Niw_25_W',
	265 => '880_S_Niw_25_W',
	266 => '883_S_Niw_25_W',
	267 => '884_S_Niw_25_W',
	268 => '885_S_Niw_25_W',
	269 => '887_S_Niw_25_W',
	270 => '888_S_Niw_25_W',
	271 => '892_S_Niw_25_W',
	272 => '895_S_Niw_25_W',
	273 => '898_S_Niw_25_W',
	274 => '899_S_Niw_25_W',
	275 => '8_S_Dem_1_W',
	276 => '904_S_Niw_25_W',
	277 => '905_S_Niw_25_W',
	278 => '906_S_Niw_25_W',
	279 => '908_S_Niw_25_W',
	280 => '909_S_Niw_25_W',
	281 => '910_S_Niw_26_W',
	282 => '912_S_Niw_26_W',
	283 => '914_S_Niw_26_W',
	284 => '915_S_Niw_26_W',
	285 => '916_S_Niw_26_W',
	286 => '917_S_Niw_26_W',
	287 => '918_S_Niw_26_W',
	288 => '919_S_Niw_26_W',
	289 => '91_S_Bor_2_W',
	290 => '920_S_Niw_26_W',
	291 => '921_S_Niw_26_W',
	292 => '923_S_Niw_26_W',
	293 => '924_S_Niw_26_W',
	294 => '928_S_Niw_26_W',
	295 => '929_S_Niw_27_W',
	296 => '930_S_Niw_27_W',
	297 => '931_S_Niw_27_W',
	298 => '932_S_Niw_27_W',
	299 => '933_S_Niw_27_W',
	300 => '934_S_Niw_27_W',
	301 => '935_S_Niw_27_W',
	302 => '937_S_Niw_27_W',
	303 => '945_S_Niw_27_W',
	304 => '946_S_Niw_27_W',
	305 => '948_S_Niw_27_W',
	306 => '949_S_Niw_27_W',
	307 => '94_S_Bor_2_W',
	308 => '950_S_Niw_27_W',
	309 => '952_S_Niw_27_W',
	310 => 'SRR8700077_black',
	311 => 'SRR8700078_black',
	312 => 'SRR8700079_black',
	313 => 'SRR8700080_black',
	314 => 'SRR8700081_black',
	315 => 'SRR8700082_black',
	316 => 'SRR8700083_black',
	317 => 'SRR8700084_black',
	318 => 'SRR8700085_black',
	319 => 'SRR8700086_black',
	320 => 'SRR8700087_red',
	321 => 'SRR8700088_red',
	322 => 'SRR8700089_red',
	323 => 'SRR8700090_red',
	324 => 'SRR8700091_red',
	325 => 'SRR8700092_red',
	326 => 'SRR8700093_red',
	327 => 'SRR8700094_red',
	328 => 'SRR8700095_red',
	329 => 'SRR8700096_red',
	330 => 'SRR8700097_red',
	);

# my %sylv = (
# 	'108_S_Bor_4_W' => 0,
# 	'149_S_Bor_3_W' => 0,
# 	'223_S_Hor_9_W' => 0,
# 	'290_S_Hor_10_W' => 0,
# 	'352_S_Dem_12_W' => 0,
# 	'397_S_Dem_13_W' => 0,
# 	'45_S_Dem_1_W' => 0,
# 	'450_S_Eva_14_W' => 0,
# 	'48_S_Dem_1_W' => 0,
# 	'480_S_Eva_14_W' => 0,
# 	'484_S_Eva_15_W' => 0,
# 	'49_S_Dem_1_W' => 0,
# 	'510_S_Eva_15_W' => 0,
# 	'534_S_Eva_16_W' => 0,
# 	'599_S_Pen_17_W' => 0,
# 	'60_S_Dem_1_W' => 0,
# 	'601_S_Pen_17_W' => 0,
# 	'61_S_Dem_1_W' => 0,
# 	'743_S_Qua_21_W' => 0,
# 	'75_S_Bor_2_W' => 0,
# 	'756_S_Qua_22_W' => 0,
# 	'766_S_Qua_22_W' => 0,
# 	'77_S_Bor_2_W' => 0,
# 	'780_S_Qua_22_W' => 0,
# 	'81_S_Bor_2_W' => 0,
# 	'812_S_Qua_23_W' => 0,
# 	'855_S_Niw_25_W' => 0,
# 	'899_S_Niw_25_W' => 0,
# 	'910_S_Niw_26_W' => 0,
# 	'952_S_Niw_27_W' => 0,
# 	);

# my %inco = (
# 	'120_S_Bor_6_W' => 0,
# 	'132_S_Bor_6_W' => 0,
# 	'139_S_Bor_6_W' => 0,
# 	'146_S_Bor_3_W' => 0,
# 	'191_S_Bor_8_W' => 0,
# 	'262_S_Hor_10_W' => 0,
# 	'293_S_Hor_10_W' => 0,
# 	'320_S_Hor_11_W' => 0,
# 	'350_S_Hor_11_W' => 0,
# 	'527_S_Eva_16_W' => 0,
# 	'612_S_Pen_18_W' => 0,
# 	'621_S_Pen_18_W' => 0,
# 	'624_S_Pen_18_W' => 0,
# 	'633_S_Pen_18_W' => 0,
# 	'679_S_Pen_19_W' => 0,
# 	'684_S_Pen_19_W' => 0,
# 	'689_S_Pen_19_W' => 0,
# 	'702_S_Pen_20_W' => 0,
# 	'739_S_Qua_21_W' => 0,
# 	'769_S_Qua_22_W' => 0,
# 	'783_S_Qua_22_W' => 0,
# 	'784_S_Qua_22_W' => 0,
# 	'794_S_Qua_23_W' => 0,
# 	'804_S_Qua_23_W' => 0,
# 	'816_S_Qua_23_W' => 0,
# 	'859_S_Niw_25_W' => 0,
# 	'887_S_Niw_25_W' => 0,
# 	'888_S_Niw_25_W' => 0,
# 	'898_S_Niw_25_W' => 0,
# 	'949_S_Niw_27_W' => 0,
# 	);

# my %mela = (
# 	'SRR8700077_black' => 0,
# 	'SRR8700078_black' => 0,
# 	'SRR8700079_black' => 0,
# 	'SRR8700080_black' => 0,
# 	'SRR8700081_black' => 0,
# 	'SRR8700082_black' => 0,
# 	'SRR8700083_black' => 0,
# 	'SRR8700084_black' => 0,
# 	'SRR8700085_black' => 0,
# 	'SRR8700086_black' => 0,
# 	'SRR8700087_red' => 0,
# 	'SRR8700088_red' => 0,
# 	'SRR8700089_red' => 0,
# 	'SRR8700090_red' => 0,
# 	'SRR8700091_red' => 0,
# 	'SRR8700092_red' => 0,
# 	'SRR8700093_red' => 0,
# 	'SRR8700094_red' => 0,
# 	'SRR8700095_red' => 0,
# 	'SRR8700096_red' => 0,
# 	'SRR8700097_red' => 0,
# 	);

my %bifa = (
	'181_S_Bor_8_W' => 0,
	'334_S_Hor_11_W' => 0,
	'528_S_Eva_16_W' => 0,
	'530_S_Eva_16_W' => 0,
	'538_S_Eva_16_W' => 0,
	'555_S_Eva_16_W' => 0,
	'625_S_Pen_18_W' => 0,
	'626_S_Pen_18_W' => 0,
	'634_S_Pen_18_W' => 0,
	'682_S_Pen_19_W' => 0,
	'719_S_Qua_21_W' => 0,
	'721_S_Qua_21_W' => 0,
	'723_S_Qua_21_W' => 0,
	'730_S_Qua_21_W' => 0,
	'733_S_Qua_21_W' => 0,
	'742_S_Qua_21_W' => 0,
	'788_S_Qua_22_W' => 0,
	);




## Keep count of differences between each haplotype comparison
my %dxy;

#my $length = 252081862; # B. sylvicola genome size
#my $length = 3118388; # Length of contig_026 for testing
#print "Sequence length = $length\n";

# Set window size
my $window_size = 20000;
my $window_start = 0;
my $window_end = $window_size;
my $window_mid = $window_start + ($window_size/2);

my $gen_window_start = 0;
my $gen_window_end = $window_size;
my $gen_window_mid = $gen_window_start + ($window_size/2);

my $prev_contig = "contig_078";	#set to be name of first contig to be looked at

my $genes_in_this_window="false";



# Get the gene coords hash for the first contig and add up the length of genes overlapping with the first window, in order to subtract this length from the window length later

my %current_contig_hash = %{ $gene_coords{$prev_contig} };

my $total_window_genic_length = 0;

foreach my $key (keys %current_contig_hash) {
	my $coords = $current_contig_hash{$key};
	my @data = split(/-/, $coords);
	my $start = $data[0];
	my $end = $data[1];
	
	if (($start >= $window_start) && ($end <= $window_end)) { # case where gene is fully contained within window
	
		my $gene_length = $end-$start;
		$total_window_genic_length = $total_window_genic_length + $gene_length;
		$genes_in_this_window="true";
		
	}
	
	elsif (($start >= $window_start) && ($start <= $window_end) && ($end >= $window_end)) { # case where a gene starts within a window but ends beyond it
		
		my $gene_length = $window_end-$start;
		$total_window_genic_length = $total_window_genic_length + $gene_length;
		$genes_in_this_window="true";
		
	}
	
	elsif (($start <= $window_start) && ($end >= $window_start) && ($end <= $window_end)) { # case where a gene starts outside a window but ends within it
		
		my $gene_length = $end-$window_start;
		$total_window_genic_length = $total_window_genic_length + $gene_length;
		$genes_in_this_window="true";
		
	}
	
	elsif (($start <= $window_start) && ($end >= $window_end)) { # case where an entire window is part of a gene
		
		$total_window_genic_length = $window_size;
		$genes_in_this_window="true";
		
	}
	
}

#print "Done with gene coords hash\n";


while (my $line = <IN3>) {
	
	chomp $line;				

	if ( $line =~ m/^\##/ ) {				
		next;
	}
	
	if ( $line =~ m/^\#CHROM/ ) { 
		
		my @headerLine = split(/\t/,$line);				
		#@names = @headerLine[ 9..333 ];					##  assign the headers (names) to the names array
		
		next;
		
	}
	
	my @elements = split(/\t/,$line); 					
	
	my $contig = $elements[0];
	
	my @split_contig = split(/_/,$contig); #comment out when running whole file
	my $contig_number = $split_contig[1];	#comment out when running whole file
	
	if ($contig_number < 78) {	#comment out when running whole file
	
		next;	#comment out when running whole file
		
	}	#comment out when running whole file


	
	my $pos = $elements[1];
	
	my $skip_snp = "false";
	
	if ($genes_in_this_window eq "true") { # no need to enter foreach loop if we know there are no genes in the current window
		# Check if SNP is within a gene
	
		%current_contig_hash = %{ $gene_coords{$contig} };
	
		foreach my $key (keys %current_contig_hash) {
			my $coords = $current_contig_hash{$key};
			my @data = split(/-/, $coords);
			my $start = $data[0];
			my $end = $data[1];
		
			if (($pos >= $start) && ($pos <= $end)) {
			
				$skip_snp = "true";
			
			}
		}

	}
	
# 	if ($skip_snp eq "true") {
# 		
# 		print $skipped_snps_fh "$line\n";
# 			
# 	}
	
	my $gen_pos = $pos + $lengths{$contig};
	
	my $ref_base = $elements[3];
	my $alt_base = $elements[4];
	
	if ($contig eq $prev_contig) {
	
		if (($pos > $window_start) && ($pos <= $window_end)) {
		
			if ($skip_snp eq "false") {
	
				for (my $i=9; $i <= 330; $i++) { # loop through all (i) of the samples in the vcf...
					my $sampleA = $samples{$i};
		
					for (my $j=$i+1; $j <= 330; $j++) { # and compare them to all other (j) samples
						my $sampleB = $samples{$j};
				
						if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
							next;
						}
				
						if ((exists $bifa{$sampleA}) && (exists $bifa{$sampleB})) { # To ensure we are only comparing samples from the same species
				
							my $comp = "$sampleA\t$sampleB";
			
							my $data_i = $elements[ $i ];		
							my @i_data = split(/:/, $data_i);
							my $geno_i = $i_data[ 0 ]; # Get the genotype data for the ith sample
		
							my $data_j = $elements[ $j ];
							my @j_data = split(/:/, $data_j);
							my $geno_j = $j_data[ 0 ]; # get the genotype data for the jth sample
			
							if (($geno_i eq "./.") || ($geno_j eq "./.")) {
								next; # skip missing data
							}
			
							else {
			
								my @i_alleles = split(/\//,$geno_i);
								my @i_alleles_sorted = sort (@i_alleles); # This will ensure we are always comparing alleles in alphabetical order
								my $i_a1 = $i_alleles_sorted[0];
								my $i_a2 = $i_alleles_sorted[1];

	
								my @j_alleles = split(/\//,$geno_j);
								my @j_alleles_sorted = sort (@j_alleles);
								my $j_a1 = $j_alleles_sorted[0];
								my $j_a2 = $j_alleles_sorted[1];
	
								if ($i_a1 != $j_a1) { # Check whether the alleles match and if not add a count to the $hap1hap1 counter
									$dxy{$comp}++;
								}
		
								if ($i_a2 != $j_a2) {
									$dxy{$comp}++;
								}
							}	
						}
				
						else {
							next;
						}
					}	
				}
		
			}
			
		}
	
		else { # We are into a new window, deal with results from previous window before continuing
		
			#print "Done with that window\n";
		
			$window_mid = $window_start + ($window_size/2);
			$gen_window_start = $window_start + $lengths{$prev_contig};
			$gen_window_end = $gen_window_start + $window_size;
			$gen_window_mid = $gen_window_start + ($window_size/2);

			my $bifa_count = 0;
			my $bifa_sum = 0;
			
			my $adj_window_size = $window_size - $total_window_genic_length; # to adjust for the length of genes within the window

			if ($adj_window_size == 0) { # Entire window is genic!
			
				print $output_fh "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\tNA\n";
				
				while( my( $key, $value ) = each %dxy ){
					$dxy{$key} = 0; # Reset dxy values for next window
				}
				
			}
			
			else {
				while( my( $key, $value ) = each %dxy ){
			
   					my $dxy = $value/$adj_window_size *100;
   
   					#print OUT "$key\t$dxy\n";
   
   					my @samples = split(/\t/, $key);
  					my $name1 = $samples[0];
   					my $name2 = $samples[1];
   
   					$bifa_count++;
   					$bifa_sum = $bifa_sum + $dxy;
   	
					$dxy{$key} = 0; # Reset dxy values for next window
		
				}


				my $bifa_pi = $bifa_sum/$bifa_count;

				print $output_fh "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\t$bifa_pi\n";
			}
			
			
			$window_start = (int($pos/$window_size)) * $window_size;	# Set new window
			$window_end = $window_start + $window_size;	
			
			
			# Check gene content of next window
			$total_window_genic_length = 0;
			$genes_in_this_window="false";

			foreach my $key (keys %current_contig_hash) {
				my $coords = $current_contig_hash{$key};
				my @data = split(/-/, $coords);
				my $start = $data[0];
				my $end = $data[1];
	
				if (($start >= $window_start) && ($end <= $window_end)) { # case where gene is fully contained within window
	
					my $gene_length = $end-$start;
					$total_window_genic_length = $total_window_genic_length + $gene_length;
					$genes_in_this_window="true";
		
				}
	
				elsif (($start >= $window_start) && ($start <= $window_end) && ($end >= $window_end)) { # case where a gene starts within a window but ends beyond it
		
					my $gene_length = $window_end-$start;
					$total_window_genic_length = $total_window_genic_length + $gene_length;
					$genes_in_this_window="true";
		
				}
	
				elsif (($start <= $window_start) && ($end >= $window_start) && ($end <= $window_end)) { # case where a gene starts outside a window but ends within it
		
					my $gene_length = $end-$window_start;
					$total_window_genic_length = $total_window_genic_length + $gene_length;
					$genes_in_this_window="true";
		
				}
	
				elsif (($start <= $window_start) && ($end >= $window_end)) { # case where an entire window is part of a gene
		
					$total_window_genic_length = $window_size;
					$genes_in_this_window="true";
		
				}
	
			}	
		
		
			if ($skip_snp eq "false") {
		
				for (my $i=9; $i <= 330; $i++) {
					my $sampleA = $samples{$i};
		
					for (my $j=$i+1; $j <= 330; $j++) {
						my $sampleB = $samples{$j};
				
						if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
							next;
						}
					
						if ((exists $bifa{$sampleA}) && (exists $bifa{$sampleB})) {
				
							my $comp = "$sampleA\t$sampleB";
			
							my $data_i = $elements[ $i ];		
							my @i_data = split(/:/, $data_i);
							my $geno_i = $i_data[ 0 ]; # Get the genotype data for the ith sample
	
							my $data_j = $elements[ $j ];
							my @j_data = split(/:/, $data_j);
							my $geno_j = $j_data[ 0 ]; # get the genotype data for the jth sample
			
							if (($geno_i eq "./.") || ($geno_j eq "./.")) {
								next; # skip missing data
							}
			
							else {
			
								my @i_alleles = split(/\//,$geno_i);
								my @i_alleles_sorted = sort (@i_alleles); # This will ensure we are always comparing alleles in alphabetical order
								my $i_a1 = $i_alleles_sorted[0];
								my $i_a2 = $i_alleles_sorted[1];

	
								my @j_alleles = split(/\//,$geno_j);
								my @j_alleles_sorted = sort (@j_alleles);
								my $j_a1 = $j_alleles_sorted[0];
								my $j_a2 = $j_alleles_sorted[1];
	
								if ($i_a1 != $j_a1) { # Check whether the alleles match and if not add a count to the $hap1hap1 counter
									$dxy{$comp}++;
								}
	
								if ($i_a2 != $j_a2) {
									$dxy{$comp}++;
								}	
							}	
						}
				
						else {
							next;
						}	
					}
				}
			}		
		}
		
	}
	
	else { # We're on to the next contig, deal with previous window before progressing
	
		#print "Done with that contig\n";
		
		$window_mid = $window_start + ($window_size/2);
		$gen_window_start = $window_start + $lengths{$prev_contig};
		$gen_window_end = $gen_window_start + $window_size;
		$gen_window_mid = $gen_window_start + ($window_size/2);

		my $bifa_count = 0;
		my $bifa_sum = 0;
		
		my $adj_window_size = $window_size - $total_window_genic_length; # to adjust for the length of genes within the window

		if ($adj_window_size == 0) { # Entire window is genic!
			
			print $output_fh "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\tNA\n";
				
			while( my( $key, $value ) = each %dxy ){
				$dxy{$key} = 0; # Reset dxy values for next window
			}
				
		}
		
		
		else {
		
			while( my( $key, $value ) = each %dxy ){
   
   				my $dxy = $value/$adj_window_size *100;
   
   				#print OUT "$key\t$dxy\n";
   
   				my @samples = split(/\t/, $key);
  				my $name1 = $samples[0];
   				my $name2 = $samples[1];
   
   				$bifa_count++;
   				$bifa_sum = $bifa_sum + $dxy;
   	
				$dxy{$key} = 0; # Reset dxy values for next window
	
			}


			my $bifa_pi = $bifa_sum/$bifa_count;

			print $output_fh "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\t$bifa_pi\n";
		}
		
		$prev_contig = $contig;
		$window_start = (int($pos/$window_size)) * $window_size;	# Set new window
		$window_end = $window_start + $window_size;		
	
		$genes_in_this_window="false";

	# Get the gene coords hash for the next contig and add up the length of genes overlapping with the first window, in order to subtract this length from the window length later

		%current_contig_hash = %{ $gene_coords{$prev_contig} };

		$total_window_genic_length = 0;

		foreach my $key (keys %current_contig_hash) {
			my $coords = $current_contig_hash{$key};
			my @data = split(/-/, $coords);
			my $start = $data[0];
			my $end = $data[1];
	
			if (($start >= $window_start) && ($end <= $window_end)) { # case where gene is fully contained within window
	
				my $gene_length = $end-$start;
				$total_window_genic_length = $total_window_genic_length + $gene_length;
				$genes_in_this_window="true";
		
			}
	
			elsif (($start >= $window_start) && ($start <= $window_end) && ($end >= $window_end)) { # case where a gene starts within a window but ends beyond it
		
				my $gene_length = $window_end-$start;
				$total_window_genic_length = $total_window_genic_length + $gene_length;
				$genes_in_this_window="true";
		
			}
	
			elsif (($start <= $window_start) && ($end >= $window_start) && ($end <= $window_end)) { # case where a gene starts outside a window but ends within it
		
				my $gene_length = $end-$window_start;
				$total_window_genic_length = $total_window_genic_length + $gene_length;
				$genes_in_this_window="true";
		
			}
	
			elsif (($start <= $window_start) && ($end >= $window_end)) { # case where an entire window is part of a gene
		
				$total_window_genic_length = $window_size;
				$genes_in_this_window="true";
		
			}
	
		}
	
	
	
		if ($skip_snp eq "false") {
		
			for (my $i=9; $i <= 330; $i++) {
				my $sampleA = $samples{$i};
	
				for (my $j=$i+1; $j <= 330; $j++) {
					my $sampleB = $samples{$j};
			
					if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
						next;
					}
				
					if ((exists $bifa{$sampleA}) && (exists $bifa{$sampleB})) {
				
						my $comp = "$sampleA\t$sampleB";
		
						my $data_i = $elements[ $i ];		
						my @i_data = split(/:/, $data_i);
						my $geno_i = $i_data[ 0 ]; # Get the genotype data for the ith sample

						my $data_j = $elements[ $j ];
						my @j_data = split(/:/, $data_j);
						my $geno_j = $j_data[ 0 ]; # get the genotype data for the jth sample
		
						if (($geno_i eq "./.") || ($geno_j eq "./.")) {
							next; # skip missing data
						}
		
						else {
		
							my @i_alleles = split(/\//,$geno_i);
							my @i_alleles_sorted = sort (@i_alleles); # This will ensure we are always comparing alleles in alphabetical order
							my $i_a1 = $i_alleles_sorted[0];
							my $i_a2 = $i_alleles_sorted[1];


							my @j_alleles = split(/\//,$geno_j);
							my @j_alleles_sorted = sort (@j_alleles);
							my $j_a1 = $j_alleles_sorted[0];
							my $j_a2 = $j_alleles_sorted[1];

							if ($i_a1 != $j_a1) { # Check whether the alleles match and if not add a count to the $hap1hap1 counter
								$dxy{$comp}++;
							}

							if ($i_a2 != $j_a2) {
								$dxy{$comp}++;
							}
						}	
					}
			
					else {
						next;
					}	
				}
			}
		}
	}	
					
}
	
# And calculate Dxy and print for last window

$window_mid = $window_start + ($window_size/2);
$gen_window_start = $window_start + $lengths{$prev_contig};
$gen_window_end = $gen_window_start + $window_size;
$gen_window_mid = $gen_window_start + ($window_size/2);

my $bifa_count = 0;
my $bifa_sum = 0;

my $adj_window_size = $window_size - $total_window_genic_length; # to adjust for the length of genes within the window

if ($adj_window_size == 0) { # Entire window is genic!
			
	print $output_fh "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\tNA\n";
					
}

else {

	while( my( $key, $value ) = each %dxy ){
   
	   	my $dxy = $value/$adj_window_size *100;
   
	   	#print OUT "$key\t$dxy\n";
   
	   	my @samples = split(/\t/, $key);
	  	my $name1 = $samples[0];
	   	my $name2 = $samples[1];
   
	   	$bifa_count++;
 	  	$bifa_sum = $bifa_sum + $dxy;
   	
		$dxy{$key} = 0; # Reset dxy values for next window

	}


	my $bifa_pi = $bifa_sum/$bifa_count;

	print $output_fh "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\t$bifa_pi\n";

}



