#!/usr/bin/perl -w

use strict; #enforce some good programming rules
use warnings;
use diagnostics;
use Data::Dumper; 

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

open IN2, "bifa_vanco_whole_genome_biallelic_snps_haploid_hets_removed.vcf" or die; # This is the vcf file containing the SNP data

open OUT1, ">>bifa_vanco_dxy_20kb_windows.txt" or die; # Write out to here

print OUT1 "contig\tcontig_start\tcontig_mid\tcontig_end\tGEN_START\tGEN_MID\tGEN_END\tdxy\n";


my %samples = ( # Set up a sample array where the key is the position of the sample in the vcf file. Names match those in the vcf file
	9 => '181_S_Bor_8_W',
	10 => '334_S_Hor_11_W',
	11 => '528_S_Eva_16_W',
	12 => '530_S_Eva_16_W',
	13 => '538_S_Eva_16_W',
	14 => '555_S_Eva_16_W',
	15 => '625_S_Pen_18_W',
	16 => '626_S_Pen_18_W',
	17 => '634_S_Pen_18_W',
	18 => '682_S_Pen_19_W',
	19 => '719_S_Qua_21_W',
	20 => '721_S_Qua_21_W',
	21 => '723_S_Qua_21_W',
	22 => '730_S_Qua_21_W',
	23 => '733_S_Qua_21_W',
	24 => '788_S_Qua_22_W',
	25 => 'SRR10590579_B.v.nearcticus',
	26 => 'SRR10590580_B.v.nearcticus',
	27 => 'SRR10590581_B.v.nearcticus',
	28 => 'SRR10590582_B.v.nearcticus',
	29 => 'SRR10590583_B.v.nearcticus',
	30 => 'SRR10590584_B.v.nearcticus',
	31 => 'SRR10590586_B.v.vancouverensis',
	32 => 'SRR10590587_B.bifarius',
	33 => 'SRR10590588_B.v.nearcticus',
	34 => 'SRR10590589_B.v.nearcticus',
	35 => 'SRR10590590_B.v.nearcticus',
	36 => 'SRR10590591_B.v.nearcticus',
	37 => 'SRR10590592_B.v.nearcticus',
	38 => 'SRR10590593_B.v.nearcticus',
	39 => 'SRR10590594_B.v.nearcticus',
	40 => 'SRR10590595_B.v.nearcticus',
	41 => 'SRR10590596_B.v.nearcticus',
	42 => 'SRR10590597_B.v.nearcticus',
	43 => 'SRR10590599_B.bifarius',
	);

my %bifa = ( # Define which samples belong to which species/population
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
	'SRR10590585_B.bifarius' => 0,
	'SRR10590587_B.bifarius' => 0,
	'SRR10590599_B.bifarius' => 0,
	);

my %vanco = (
	'SRR10590579_B.v.nearcticus' => 0,
	'SRR10590580_B.v.nearcticus' => 0,
	'SRR10590581_B.v.nearcticus' => 0,
	'SRR10590582_B.v.nearcticus' => 0,
	'SRR10590583_B.v.nearcticus' => 0,
	'SRR10590584_B.v.nearcticus' => 0,
	'SRR10590588_B.v.nearcticus' => 0,
	'SRR10590589_B.v.nearcticus' => 0,
	'SRR10590590_B.v.nearcticus' => 0,
	'SRR10590591_B.v.nearcticus' => 0,
	'SRR10590592_B.v.nearcticus' => 0,
	'SRR10590593_B.v.nearcticus' => 0,
	'SRR10590594_B.v.nearcticus' => 0,
	'SRR10590595_B.v.nearcticus' => 0,
	'SRR10590596_B.v.nearcticus' => 0,
	'SRR10590597_B.v.nearcticus' => 0,
	);


## Keep count of differences between each haplotype comparison
my %dxy;

#my $length = 252081862; # B. sylvicola genome size
#print "Sequence length = $length\n";

# Set window size
my $window_size = 20000;
my $window_start = 0;
my $window_end = $window_size;
my $window_mid = $window_start + ($window_size/2);

my $gen_window_start = 0;
my $gen_window_end = $window_size;
my $gen_window_mid = $gen_window_start + ($window_size/2);

my $prev_contig = "contig_001";


while (my $line = <IN2>) {
	
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
# 	my @split_contig = split(/_/,$contig); #comment out when running whole file
# 	my $contig_number = $split_contig[1];	#comment out when running whole file
# 	
# 	if ($contig_number > 25) {	#comment out when running whole file
# 	
# 		last;	#comment out when running whole file
# 		
# 	}	#comment out when running whole file
	
	
	my $pos = $elements[1];
	my $gen_pos = $pos + $lengths{$contig};
	
	my $ref_base = $elements[3];
	my $alt_base = $elements[4];
	
	if ($contig eq $prev_contig) {
	
		if (($pos > $window_start) && ($pos <= $window_end)) {
	
			for (my $i=9; $i <= 43; $i++) { # loop through all (i) of the samples in the vcf...
				my $sampleA = $samples{$i};
		
				for (my $j=$i+1; $j <= 43; $j++) { # and compare them to all other (j) samples
					my $sampleB = $samples{$j};
				
					if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
						next;
					}
				
					if ((exists $bifa{$sampleA}) && (exists $vanco{$sampleB})) { # To ensure we are only comparing samples between the two groups
				
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
	
		else { # We are into a new window, deal with results from previous window before continuing
		
			$window_mid = $window_start + ($window_size/2);
			$gen_window_start = $window_start + $lengths{$prev_contig};
			$gen_window_end = $gen_window_start + $window_size;
			$gen_window_mid = $gen_window_start + ($window_size/2);

			my $bifa_vanco_count = 0;
			my $bifa_vanco_sum = 0;
	

			while( my( $key, $value ) = each %dxy ){
   
   				my $dxy = $value/$window_size;
   
   				#print OUT "$key\t$dxy\n";
   
   				my @samples = split(/\t/, $key);
  				my $name1 = $samples[0];
   				my $name2 = $samples[1];
   
   				$bifa_vanco_count++;
   				$bifa_vanco_sum = $bifa_vanco_sum + $dxy;
   	
				$dxy{$key} = 0; # Reset dxy values for next window
		
			}


			my $bifa_vanco_dxy = $bifa_vanco_sum/$bifa_vanco_count;

			print OUT1 "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\t$bifa_vanco_dxy\n";
		
			$window_start = (int($pos/$window_size)) * $window_size;	# Set new window
			$window_end = $window_start + $window_size;		
		
		
			for (my $i=9; $i <= 43; $i++) {
				my $sampleA = $samples{$i};
		
				for (my $j=$i+1; $j <= 43; $j++) {
					my $sampleB = $samples{$j};
				
					if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
						next;
					}
					
					if ((exists $bifa{$sampleA}) && (exists $vanco{$sampleB})) {
				
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
	
	else { # We're on to the next contig, deal with previous window before progressing
		
		$window_mid = $window_start + ($window_size/2);
		$gen_window_start = $window_start + $lengths{$prev_contig};
		$gen_window_end = $gen_window_start + $window_size;
		$gen_window_mid = $gen_window_start + ($window_size/2);

		my $bifa_vanco_count = 0;
		my $bifa_vanco_sum = 0;


		while( my( $key, $value ) = each %dxy ){
   
   			my $dxy = $value/$window_size;
   
   			#print OUT "$key\t$dxy\n";
   
   			my @samples = split(/\t/, $key);
  			my $name1 = $samples[0];
   			my $name2 = $samples[1];
   
   			$bifa_vanco_count++;
   			$bifa_vanco_sum = $bifa_vanco_sum + $dxy;
   	
			$dxy{$key} = 0; # Reset dxy values for next window
	
		}


		my $bifa_vanco_dxy = $bifa_vanco_sum/$bifa_vanco_count;

		print OUT1 "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\t$bifa_vanco_dxy\n";
	
		$prev_contig = $contig;
		$window_start = (int($pos/$window_size)) * $window_size;	# Set new window
		$window_end = $window_start + $window_size;		
	
	
		for (my $i=9; $i <= 43; $i++) {
			my $sampleA = $samples{$i};
	
			for (my $j=$i+1; $j <= 43; $j++) {
				my $sampleB = $samples{$j};
			
				if ($sampleA eq $sampleB) { # avoids comparing samples to themselves
					next;
				}
				
				if ((exists $bifa{$sampleA}) && (exists $vanco{$sampleB})) {
			
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
	
# And calculate Dxy and print for last window

$window_mid = $window_start + ($window_size/2);
$gen_window_start = $window_start + $lengths{$prev_contig};
$gen_window_end = $gen_window_start + $window_size;
$gen_window_mid = $gen_window_start + ($window_size/2);

my $bifa_vanco_count = 0;
my $bifa_vanco_sum = 0;


while( my( $key, $value ) = each %dxy ){
   
   	my $dxy = $value/$window_size;
   
   	#print OUT "$key\t$dxy\n";
   
   	my @samples = split(/\t/, $key);
  	my $name1 = $samples[0];
   	my $name2 = $samples[1];
   
   	$bifa_vanco_count++;
   	$bifa_vanco_sum = $bifa_vanco_sum + $dxy;
   	
	$dxy{$key} = 0; # Reset dxy values for next window

}


my $bifa_vanco_dxy = $bifa_vanco_sum/$bifa_vanco_count;

print OUT1 "$prev_contig\t$window_start\t$window_mid\t$window_end\t$gen_window_start\t$gen_window_mid\t$gen_window_end\t$bifa_vanco_dxy\n";




