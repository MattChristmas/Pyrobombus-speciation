#!/usr/bin/perl -w

use strict; #enforce some good programming rules
use warnings;
use diagnostics;
use Data::Dumper; 
use IO::Handle;

open IN1, "bsyl_contig_lengths.txt" or die; # This file contains two columns, listing the contig names and lengths for the assembly

my %lengths;

while (my $line = <IN1>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);	
	
	my $contig = $data[0];
	my $length = $data[1];
	
	#print "$contig\t$length\n";
	
	$lengths{$contig} = $length;
	
# 	while( my( $key, $value ) = each %lengths ){
#     print "$key: $value\n";
# }
									
}

open IN2, "sylv_genome_contig_cum_lengths.txt" or die; #This file contains two columns, listing the contig names and cumulative lengths for the assembly

my %cumul_lengths;

while (my $line = <IN2>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);	
	
	my $contig = $data[0];
	my $length = $data[1];
	
	#print "$contig\t$length\n";
	
	$cumul_lengths{$contig} = $length;
	
# 	print Dumper(\%lengths);
									
}

open OUT, ">>Bsyl_whole_genome_LDhat_stat_out_bpen_1_20kb_windows.txt";
print OUT "Contig\tWindow_start\tGEN_MID\trho_per_kb\n";

open IN3, "Bsyl_whole_genome_LDhat_stat_out_bpen_1.txt" or die; # This is the output file from LDHat stat

my $header=<IN3>;

my $window_size=20; # in kbp
my $adj_window_size=20;
my $prev_position=1;
my $prev_rho=0;
my $window_rho=0;
my $window_distance=0;
my $current_contig="contig_000";
my $window_counter=0; 

while (my $line = <IN3>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);
	
	my $contig = $data[0];
	#print "$current_contig\n";
	my $pos = $data[1]; # value is in kbp
	my $mean_rho = $data[2];
	
	if ($contig eq $current_contig) {
	
		if (int($pos/$window_size) == $window_counter) { # we are still in the same window
			
			my $physical_distance = $pos-$prev_position;
			my $interval_rho = $physical_distance * $prev_rho;
			$window_rho=$window_rho+$interval_rho;
			$prev_position=$pos;
			$prev_rho=$mean_rho;
		
		}
		
		else { # We are into a new window, print out previous before continuing
		
			my $remainder = $window_size*($window_counter+1)-$prev_position;
			my $physical_distance = $pos-$prev_position;
			my $interval_rho = $physical_distance * $prev_rho;
			my $remainder_rho = $remainder/$physical_distance*$interval_rho; # gives rho for the remaining part of the previous window
			$window_rho = $window_rho + $remainder_rho;
				
			my $window_start= $window_counter*$window_size;
			my $window_rho_per_kbp = $window_rho/$adj_window_size;
			my $GEN_MID = ($window_start*1000+10000) + $cumul_lengths{$current_contig};
			
			
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$window_rho_per_kbp\n";
		
			
			# Adjust values of counters for current line
			$window_counter=$window_counter+1;
			my $start_length = $pos-($window_size*$window_counter);
			$window_rho = $start_length/$physical_distance*$interval_rho;
			$adj_window_size=$window_size;
			$prev_position=$pos;
			$prev_rho=$mean_rho;
				
		
		}
	
	}
	
	else { 
	
		
	
		if ($current_contig eq "contig_000") { # Then we are looking at the first entry, set starting values
		
			$current_contig= $contig;
			$prev_position = $pos;
			$adj_window_size = 20-$pos; #To take account of base pairs at start of contig where we have no info on rho
			$prev_rho=$mean_rho; # Don't do anything with this until we get to the next line
			$window_counter=0;
			
		
		}
		
		else { # We are into a new contig, need to take account of end of last contig, and start of this one
	
			my $window_start = $window_counter*$window_size;
			my $window_end = $lengths{$current_contig};
			my $window_length = $window_end-$window_start;
			my $window_rho_per_kbp = $window_rho/$window_length*1000;
			my $GEN_MID = ($window_start*1000+10000) + $cumul_lengths{$current_contig};
			
			
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$window_rho_per_kbp\n";
			
			# Adjust values of counters for current line
			$current_contig= $contig;
			$prev_position = $pos;
			$adj_window_size = 20-$pos; #To take account of base pairs at start of contig where we have no info on rho
			$prev_rho=$mean_rho; # Don't do anything with this until we get to the next line
			$window_rho=0;
			$window_counter=0;
		
		}
	
	}
	
}
	
	
	