#! /usr/bin/perl -w

#use strict; #enforce some good programming rules
use warnings;
use diagnostics;

open OUT1, ">>average_dxy_from_IoD_midpoints_5Mbp.txt" or die;
#open OUT2, ">>windows_within_1mb_downstream_pbs_peaks.txt" or die;

#print OUT1 "Distance\tdxy\tspecies_comparison\n";


open my $fh1, '<', "BsBi_large_IoDs_midpoints.txt" or die; # This file contains two columns, 1st is the chromosome, 2nd is the IoD midpoint position
open my $fh2, '<', "dxy_genome_coords.txt" or die; # This is a file containing dxy values in 20kbp windows across the genome

my $header = <$fh2>;
my @dxy_data = <$fh2>;


my $window_size=20000;
my $max_distance = 5000000; # This is the distance we want to move away from the dxy peak, both upstream and downstream

my %upstream_hash; # For storing the dxy values upstream of the peak...
my %downstream_hash; # ...and downstream
my %upstream_counter; # for keeping count of the number of observations made at each distance in order to calculate averages at the end
my %downstream_counter;

while (my $line = <$fh1>) {

	
	chomp $line;
	
	my @data = split "\t", $line;
	
	my $peak_contig = $data[0];
	my $peak_start = $data[1];

	
	for (my $i=$peak_start; $i <= $peak_start+$max_distance; $i=$i+$window_size) { # To move from 0 to 5 mb away from the peak downstream, 20 kbp at a time
	
		#print "$i\n";
	
		for my $window (@dxy_data) {
			
			chomp $window;
	
			my @elements = split "\t", $window;
			my $contig = $elements[0];
			my $start = $elements[1];
			my $dxy = $elements[3]; # CHANGE FOR EACH SPECIES COMPARISON
		
			next if ($peak_contig ne $contig);
			
			if (($start >= ($i - $window_size/2)) && ($start <= ($i + $window_size/2))) {
		
				$current_i = $i-$peak_start; # To get the distance from the peak in bp
				
				if (exists ($downstream_hash{$current_i})) {
				
					$downstream_hash{$current_i} = $downstream_hash{$current_i} + $dxy; # to save the dxy of this window to the downstream hash
					$downstream_counter{$current_i}++;
				
				}
				
				else {
				
					$downstream_hash{$current_i} = $dxy; # to save the dxy of this window to the downstream hash
					$downstream_counter{$current_i}=1;
				
				}
			}
			
		}
		
	}
	
	for (my $j=$peak_start; $j >= $peak_start-$max_distance; $j=$j-$window_size) { # Then do the same as above but moving upstream from the peak
	
		#print "$i\n";
	
		for my $window (@dxy_data) {
			
			chomp $window;
	
			my @elements = split "\t", $window;
			my $contig = $elements[0];
			my $start = $elements[1];
			my $dxy = $elements[3]; # CHANGE FOR EACH SPECIES COMPARISON
		
			next if ($peak_contig ne $contig);
			
			if (($start >= ($j - $window_size/2)) && ($start <= ($j + $window_size/2))) {
		
				$current_j = $peak_start-$j;
				
				if (exists ($upstream_hash{$current_j})) {
				
					$upstream_hash{$current_j} = $upstream_hash{$current_j} + $dxy; # to save the dxy of this window to the upstream hash
					$upstream_counter{$current_j}++;
				
				}
				
				else {
				
					$upstream_hash{$current_j} = $dxy; # to save the dxy of this window to the upstream hash
					$upstream_counter{$current_j}=1;
				
				}
				
			}
			
		}
		
	}
}
	
for (my $x=0; $x <= $max_distance; $x=$x+$window_size) { # Now to get the average dxys per distance from peak
	
	if ((exists $downstream_hash{$x}) && (exists $upstream_hash{$x})) { # If there is an upstream and a downstream dxy value at $x distance from the peak...
	
		my $total_obs = $upstream_counter{$x} + $downstream_counter{$x};
		my $total_dxy = $downstream_hash{$x} + $upstream_hash{$x};
		my $av_dxy = $total_dxy/$total_obs;
		print OUT1 "$x\t$av_dxy\tBs_Bi\n"; # ...and print
			
	}
		
	elsif (exists $downstream_hash{$x}) { # else, if we only have a downstream value at that distance (perhaps this peak is close to the start of a contig so not much upstream of it)...
		
		my $total_obs = $downstream_counter{$x};
		my $total_dxy = $downstream_hash{$x};
		my $av_dxy = $total_dxy/$total_obs;
		print OUT1 "$x\t$av_dxy\tBs_Bi\n"; # ...and print
		
	}
		
	elsif (exists $upstream_hash{$x}) { # else, if we only have an upstream value at that distance...
		
		my $total_obs = $upstream_counter{$x};
		my $total_dxy = $upstream_hash{$x};
		my $av_dxy = $total_dxy/$total_obs;
		print OUT1 "$x\t$av_dxy\tBs_Bi\n"; # ...and print
		
	}

	
}

			

		
