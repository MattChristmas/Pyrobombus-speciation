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



open OUT, ">>repeat_content_20kb_windows.txt" or die;

print OUT "Contig\tWindow_start\tGEN_MID\tprop_repeat\n";


open IN3, "repeatmasker_single_records.gff" or die; # This is a file containing the output from 'parse_repeatmasker_out.pl'

my $window_size=20000;
my $current_contig="contig_447"; # set to be first contig in in file
my $window_number=0;
my $total_repeat_length=0;


while (my $line = <IN3>) {

	chomp $line;
	
	my @data = split "\t", $line;
	
	my $contig = $data[0];
	my $position_1 = $data[3];
	my $position_2 = $data[4];
	my $length = abs($position_2-$position_1);
	
	if ($contig eq $current_contig) { # we are on the same contig
	
		my $start; # this is in case position 2 is less than position 1
		my $stop;
	
		if ($position_2 > $position_1) {
			$start=$position_1;
			$stop=$position_2;
		}
		else {
			$start=$position_2;
			$stop=$position_1;
		}
		
		if ((int($start/$window_size) == $window_number) && (int($stop/$window_size) == $window_number)) { # we are entirely in the same window
		
			#print "Adding up repeat length for $current_contig window $window_number\n";
				
			$total_repeat_length=$total_repeat_length+$length;
			
		}
			
		elsif (int($start/$window_size) == $window_number) { # the repeat starts in the current window but continues into the next window (and maybe beyond), need to find out how far it goes...
			
			# First deal with previous window
				
			my $end_window = int($stop/$window_size);
				
			my $window_start = $window_number*$window_size;
			my $window_end = $window_start+$window_size;
			my $remainder = $window_end - $start;
			$total_repeat_length = $total_repeat_length + $remainder;
			my $prop_repeat = $total_repeat_length/$window_size;
			my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
			
			if ($end_window-$window_number==1) { # We are into the next immediate window
			
				$window_number++;
				$total_repeat_length = $length-$remainder;
			}
				
			else { # We have moved through more than one window, need to print for those windows
			
				for (my $i = ($window_number+1); $i < $end_window; $i++) {
				
					$window_start = $i*$window_size;
					$GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
					print OUT "$current_contig\t$window_start\t$GEN_MID\t1\n"; # All these windows are entirely repeats so proportion=1
				
				}
				
				
				$window_number=$end_window;
				$total_repeat_length = $stop-($end_window*$window_size); # To get the remaining distance between the end of the last window and the
																				# end of this section
					
			}
			

		}
		
		else { # The next repeat is in another window entirely, deal with last window then find out where we are...
		
			my $window_start = $window_number*$window_size;
			my $prop_repeat = $total_repeat_length/$window_size;
			my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
			
			$window_number++;
			
			if ((int($start/$window_size) == $window_number) && (int($stop/$window_size) == $window_number)) { # the next repeat is contained within the next window
				
				$total_repeat_length=$length;
			}
			
			elsif (int($start/$window_size) == $window_number) { # the repeat starts in the current window but continues into the next window (and maybe beyond), need to find out how far it goes...
			
				# First deal with previous window
				
				my $end_window = int($stop/$window_size);
				
				my $window_start = $window_number*$window_size;
				my $window_end = $window_start+$window_size;
				my $remainder = $window_end - $start;
				$total_repeat_length = $remainder;
				my $prop_repeat = $total_repeat_length/$window_size;
				my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
				print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
			
				if ($end_window-$window_number==1) { # We are into the next immediate window
			
					$window_number++;
					$total_repeat_length = $length-$remainder;
				}
				
				else { # We have moved through more than one window, need to print for those windows
			
					for (my $i = ($window_number+1); $i < $end_window; $i++) {
				
					$window_start = $i*$window_size;
					$GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
					print OUT "$current_contig\t$window_start\t$GEN_MID\t1\n"; # All these windows are entirely repeats so proportion=1
				
					}
				
				
					$window_number=$end_window;
					$total_repeat_length = $stop-($end_window*$window_size); # To get the remaining distance between the end of the last window and the
																				# end of this section
					
				}
			}
			
			else { # We have moved through (multiple) windows where there are no repeats, need to print those
			
				my $start_window = int($start/$window_size);
				
				for (my $i = ($window_number); $i < $start_window; $i++) {
				
					my $window_start = $i*$window_size;
					$GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
					print OUT "$current_contig\t$window_start\t$GEN_MID\t0\n"; # All these windows do not contain repeats so proportion=0
				
				}
				
				$window_number=$start_window;
				$total_repeat_length = $length;
				
			}	

		}	
			
	}
			
					
	else { # We are onto a new contig, print last window of previous contig before resetting for new one
	
		my $contig_length = $lengths{$current_contig};	
		my $window_start = $window_number*$window_size;
		my $window_end = $window_start+$window_size;
		my $end_window = int($contig_length/$window_size);
		
		if ($contig_length>$window_end) { # The contig continues beyond the end of this window
			
			my $prop_repeat = $total_repeat_length/$window_size;
			my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
			
			for (my $i = ($window_number+1); $i < $end_window; $i++) {
				
				$window_start = $i*$window_size;
				$GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
				print OUT "$current_contig\t$window_start\t$GEN_MID\t0\n"; # All these windows do not contain repeats so proportion=0
				
			}
			
		}
		
		else { # The contig stops short of the end of the window, need to adjust the length
		
			my $adjusted_window_size = $contig_length - $window_start;
			my $prop_repeat = $total_repeat_length/$adjusted_window_size;
			my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
			
			if ($prop_repeat > 1) { # Not sure why this is sometimes true at the end of contigs...
			
				$prop_repeat=1; 
			
			}
					
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
		
		}
		
		$current_contig=$contig;
		$window_number=0;
		$total_repeat_length=$length;
		
		
	}
			
}		

# Print out last window

my $contig_length = $lengths{$current_contig};	
my $window_start = $window_number*$window_size;
my $window_end = $window_start+$window_size;
my $end_window = int($contig_length/$window_size);
		
if ($contig_length>$window_end) { # The contig continues beyond the end of this window
			
	my $prop_repeat = $total_repeat_length/$window_size;
	my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
	print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
			
	for (my $i = ($window_number+1); $i < $end_window; $i++) {
				
		$window_start = $i*$window_size;
		$GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
		print OUT "$current_contig\t$window_start\t$GEN_MID\t0\n"; # All these windows do not contain repeats so proportion=0
				
	}
			
}
		
else { # The contig stops short of the end of the window, need to adjust the length
		
	my $adjusted_window_size = $contig_length - $window_start;
	my $prop_repeat = $total_repeat_length/$adjusted_window_size;
	my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
	print OUT "$current_contig\t$window_start\t$GEN_MID\t$prop_repeat\n";
		
}
			
			
	