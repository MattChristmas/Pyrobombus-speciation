#! /usr/bin/perl -w

#use strict; #enforce some good programming rules
use warnings;
use diagnostics;

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

open IN2, "sylv_genome_contig_cum_lengths.txt" or die;  #This file contains two columns, listing the contig names and cumulative lengths for the assembly

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

open OUT, ">>genmap_k150_e2_20kb_windows.txt" or die; 


print OUT "Contig\tWindow_start\tGEN_MID\tmappability\n";


open my $fh1, '<', "genmap_out.bed" or die; # This is the output file from GenMap

my $window_size=20000;
my $current_contig="contig_001";
my $window_number=0;
my $total_distance=0;
my $total_mappability=0;


while (my $line = <$fh1>) {

	chomp $line;
	
	my @data = split "\t", $line;
	
	my $contig = $data[0];
	my $position_1 = $data[1];
	my $position_2 = $data[2];
	my $mappability = $data[4];
	my $distance = $position_2-$position_1;
	my $summed_mappability = $mappability*$distance;
	
	
	if ($contig eq $current_contig) { # we are on the same contig
	
		if ((int($position_2/$window_size)) == $window_number) { # we are in the same window as before

			$total_distance = $total_distance + $distance;
			$total_mappability = $total_mappability + $summed_mappability;
			
		}
		
		else { # We are into another window, output previous window mappability score...
		
			# The first part of this region belongs in the previous window, first calculate how much:
			my $window_start = $window_number*$window_size;
			my $window_end = $window_start+$window_size;
			my $remainder = $window_end - $position_1;
			my $remainder_mappability = $remainder*$mappability;
			
			$total_distance = $total_distance + $remainder;
			$total_mappability = $total_mappability + $remainder_mappability;
			
			my $average_mappability = $total_mappability/$total_distance;
			
			#print "Window start = $window_start\n";
			my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
			print OUT "$current_contig\t$window_start\t$GEN_MID\t$average_mappability\n";
			
			# And then figure out which window we are in...
			
			my $end_window = int($position_2/$window_size);
			
			if ($end_window-$window_number==1) { # We are into the next immediate window
			
				$window_number++;
				$total_distance = $distance-$remainder;
				$total_mappability = $total_distance*$mappability;
			}
			
			else { # We have moved through more than one window, where all bases had the same mappability
			
				for (my $i = ($window_number+1); $i < $end_window; $i++) {
				
					$window_start = $i*$window_size;
					$GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
					
					print OUT "$current_contig\t$window_start\t$GEN_MID\t$mappability\n";
				
				}
				
				
				$window_number=$end_window;
				$total_distance = $position_2-($end_window*$window_size); # To get the remaining distance between the end of the last window and the
																				# end of this section
				$total_mappability = $total_distance*$mappability;
			
			}
			
		}
	
	}
	
	else { # We are on to a new contig, print out last part of previous contig first
			
		
		my $contig_length= $lengths{$current_contig};
		
		
		my $average_mappability = $total_mappability/$total_distance;
			
		my $window_start = $window_number*$window_size;
		my $GEN_MID = $cumul_lengths{$current_contig} + ($window_start) + ($window_size/2);
		print OUT "$current_contig\t$window_start\t$GEN_MID\t$average_mappability\n";
			
			
		$current_contig=$contig;
			
		# Need to check how many windows this first line counts for
		my $end_window = int($position_2/$window_size);
			
		if ($end_window==0) { # This first line is contained within the first 20kb window
			
			$window_number=0;
			$total_distance = $distance;
			$total_mappability = $summed_mappability;
		}
			
		else { # We have moved through more than one window, where all bases had the same mappability
			
			for (my $i = 0; $i < $end_window; $i++) {
				
				$window_start = $i*$window_size;
				$GEN_MID = $lengths{$current_contig} + ($window_start) + ($window_size/2);
					
				print OUT "$current_contig\t$window_start\t$GEN_MID\t$mappability\n";
				
			}
			
			
			$window_number=$end_window;
			$total_distance = $position_2-(($end_window-1)*$window_size); # To get the remaining distance between the end of the last window and the
																				# end of this section
			$total_mappability = $total_distance*$mappability;
		}
			
	
	}
	
}


#Print out final window	
my $average_mappability = $total_mappability/$total_distance;
			
my $window_start = $window_number*$window_size;
my $GEN_MID = $lengths{$current_contig} + ($window_start) + ($window_size/2);
print OUT "$current_contig\t$window_start\t$GEN_MID\t$average_mappability\n";
			

		
