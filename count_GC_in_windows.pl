#! /usr/bin/perl -w

use strict; #enforce some good programming rules
use warnings;
use diagnostics;

# open OUT1, ">>Bsyl_contig_single_lines.fasta" or die; # If uncommented, this section will convert a genome.fasta to place sequence for each contig/chromosome onto single lines
# open IN1, "Bombus.sylvicola.genome.assembly.v1.fasta";
# 
# my $first_contig = <IN1>;
# print OUT1 "$first_contig";
# 
# while (my $line = <IN1>) {
#  
# 	chomp $line;
# 
# 	if ($line =~ m/\>/) {
# 		print OUT1 "\n$line\n";
# 	}
# 	
# 	else {
# 		print OUT1 $line;
# 	}
# }
# 
# close IN1;

open IN2, "sylv_genome_contig_cum_lengths.txt" or die; # This is a file containing the cumulative lengths of all the contigs in the assembly

my %lengths;

while (my $line = <IN2>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);	
	
	my $contig = $data[0];
	my $length = $data[1];
	
	#print "$contig\t$length\n";
	
	$lengths{$contig} = $length;
	
# 	print Dumper(\%lengths);
									
}

close IN2;

open OUT2, ">>Bsyl_gc_content_20kb_windows.txt" or die;

print OUT2 "Contig\tWindow\tcontig_start\tcontig_mid\tgenome_start\tgenome_mid\tGC_prop\n";

open IN3, "Bsyl_contig_single_lines.fasta" or die;

my $window_size = 20000;
my $length;
my $contig_name = "undef";
my $window_number = 1;


my %GC_hash; #set up a hash where we will store the GC proportions for each window

while (my $line = <IN3>) {
 
	chomp $line;
	
	
	if ($line =~ m/\>/) {
		
		if ($contig_name eq "undef") { #first line of the file this will be true
			$contig_name = $line;
			$window_number = 1;
		}
		
		else {
		
			#print "$contig_name\n"; #for testing
		
			foreach my $key (sort {$a<=>$b} keys %GC_hash) {
			
				my $contig_start = ($key - 1) * $window_size;
				my $contig_mid = $contig_start + ($window_size/2);
				
				my $genome_start = $contig_start + $lengths{$contig_name};
				my $genome_mid = $genome_start + ($window_size/2);
				
   				print OUT2 "$contig_name\t$key\t$contig_start\t$contig_mid\t$genome_start\t$genome_mid\t$GC_hash{$key}\n";
		}
		
		$contig_name = $line;
		$window_number = 1;
		%GC_hash = (); # Reset the hash for the next contig
		}
	}
	else {
	
		$length = length($line);
		
		#print "$length\n"; #for testing
	
		while (my $window = substr( $line, 0, $window_size )) {
			my $gc_count = 0;
			my $base_count = length($window); # This is to account for the fact the last window won't be exactly the specified window size when calculating proportion of GC
			my @characters = split //, $window; # this will split each of the characters (nucleotides)

			foreach my $character (@characters) {
				if ($character =~ /g|c/i) { # The 'i' match modifier makes it case insensitive
					$gc_count++;
				}
			}
			
		my $GC_prop = $gc_count/$base_count;
		
		$GC_hash{$window_number} = $GC_prop; 

		substr( $line, 0, $window_size ) = ""; #This removes the first 'x' number of bases from the string so next time through the loop we're looking at the next window
		
		$window_number++;
		}
	}
	
}

