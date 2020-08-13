#!/usr/bin/perl -w

use strict; #enforce some good programming rules
use warnings;
use diagnostics;
use Data::Dumper; 
use IO::Handle;

open IN1, "repeatmasker.gff" or die; # This is the output from repeatmasker

open OUT, ">>repeatmasker_single_records.gff" or die; # Here we'll store only single records for each repeat (i.e. remove surplus info from the file)

my $previous_start =0;
my $previous_stop=0;
my $previous_contig=0;

while (my $line = <IN1>) {
	
	chomp $line;				

	my @data = split(/\t/,$line);
	
	my $contig = $data[0];	
	
	my $start = $data[3];
	my $stop = $data[4];
	
	if ($contig eq $previous_contig) {
	
		if ($start == $previous_start) {
			next;
		}
		
		else {
		
		print OUT "$line\n";
		
		$previous_contig=$contig;
		$previous_start=$start;
		$previous_stop=$stop;
		
		}
	
	
	}
	
	else {
	
		print OUT "$line\n";
		
		$previous_contig=$contig;
		$previous_start=$start;
		$previous_stop=$stop;
	
	}

									
}