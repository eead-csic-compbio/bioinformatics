#!/usr/bin/perl

# Reads output from bedtools getfasta a produces a 1-based TSV file matching old to new coordinates
# by parsing the FASTA header (which are 0-based). Note it prints strand as well. 

use strict;
use warnings;

my $getfasta_input = $ARGV[0] || die "# usage: $0 <get_fasta output>\n";

my ($orig,$MorexV3,$oldchr,$oldsta,$oldend,$chr,$sta,$end);

open(FA,"<",$getfasta_input) || die "# ERROR: cannot open $getfasta_input\n";
while(<FA>) {
	#>T::3:127834508-127834609::chr3H_LR890098.1:109311845-109311846
	#T
	if(/^>(\w)::([^:]+):(\d+)-(\d+)::([^:]+):(\d+)-(\d+)/) {
		($orig,$oldchr,$oldsta,$oldend,$chr,$sta,$end) = ($1,$2,$3,$4,$5,$6,$7);
		#edit chr name to match INSDC contig/scaffold name
		$chr =~ s/chr\dH_//;
	} else {
		$MorexV3 = substr($_,0,1);	

		if($orig eq $MorexV3) {
			printf("%s\t%d\t%s\t%s\t%d\t%s\t%s\n",
				$oldchr,$oldend,'+',$chr,$end,$orig,$MorexV3);
		} else {
			if(($orig eq 'A' && $MorexV3 eq 'T') ||
				($orig eq 'T' && $MorexV3 eq 'A') ||
				($orig eq 'C' && $MorexV3 eq 'G') ||
				($orig eq 'G' && $MorexV3 eq 'C')) {
				
				printf("%s\t%d\t%s\t%s\t%d\t%s\t%s\n",
                                	$oldchr,$oldend,'-',$chr,$end,$orig,$MorexV3);
			}
		}
	}
}
close(FA);
