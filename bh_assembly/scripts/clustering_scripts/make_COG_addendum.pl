#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by run_clustering.pl via cogtriangle_run_C.sh | Submits nothing       #
# Adds entries for every protein onto the initial clustering output               #
###################################################################################

open IN, $ARGV[0];

my %CDS;

foreach (<IN>) {
	chomp;
	my @data = split(/\,/,$_);
	$CDS{$data[0]} = $data[1];
}

close IN;

open ADD, "> addendum.out";

open CLU, $ARGV[1];

foreach (<CLU>) {
	chomp;
	my @data = split(/\,/,$_);
	$CDS{$data[0]} = 0;
}

close CLU;

foreach my $cds (sort keys %CDS) {
	if (defined($CDS{$cds}) && $CDS{$cds} ne '0') {
		print ADD "$cds,$CDS{$cds},$cds,,,,,\n";
	}
}

print STDERR "Completed make_COG_addendum!\n";
