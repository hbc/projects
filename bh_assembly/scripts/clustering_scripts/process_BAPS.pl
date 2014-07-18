#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by orthologue_alignment_processing.pl | Submits nothing               #
# Runs BAPS on the core genome alignment and processes the results into a Figtree #
# readable annotation file                                                        #
###################################################################################

open INFO, $ARGV[0] or die;

my %link;

foreach (<INFO>) {
	chomp;
	my @data = split(/\s+/,$_);
	$link{$data[0]} = $data[1];
}

close INFO;

open BAPS, $ARGV[1] or die;

my $in_cluster = 0;
my %cluster;
my $n;

foreach (<BAPS>) {
	chomp;
	if (/^Cluster /) {
		$in_cluster = 1;
		my @data = split(/\s+/,$_);
		shift(@data);
		$n = shift(@data);
		$n =~ s/\://g;
		foreach my $entry (@data) {
			$entry =~ s/\,|\{|\}//g;
			if (defined($entry) && length($entry) > 0) {
				push(@{$cluster{$n}},$entry);
			}
		}
	}
	if (/\}/) {
		unless (/\{/) {
			my @data = split(/\s+/,$_);
			foreach my $entry (@data) {
				$entry =~ s/\,|\}|\{//g;
				if (defined($entry) && length($entry) > 0) {
					push(@{$cluster{$n}},$entry);
				}
			}
		}
		$in_cluster = 0;
	}
	if ($in_cluster == 1) {
		unless (/\{/) {
			my @data = split(/\s+/,$_);
			foreach my $entry (@data) {
				$entry =~ s/\,|\}|\{//g;
				if (defined($entry) && length($entry) > 0) {
					push(@{$cluster{$n}},$entry);
				}
			}
		}
	}
}

open BOUT, "> baps_clustering_annotation.tab";

print BOUT "Strain\tBAPS\n";

foreach $n (sort keys %cluster) {
	foreach my $entry (@{$cluster{$n}}) {
		print BOUT "$link{$entry}\t$n\n";
	}
}

print STDERR "Completed process_BAPs.pl\n";
