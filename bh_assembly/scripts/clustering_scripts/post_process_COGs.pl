#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by clustering_checkpoint.pl | Submits nothing                         #
# Looks for proteins that are not clustered by COGtriangles and assigns them into #
# pairs, singletons or assigns that as unclusterable based on self-BLAT           #
###################################################################################

# set up script_dir variable

my $script_dir = $ENV{SCRIPT_DIR};

# note reference sequence

my $refnum= $ARGV[3];
my $refstring = "seq"."$refnum"."orf";

# use original query CSV file to find all genes, as some won't appear in clustering or BLAT output

open CSV, $ARGV[2] or die;

my %unclustered;

while (my $line = <CSV>) {
	chomp $line;
	my @data = split(/\,/,$line);
	$unclustered{$data[0]} = 1;
}

close CSV;

# process COGcognitor output file to output CDSs assigned to multiple clusters or none

open IN, $ARGV[0] or die;

my %COG_hit;
my $blank = "X";

print STDERR "Reading in COGcognitor output...";

while (my $line = <IN>) {
	chomp $line;
	my @data = split(/\,/,$line);
	if (defined($data[5]) && length($data[5]) > 0 && $data[5] ne '-1') {
		$COG_hit{$data[0]}{$data[5]} = $data[4];
	} else {
		$COG_hit{$data[0]}{$blank} = $data[4];
	}
}

close IN;

print STDERR "done\nNow removing duplicates...";

# print out CDSs correctly assigned to one cluster, find best scoring hit for those
# assigned to multiple clusters

open OUT, "> all.strains.filtered.cls.csv";
open REF, "> reference.strain.filtered.cls.csv";

foreach my $CDS (sort keys %COG_hit) {
	if (scalar(keys %{$COG_hit{$CDS}}) == 1) {
		foreach my $clu (keys %{$COG_hit{$CDS}}) {
			unless ($clu eq $blank) {
				if ($CDS =~ /$refstring/) {
					print REF "$CDS,$clu\n";
				} else {
					print OUT "$CDS,$clu\n";
				}
				# added in to check on problem
				if (defined($unclustered{$CDS})) {
					delete($unclustered{$CDS});
				}
			}
		}
	} else {
		my $best_hit = 0;
		my $best_clu =$blank;
		foreach my $clu (keys %{$COG_hit{$CDS}}) {
			if ($COG_hit{$CDS}{$clu} > $best_hit) {
				$best_hit = $COG_hit{$CDS}{$clu};
				$best_clu = $clu;
			}
		}
		unless ($best_clu eq $blank) {
			if ($CDS =~ /$refstring/) {
				print REF "$CDS,$best_clu\n";
			} else {
				print OUT "$CDS,$best_clu\n";
			}
			if (defined($unclustered{$CDS})) {
				delete($unclustered{$CDS});
			}
		}
	}
}

undef(%COG_hit);

my $hitnum = 99999;
my %first_BLA_hit;
my %second_BLA_hit;
my %BLA_hit_pair;
my $num_hits = `wc -l $ARGV[1] | awk '{print \$1}'`;
chomp $num_hits;
my %self_blast;
my $threshold = 0.05/$num_hits; # crude Bonferroni correction - conservative for pairing

print STDERR "done\nNow identifying BLAST hits involving unclustered proteins...";

open BLA, $ARGV[1] or die;

# search BLAST file for hits involving unclustered CDSs - find the top two BLAST scores for each unclustered CDS

while (my $line = <BLA>) {
	my @data = split(/\s+/,$line);
	if ($data[0] eq $data[1] && defined($unclustered{$data[0]}) && $unclustered{$data[0]} == 1) {
		unless (defined($self_blast{$data[0]}) && $self_blast{$data[0]} < $data[10]) { 
			$self_blast{$data[0]} = $data[10];
		}
	}
	if (defined($unclustered{$data[0]}) && $unclustered{$data[0]} == 1 && $data[1] ne $data[0] && $data[10] < $threshold && defined($unclustered{$data[1]}) && $unclustered{$data[1]} == 1) {
		unless (defined($first_BLA_hit{$data[0]}) && $data[10] > $first_BLA_hit{$data[0]}) {
			$second_BLA_hit{$data[0]} = $first_BLA_hit{$data[0]};
			$first_BLA_hit{$data[0]} = $data[10];
			$BLA_hit_pair{$data[0]} = $data[1];
		}
	}
}

close BLA;

print STDERR "done\nNow pairing proteins...\n";

my $best_hit;
my @best_hit;
my $recip_best_hit;
my @recip_best_hit;
my %now_clustered;

# identify cases where there is a reciprocal best match between unclustered
# CDSs that exceeds the significance threshold

foreach my $CDS (sort keys %unclustered) {
	if (!(defined($self_blast{$CDS})) || $self_blast{$CDS} >= $threshold) {
		unless (defined($now_clustered{$CDS}) && $now_clustered{$CDS} == 1) {
			if ($CDS =~ /$refstring/) {
				print REF "$CDS,UNCLS\n";
			} else {
				print OUT "$CDS,UNCLS\n";
			}
			$now_clustered{$CDS} = 1;
		}
		if (defined($BLA_hit_pair{$CDS})) {
			unless (defined($now_clustered{$BLA_hit_pair{$CDS}}) && $now_clustered{$BLA_hit_pair{$CDS}} == 1) {
				if ($BLA_hit_pair{$CDS} =~ /$refstring/) {
					print REF "$BLA_hit_pair{$CDS},UNCLS\n";
				} else {
					print OUT "$BLA_hit_pair{$CDS},UNCLS\n";
				}
				$now_clustered{$BLA_hit_pair{$CDS}} = 1;
			}
		}
	} else {
		unless (defined($now_clustered{$CDS}) && $now_clustered{$CDS} == 1) {
			if (!(defined($second_BLA_hit{$CDS})) || $first_BLA_hit{$CDS} == $second_BLA_hit{$CDS}) {
				if (defined($first_BLA_hit{$CDS}) && defined($BLA_hit_pair{$BLA_hit_pair{$CDS}}) && $BLA_hit_pair{$BLA_hit_pair{$CDS}} eq $CDS && (!(defined($second_BLA_hit{$CDS})) || $first_BLA_hit{$CDS} != $second_BLA_hit{$CDS}) && !(defined($now_clustered{$BLA_hit_pair{$CDS}}) && $now_clustered{$BLA_hit_pair{$CDS}} == 1)) {
					print STDERR "Paired proteins $CDS and $BLA_hit_pair{$CDS},E values $first_BLA_hit{$CDS}\n";
					if ($CDS =~ /$refstring/) {
						print REF "$CDS,CLS"."$hitnum\n";
					} else {
						print OUT "$CDS,CLS"."$hitnum\n";
					}
					if ($BLA_hit_pair{$CDS} =~ /$refstring/) {
						print REF "$BLA_hit_pair{$CDS},CLS"."$hitnum\n";
					} else {
						print OUT "$BLA_hit_pair{$CDS},CLS"."$hitnum\n";
					}
					$now_clustered{$CDS} = $now_clustered{$BLA_hit_pair{$CDS}} = 1;
					$hitnum--;
				} else {
					if ($CDS =~ /$refstring/) {
						print REF "$CDS,CLS"."$hitnum\n";
					} else {
						print OUT "$CDS,CLS"."$hitnum\n";
					}
					$now_clustered{$CDS} = 1;
					$hitnum--;
				}
			} else {
				if (defined($unclustered{$BLA_hit_pair{$CDS}}) && $BLA_hit_pair{$BLA_hit_pair{$CDS}} eq $CDS && $first_BLA_hit{$BLA_hit_pair{$CDS}} != $second_BLA_hit{$BLA_hit_pair{$CDS}}) {
					if (defined($now_clustered{$BLA_hit_pair{$CDS}}) && $now_clustered{$BLA_hit_pair{$CDS}} == 1) {
						if ($CDS =~ /$refstring/) {
							print REF "$CDS,CLS"."$hitnum\n";
						} else {
							print OUT "$CDS,CLS"."$hitnum\n";
						}
						$now_clustered{$CDS} = 1;
						$hitnum--;
					} else {
						if ($CDS =~ /$refstring/) {
							print REF "$CDS,CLS"."$hitnum\n";
						} else {
							print OUT "$CDS,CLS"."$hitnum\n";
						}
						if ($BLA_hit_pair{$CDS} =~ /$refstring/) {
							print REF "$BLA_hit_pair{$CDS},CLS"."$hitnum\n";
						} else {
							print OUT "$BLA_hit_pair{$CDS},CLS"."$hitnum\n";
						}
						$now_clustered{$CDS} = $now_clustered{$BLA_hit_pair{$CDS}} = 1;
						print STDERR "Paired proteins $CDS and $BLA_hit_pair{$CDS}, E value $first_BLA_hit{$CDS}\n";
						$hitnum--;
					}
				} else {
					if ($CDS =~ /$refstring/) {
						print REF "$CDS,CLS"."$hitnum\n";
					} else {
						print OUT "$CDS,CLS"."$hitnum\n";
					}
					$hitnum--;
				}
			}
		}
	}
}

close OUT;
close REF;

print STDERR "done\nCompleted post_process_COGs.pl\n";
