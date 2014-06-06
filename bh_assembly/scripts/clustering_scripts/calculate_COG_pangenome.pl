#! /usr/bin/perl -w
use strict;
use warnings;

###################################################################################
# Submitted by extract_annotate_orthologues.pl | Submits nothing                  #
# Calculates the core and pangenome curves and fits based on the presence/absence #
# matrix                                                                          #
###################################################################################


# set up script_dir variable

#JH my $script_dir = $ENV{SCRIPT_DIR};
my $script_dir = "~/consults/bh_assembly/scripts/clustering_scripts/";

my $current_dir = `pwd | tr -d "\n"`;
my $out_dir = "$current_dir"."/pangenome/";
system "mkdir $out_dir";

# copy over R script

system "cp $script_dir/odds_ratio.R .";

my $exclude;

if (defined($ARGV[2]) && length($ARGV[2]) > 0) {
	$exclude = $ARGV[2];
}

open INFO, $ARGV[0] or die print STDERR "Cannot open information file $ARGV[0]\n";
print STDERR "Reading strain information file\n";

my %translate;

foreach (<INFO>) {
	chomp;
	my @data = split(/\s+/,$_);
	unless (defined($exclude) && $data[0] == $exclude) {
		my $seq = "seq".$data[0];
		$data[1] =~ s/.fa|.dna//g;
		$translate{$seq} = $data[1];
	}
}

close INFO;

my %cluster_count;
my %group_count;
my %strain_clusters;
my $num = 0;

open CLS, $ARGV[1] or die print STDERR "Cannot open clusters file $ARGV[1]\n";
print STDERR "Reading cluster information file\n";

open UCL, "> $out_dir"."/unclustered_CDS.out";

while (my $line = <CLS>) {
	chomp $line;
	my @data = split(/\,/,$line);
	my @strain = split(/orf/,$data[0]);
	if (defined($translate{$strain[0]})) {
		if (defined($data[1])) {
			push(@{$strain_clusters{$translate{$strain[0]}}},$data[1]);
		} else {
			print UCL "$strain[0]\t$data[0]\n";
			push(@{$strain_clusters{$translate{$strain[0]}}},"UNIQ.$num");
			$num++;
		}
	}
}

close CLS;

close UCL;

# calculate pangenome

my @unused;
my %seen;
my %always_seen;
my %in_strain;
my %core;
my %pan;
my %specific;
my %shared;
my $count = 0;

print STDERR "Starting pangenome calculation\n";

my $max = scalar(keys %strain_clusters);
while ($count <= 100) {
	@unused = keys %strain_clusters;
	my $seqcount = 1;
	while ($seqcount <= $max) {
		my $rand = int(rand(scalar(@unused)));
		$shared{$seqcount}{$count} = 0;
		$specific{$seqcount}{$count} = 0;
		foreach my $cluster (@{$strain_clusters{$unused[$rand]}}) {
			if ($seqcount == 1) {
				$always_seen{$cluster} = 1;
			}
			if (defined($always_seen{$cluster}) && $always_seen{$cluster} == 1) {
				$shared{$seqcount}{$count}++;
			}
			unless (defined($seen{$cluster}) && $seen{$cluster} == 1) {
				$specific{$seqcount}{$count}++;
				$seen{$cluster} = 1;
			}
			$in_strain{$cluster} = 1;
		}
		foreach my $cluster (keys %always_seen) {
			unless (defined($in_strain{$cluster}) && $in_strain{$cluster} == 1) {
				delete($always_seen{$cluster});
			}
		}
		$core{$seqcount}{$count} = scalar(keys %always_seen);
		if ($seqcount > 1) {
			$pan{$seqcount}{$count} = $pan{($seqcount-1)}{$count}+$specific{$seqcount}{$count};
		} else {
			$pan{$seqcount}{$count} = $specific{$seqcount}{$count};
		}
		splice(@unused,$rand,1);
		undef(%in_strain);
		$seqcount++;
	}
	undef(%seen);
	undef(%always_seen);
	$count++;
}

print STDERR "Calculation complete; printing data\n";

# print pangenome
 
open CORE, "> $out_dir"."/core_genome.txt";
open PAN, "> $out_dir"."/pan_genome.txt";
open SPEC, "> $out_dir"."/strain_specific.txt";
open SHARE, "> $out_dir"."/shared_genome.txt";

foreach my $seqcount (sort {$a<=>$b} keys %pan) {
	print PAN "$seqcount";
	print CORE "$seqcount";
	print SPEC "$seqcount";
	print SHARE "$seqcount";
	foreach my $count (sort {$a<=>$b} keys %{$pan{$seqcount}}) {
		print PAN ",$pan{$seqcount}{$count}";
		print CORE ",$core{$seqcount}{$count}";
		print SPEC ",$specific{$seqcount}{$count}";
		print SHARE ",$shared{$seqcount}{$count}";
	}
	print PAN "\n";
	print CORE "\n";
	print SPEC "\n";
	print SHARE "\n";
}

close PAN;
close CORE;
close SPEC;

chdir "$out_dir";
system "R CMD BATCH $script_dir/pangenome_analysis.R";
