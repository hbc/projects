#! /usr/bin/perl
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::AlignIO::largemultifasta;
use strict;
use warnings;

# parse input options

my $aln_file;
my $aln_format = "Fasta";
my $ref_format = "EMBL";
my $outfile = "summarised_snps.out";
my $ref_file;
my $help = 0;

GetOptions (
        "a=s" => \$aln_file,
	"f=s" => \$aln_format,
	"n=s" => \$ref_format,
	"o=s" => \$outfile,
	"r=s" => \$ref_file,
        "h" => \$help
);

if ($help == 1) {
	croak();
}

unless (-e $aln_file) {
	print STDERR "Cannot open alignment file $aln_file!\n" and croak();
}

# parse alignment

my %aln;
my $length;
my %names;

my $alnObject = Bio::SeqIO->new(-file => $aln_file, -format => $aln_format);

my $alength = 0;
while (my $seq = $alnObject->next_seq) {
	$names{$seq->id} = 1;
	if ($seq->length > $alength) {
		$alength = $seq->length;
	}
}

my $offset = 1;
my $window = 99999;
my %polymorphic;

while (($offset+$window) <= $alength) {
	$alnObject = Bio::SeqIO->new(-file => $aln_file, -format => $aln_format);
	while (my $seq = $alnObject->next_seq) {
		my $name = $seq->id;
        	my $aln_seq = $seq->subseq($offset,($offset+$window));
		@{$aln{$name}} = split(//,$aln_seq);
	}
	extract_poly($offset);
	undef(%aln);
	$offset = $offset+$window+1;
}


if (($alength-$offset) > 1) {
	while (my $seq = $alnObject->next_seq) {
		my $name = $seq->id;
        	my $aln_seq = $seq->subseq($offset,$length);
		@{$aln{$name}} = split(//,$aln_seq);
	}
	extract_poly($offset);
	undef(%aln);
	$offset = $offset+$window+1;
}

# check if changes cause alterations in the protein coding sequences

my $cds_count = 0;
my %cds_start;
my %cds_end;
my %cds_name;

if (defined($ref_file) && length($ref_file) > 1) {
	my $annotationObject = Bio::SeqIO->new(-file => $ref_file, -format => $ref_format);
	my $annotation = $annotationObject->next_seq;
	foreach my $cds ($annotation->get_SeqFeatures) {
		if ($cds->primary_tag eq 'CDS' && !($cds->has_tag('pseudo'))) {
			$cds_name{$cds_count} = get_name($cds);
			$cds_start{$cds_count} = $cds->location->start;
			$cds_end{$cds_count} = $cds->location->end;
			$cds_count++;
		}
	}
	
}

open OUT, "> $outfile" or print STDERR "Cannot open output file $outfile!\n" and croak();
open ALN, "> $outfile.aln" or print STDERR "Cannot open output file $outfile.aln!\n" and croak();

if (scalar(keys %polymorphic) == 0) {
	print OUT "No polymorphic sites found!\n";
	close OUT;
	exit(0);
}

my $refsite = (keys %polymorphic)[0];

# print alignment

my $taxon_count = scalar(keys %{$polymorphic{$refsite}});
my $alignment_length = scalar(keys %polymorphic);

print ALN "   $taxon_count $alignment_length\n";

# print summary out file

print OUT "Position";

foreach my $name (sort keys %{$polymorphic{$refsite}}) {
	print ALN "$name\t";
	print OUT "\t$name";
	foreach my $site (sort {$a<=>$b} keys %polymorphic) {
		print ALN "$polymorphic{$site}{$name}"
	}
	print ALN "\n";
}

if (defined($ref_file) && length($ref_file) > 1) {
	print OUT "\tAffectedGene";
}

print OUT "\n";

# print one line per polymorphic site

foreach my $site (sort {$a<=>$b} keys %polymorphic) {
	my $pos = $site;
	print OUT "$pos";
	foreach my $name (sort keys %{$polymorphic{$site}}) {
		print OUT "\t$polymorphic{$site}{$name}"
	}
	if (defined($ref_file) && length($ref_file) > 1) {
		my $affectedGene = '-';
		foreach $cds_count (keys %cds_start) {
			if ($pos >= $cds_start{$cds_count} && $pos <= $cds_end{$cds_count}) {
				$affectedGene = $cds_name{$cds_count};
			}
		}
		print OUT "\t$affectedGene";
	}
	print OUT "\n";
}

close OUT;

# subroutines

sub croak {
	print STDERR "\nsummarise_snps_from_alignment.pl -o [output file name (default: summarised_snps.out)] -a [alignment file] -f [alignment format (default: Fasta)] -r [reference annotation] -n [reference annotation format (default: EMBL)]\n\n";
	exit(0);
}

sub get_name {
	my $cds = shift;
	my @quals = ('locus_tag','gene','systematic_id','primary_name');
	my $name;
	foreach my $q (@quals) {
		if ($cds->has_tag($q) && !(defined($name))) {
			foreach my $v ($cds->get_tag_values($q)) {
				$name = $v;
			}
		}
	}
	return($name);
}

sub extract_poly {
	my $i = shift;
	my $n = 0;
	my %bases;
	foreach my $name (keys %aln) {
		for ($n = 0; $n <= $#{$aln{$name}}; $n++) {
			my $base = uc($aln{$name}[$n]);
			unless ($base eq 'N' || $base eq 'n' || $base eq '-') {
				$bases{$n}{$base}++;
			}
		}
	}
	foreach my $n (keys %bases) {
		if (scalar(keys %{$bases{$n}}) > 1) {
			foreach my $name (keys %aln) {
				my $base = uc($aln{$name}[$n]);
				$polymorphic{$n+$i}{$name} = $base;
			}
		}
	}
}
