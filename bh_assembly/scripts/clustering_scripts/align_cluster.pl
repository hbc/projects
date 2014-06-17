#! /usr/bin/perl -w
use strict;
use Bio::SeqIO;
use DBI;
use Bio::AlignIO;
#use lib "~croucher/bin/Pfam/PfamScan.pm";
#use lib "$ENV{'HOME'}/bin/Pfam";
use Bio::Pfam::Scan::Seq;
use Bio::Pfam::Scan::PfamScan;

###################################################################################
# Submitted by extract_annotate_orthologues.pl | Submits nothing                  #
# Extracts the protein and DNA sequences, runs an alignment, extracts polymorphic #
# sites. Also runs analysis tools signalP, TMHMM, PFAM, makes a summary report    #
###################################################################################

# genetic code

my %genetic_code = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # Leucine
'TTG' => 'L', # Leucine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '_', # Stop
'TAG' => '_', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '_', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # Leucine
'CTC' => 'L', # Leucine
'CTG' => 'L', # Leucine
'CTT' => 'L', # Leucine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # Glutamine
'CAG' => 'Q', # Glutamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # Isoleucine
'ATC' => 'I', # Isoleucine
'ATT' => 'I', # Isoleucine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # Glutamic Acid
'GAG' => 'E', # Glutamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G',  # Glycine
# two base inference with ambiguous bases
'GGN' => 'G',	# Glycine
'CTN' => 'L',	# Leucine
'CGN' => 'R',	# Arginine
'TCN' => 'S',	# Serine
'CCN' => 'P',	# Proline
'GTN' => 'V',	# Valine
'GCN' => 'A',	# Alanine
'ACN' => 'T',	# Threonine
# two base inferences
'GG' => 'G',	# Glycine
'CT' => 'L',	# Leucine
'CG' => 'R',	# Arginine
'TC' => 'S',	# Serine
'CC' => 'P',	# Proline
'GT' => 'V',	# Valine
'GC' => 'A',	# Alanine
'AC' => 'T'	# Threonine
);

# set up variables

my @editable = @ARGV;
my $reference = shift(@editable);
my $refstring = shift(@editable);
my $core = shift(@editable);
my $SCcount = shift(@editable);
my $cluster_name = shift(@editable);
my $working_dir = shift(@editable);
my @orthos = @editable;
my $non_ref = 0;
my @fetch;

# change directory

chdir "$working_dir";

# connect to SQL database

my $dbfile = "../protein_sequences.db";
my $db = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
my $DNA_dbfile = "../DNA_sequences.db";
my $DNA_db = DBI->connect("dbi:SQLite:dbname=$DNA_dbfile","","");

# retrieve sequences

my $names;
my $sequences;

open OUT, "> cluster.$SCcount.seq";
open DNA, "> cluster.$SCcount.DNA.seq";

foreach my $ortho (sort @orthos) {
	unless ($ortho =~ /$refstring/) {
		$non_ref++;
		push(@fetch,$ortho);
	}
}

$" = "', '";

if (scalar(@fetch) > 0) {	
	my $command = $db->prepare("SELECT * FROM Proteins WHERE Label in ('@fetch')");
	$command->execute();
	$command->bind_columns(\$names, \$sequences);
	while ($command->fetch()) {
		print OUT ">$names\n$sequences\n";
	}
	$command = $DNA_db->prepare("SELECT * FROM DNA WHERE Label in ('@fetch')");
	$command->execute();
	$command->bind_columns(\$names, \$sequences);
	while ($command->fetch()) {
		print DNA ">$names\n$sequences\n";
	}
}

close OUT;
close DNA;

$" = " ";

# hashes for calculating polymorphic sites and sequence lengths

my %prot_aln_poly;
my %dna_aln_poly;
my $prot_site = 1;
my $dna_site = 1;
my $total_dna_length = 0;
my %protein_lengths;
my %prot_aln;
my $prot_aln;
my $aln;

if ($non_ref > 1) {									# only align sequences where there are multiple representatives
	# align protein sequences

	system "muscle -in cluster."."$SCcount".".seq -out cluster."."$SCcount".".aln";

	# parse protein alignment

	$prot_aln = Bio::AlignIO->new(-file => "cluster.$SCcount.aln", -format => 'Fasta');

} else {
	
	# use single entry as "alignment"
	
	$prot_aln = Bio::AlignIO->new(-file => "cluster.$SCcount.seq", -format => 'Fasta');
	
}

$aln = $prot_aln->next_aln();

foreach my $seq ($aln->each_seq()) {
	my $name = $seq->id;
	my $aln_seq = $seq->seq;
	$prot_aln{$name} = $aln_seq;
	my @amino_acids = split(//,$aln_seq);
	$protein_lengths{$name} = scalar(@amino_acids);
	$prot_site = 1;
	foreach my $aa (@amino_acids) {
		unless ($aa eq '-' || $aa eq 'X' || $aa eq '?') {
			$prot_aln_poly{$prot_site}{$aa} = 1;
		}
		$prot_site++;
	}
}

# find sequence of median length and run Pfam scan

my $num_entries = scalar(keys %protein_lengths);

my $counter = 1;
my $pfam_domains;
my $signalP;
my $tmhmm;

foreach my $name (sort {$protein_lengths{$a}<=>$protein_lengths{$b}} keys %protein_lengths) {
	if ($counter == int(0.5*$num_entries) || $num_entries == 1) {
		open PFS, "> cluster.$SCcount.pfam.seq";
		my $seq_samp = $prot_aln{$name};
		$seq_samp =~ s/-//g;
		print PFS ">$name\n$seq_samp\n";
		close PFS;
		my $hostname = `hostname | tr -d "\n"`;
		until (defined($pfam_domains) && length($pfam_domains) > 0) {
			
			if (-e "output.$SCcount.pfam") {
				system "rm output.$SCcount.pfam";
			}
			
			my $dir;
			my @hmmlib = ("Pfam-A.hmm");
			my ($e_seq, $e_dom, $b_seq,$b_dom, $outfile) = (0.05, 0.05, 0, 0, "output.$SCcount.pfam");
			
			$dir = ("/n/home10/croucher/bin/Pfam/");
			
			if (-e "$dir/$hmmlib[0]") {
				print STDERR "Found file $hmmlib[0] in $dir!\n";
			}
			
			# build PFAM object
			my $ps = Bio::Pfam::Scan::PfamScan->new(
			  -fasta        => "cluster.$SCcount.pfam.seq",
			  -hmmlib       => \@hmmlib,
			  -e_seq        => $e_seq,
			  -e_dom        => $e_dom,
			  -b_seq        => $b_seq,
			  -b_dom        => $b_dom,
			  -dir		=> $dir
			);
			
			# run PFAM search
			
			$ps->search;

			$ps->write_results( $outfile, $e_seq, $e_dom, $b_seq, $b_dom );
			
			if (-e "output.$SCcount.pfam") {
				$pfam_domains = `grep PF output.$SCcount.pfam | grep -v '#' | awk '{print \$7}' | sort | uniq | tr "\n" ","`;
				
				unless (defined($pfam_domains) && length($pfam_domains) > 0) {
					$pfam_domains = "-";
				}
			}
			
		}
		
		$signalP = `signalp -t gram+ -f summary cluster.$SCcount.pfam.seq | grep Prediction | awk -F ': ' '{print \$2}' | sed "s/ /_/g"`;
		chomp $signalP;
		$tmhmm = `tmhmm cluster.$SCcount.pfam.seq 2> /dev/null | grep -c 'TMhelix'`;
		chomp $tmhmm;
		system "rm cluster.$SCcount.pfam.seq";
		last;
	}
	$counter++;
}

# parse DNA sequences

my %DNA_seq;

my $dna_seq = Bio::SeqIO->new(-file => "cluster.$SCcount.DNA.seq", -format => 'Fasta');

while (my $seq = $dna_seq->next_seq) {
	my $name = $seq->id;
	my $aln_seq = $seq->seq;
	$DNA_seq{$name} = $aln_seq;
	$total_dna_length+=length($aln_seq);
}

# produce back translated alignment

my %DNA_aln;
my $max_length = 0;

if ($non_ref > 1) {
	foreach my $name (keys %DNA_seq) {
		my @prot_aln = split(//,$prot_aln{$name});
		my @dna_seq = unpack("(A3)*",$DNA_seq{$name});
		
		# need to remove the final stop codon when present
		if (defined($genetic_code{$dna_seq[$#dna_seq]}) && $genetic_code{$dna_seq[$#dna_seq]} eq '_') {
			pop(@dna_seq);
			$total_dna_length = $total_dna_length-3;
		}
		
		$DNA_aln{$name} = "";
		my $n = 0;
		my $aa;
		my $codon;

		while ($#prot_aln > -1 && $#dna_seq > -1) {
			$aa = shift(@prot_aln);
			# added
			unless (defined($aa) && $aa eq '_') {
				if ($aa eq '-' || $aa eq 'X') {
					$DNA_aln{$name}.="---";
				} else {
					$codon = shift(@dna_seq);
					unless (defined($genetic_code{$codon}) && $genetic_code{$codon} eq '_') {
						while ($codon =~ /N/ && !(defined($genetic_code{$codon}))) {
							$codon = shift(@dna_seq);
						}
						if ($genetic_code{$codon} eq $aa) {
							$DNA_aln{$name}.=$codon;
						} else {
							my $length_thing = scalar(@dna_seq);
							my @whole = unpack("(A3)*",$DNA_seq{$name});
							die print STDERR "Doesn't work:\nDNA: @dna_seq\nProtein: @prot_aln\nname $name codon $codon, aa $aa, should be $genetic_code{$codon}!\nDNA:\t$DNA_seq{$name}\nAA:$prot_aln{$name}\t\n------\nUsed to be: @whole\n";
						}
					}
				}
			}
		}
		
		if (length($DNA_aln{$name}) > $max_length) {
			$max_length = length($DNA_aln{$name});
		}
		
	}

	foreach my $name (keys %DNA_aln) {
		my $pad = $max_length - length($DNA_aln{$name});
		$DNA_aln{$name}.="-"x$pad;
		my @dna_bases = split(//,$DNA_aln{$name});
		$dna_site = 1;
		foreach my $base (@dna_bases) {
			unless ($base eq '-' || $base eq 'X' || $base eq '?') {
				$dna_aln_poly{$dna_site}{$base} = 1;
			}
			$dna_site++;
		}
	}
} else {
	$max_length = $total_dna_length;
}

# identify polymorphic sites in DNA and protein sequences

open REP, "> cluster.$SCcount.report" or die;

my $mean_length = $total_dna_length/$num_entries;
my $rounded_mean = sprintf("%.2f",$mean_length);

my $dna_polymorphic_sites = 0;
my $protein_polymorphic_sites = 0;

if ($non_ref > 1) {
	foreach $dna_site (keys %dna_aln_poly) {
		if (scalar(keys %{$dna_aln_poly{$dna_site}}) > 1) {
			$dna_polymorphic_sites++;
		}
	}

	foreach $prot_site (keys %prot_aln_poly) {
		if (scalar(keys %{$prot_aln_poly{$prot_site}}) > 1) {
			$protein_polymorphic_sites++;
		}
	}
}

# print out summary report

print REP "$cluster_name\t$num_entries\t$rounded_mean\t$max_length\t$dna_polymorphic_sites\t$protein_polymorphic_sites\t";


if (defined($signalP) && length($signalP) > 1) {
	print REP "$signalP\t";
} else {
	print REP "-\t";
}

if (defined($tmhmm) && length($tmhmm) > 1) {
	print REP "$tmhmm\t";
} else {
	print REP "0\t";
}

if (defined($pfam_domains) && length($pfam_domains) > 1) {
	print REP "$pfam_domains\n";
} else {
	print REP "-\n";
}

# use core genes to generate a phylogeny, mark these up relative to the reference genome

if ($core ne "0") {

	# print out DNA alignment

	open ALN, "> cluster.$SCcount.DNA.aln";

	foreach my $name (keys %DNA_aln) {
		unless ($name =~ /$refstring/) {
			print ALN ">$name\n$DNA_aln{$name}\n";
		}
	}

	close ALN;

	if ($core ne "-") {
		my $infile = Bio::SeqIO->new(-file => $reference, -format => "EMBL");
		my $annotation = $infile->next_seq();

		open TAB, "> partial.$SCcount.core_genome.tab" or die print STDERR "Cannot print features to tab file core_genome.tab\n";

		for my $gene ($annotation->get_SeqFeatures) {
			if ($gene->primary_tag eq "CDS") {
				for my $name ($gene->get_tag_values('label')) {
					if ($name eq $core) {
						my $start = $gene->location->start;
						my $end = $gene->location->end;
						if ($gene->location->strand == 1) {
							print TAB "FT   CDS             $start..$end\n";
						} else {
							print TAB "FT   CDS             complement($start..$end)\n";
						}
						print TAB "FT                   /note=Cluster $SCcount\nFT                   /gene=$cluster_name\n";
					}
				}
			}
		}

		close TAB;
	}
}

# if this process completes, remove as many files as possible

if (-e "cluster.$SCcount.report") {
	if ($core ne "0" && -e "partial.$SCcount.core_genome.tab") {
		system "rm cluster.$SCcount.*seq cluster.$SCcount.aln"; #JH mod to remove period before wildcard, can't find file otherwise
	} elsif ($core eq "0") {
		system "rm cluster.$SCcount.*seq cluster.$SCcount.*aln";#JH mod to remove period before wildcard, can't find file otherwise
	}
}
