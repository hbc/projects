#! /usr/bin/perl -w
use strict;
use Bio::SeqIO;

###################################################################################
# Submitted by check_COG_analysis.pl | Submits RAxML and BAPS jobs                #
# Creates the whole genome alignment in FASTA and BAPS formats                    #
###################################################################################

# set up script_dir variable

my $script_dir = $ENV{SCRIPT_DIR};

# first concatenate COG reports and clear up files

if (-e "all.strains.COG.report") {
	system "rm all.strains.COG.report";
}

system "echo 'Name	#_entries       Mean_length     Aln_length      DNA_ps  AA_ps   SigPep  #TMhelix        PFAM' > all.strains.COG.report";
system "cat aligndir*/cluster.*.report >> all.strains.COG.report";
system "rm aligndir*/cluster.*.report";

# then concatenate all tab file components into a large tab file

if (-e "core_genome.tab") {
	system "rm core_genome.tab";
}

system "cat aligndir*/partial.*.core_genome.tab > core_genome.tab";
system "rm aligndir*/partial.*.core_genome.tab";

# set up other variables

my $count = $ARGV[0];
my $reference = $ARGV[1];
$reference =~ s/embl/dna/g;
my $refnum = $ARGV[2];
my $n = 1;
my %align;
my $strain;
my %tracking;
my %rev_tracking;
my $max_index = 0;
my $dir = `pwd`;
chomp $dir;
my $range = 1000;
my $jobid = int(rand($range));

open DATA, 'strain.info';

foreach (<DATA>) {
	chomp;
	my @data = split(/\t/,$_);
	$tracking{$data[0]} = $data[1];
	$rev_tracking{$data[1]} = $data[0];
	if ($data[0] > $max_index) {
		$max_index = $data[0];
	}
}

close DATA;

# process core_genome.tab to get the list of core genome genes

my @core_genes;

open TAB, 'core_genome.tab' or die print STDERR "Could not open core_genome.tab - may not be any core clusters!\n";

foreach (<TAB>) {
	chomp;
	if (/Cluster/) {
		my @data = split(/\s+/,$_);
		push(@core_genes,$data[2]);
	}
}

# process alignments

foreach $n (@core_genes) {
	my $file_path = `ls aligndir_*/cluster.$n.DNA.aln | tr -d "\n"`;
	my $head = `head -n 1 $file_path`;
	chomp $head;
	my $seqs = `grep -c ">" $file_path`;
	chomp $seqs;
	if (defined($file_path) && length($file_path) > 0 && -e $file_path) {
		open IN, $file_path or die print STDERR "Cannot open $file_path!\n";
		if ($head =~ /^>/ && $seqs == ($count-1)) {
			print STDERR "Alignment $n read correctly!\n";
			foreach (<IN>) {
				chomp;
				if (/^>/) {
					my @name = split(/orf/,$_);
					$strain = $name[0];
					$strain =~ s/\>seq//g;
				} else {
					if (defined($align{$strain})) {
						$align{$strain} .= $_;
					} else {
						$align{$strain} = $_;
					}
				}
			}
			close IN;
			system "rm $file_path";
		} else {
			print STDERR "Alignment $n is not correctly formatted\n";
			close IN;
		}
	} else {
		print STDERR "Could not open alignment $n\n";
	}
	$n++;
}
$n = 1;

open OUT, "> core_genes.aln";

# print out concatenated protein alignment for all strains

foreach my $index (sort keys %tracking) {
	unless ($index == $refnum) {
		my $strain = $tracking{$index};
		print OUT ">$strain\n$align{$index}\n";
	}
}

close OUT;

# run BAPS analysis

my $pseudo_ref = `head -n 1 core_genes.aln | sed 's/>//g' | tr -d "\n"`;

system "~croucher/Scripts/summarise_snps_from_alignment.pl  -a core_genes.aln -f fasta -o core_genes.out";

open COR, "core_genes.out.aln" or die print STDERR "Could not open core_genes.out.aln!\n";

open BAP, "> core_genes.baps";

my $header = 0;

foreach (<COR>) {
	if ($header == 0) {
		$header = 1;
	} else {
		my @data = split(/\s+/,$_);
		$data[1] = uc($data[1]);
		$data[1] =~ s/A/1\t/g;
		$data[1] =~ s/C/2\t/g;
		$data[1] =~ s/G/3\t/g;
		$data[1] =~ s/T/4\t/g;
		$data[1] =~ s/-/-9\t/g;
		$data[1] =~ s/N/-9\t/g;
		$data[1].=$rev_tracking{$data[0]};
		print BAP "$data[1]\n";
	}
}

close BAP;

close COR;

open RUN, "> core_baps.runfile";

my $max_baps = int($max_index/5); # on average, each BAPS group should contain at least five strains

print RUN "datafile('$dir/core_genes.baps')\n";
print RUN "mixturetype('mix')\n";
print RUN "initialk($max_baps)\n";	
print RUN "fixedk('no')\n";
print RUN "datatype('numeric')\n";
print RUN "outputmat('$dir/baps_clustering.mat')\n";

close RUN;

system "bsub -J $jobid"."BAPS -o baps.o -e baps.e -M 10000000 -R 'select[mem>10000] rusage[mem=10000]' -q long_serial 'module load bio/BAPS;run_baps6.sh /n/sw/matlab-2010a/MATLAB_Compiler_Runtime/v713/ core_baps.runfile'";
system "bsub -w \"ended($jobid"."BAPS)\" -o baps.o -e baps.e $script_dir/process_BAPS.pl strain.info baps_clustering.mat.txt";

# run phylogenetic analysis
my $raxrand = int(rand(10000));

system "bsub -M 2000000 -R 'select[mem>4000] rusage[mem=4000]' -o rax.o -e rax.e -q long_serial /n/sw/odyssey-apps/RAxML-7.0.4/bin/raxmlHPC -s core_genes.out.aln -m GTRGAMMA -n all.strains -p $raxrand";

# tidy up
system "rm *.aa seq.*.csv";
system "rm -rf aligndir_* tmp all-edges.txt cog-edges.txt cogtriangle_run.sh concatenation_script.sh LSE.csv all.strains.icm cogtriangle_run*.sh";# all.strains.csv
system "for f in log.* Slog.*; do echo \$f >> all.log; cat \$f >> all.log; rm \$f; done";
system "for f in err.* Serr.*; do echo \$f >> all.err; cat \$f >> all.err; rm \$f; done";
system "rm JobArray.* FilterBlast* BlastJob* run_baps5.sh";
system "rm cognitor.log addendum.out conflict.txt odds_ratio.R";

if (-e "previous.blast_outputs.log") {
	system "rm previous.blast_outputs.log";
}

if (-e "blast_outputs.log") {
	system "rm blast_outputs.log";
}

print STDERR "Completed orthologue_alignment_processing.pl\n";
