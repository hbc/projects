#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by run_clustering.pl | Submits BLAT jobs, COGtriangles runs and       #
# extract_annotate_orthologues.pl                                                 #
# Checks whether gene predictions have run on each genome appropriately           #
###################################################################################

# set up script_dir variable

my $script_dir = $ENV{SCRIPT_DIR};

# process input

my $count = $ARGV[0];
my $jobid = $ARGV[1];
my $assembly = $ARGV[2];
my $refnum = $ARGV[3];
my $reference = $ARGV[4];

my $dir = `pwd`;
chomp $dir;

# read in strain names

my %tracking;

open IN, "strain.info";

foreach (<IN>) {
	chomp;
	my @data = split(/\s+/,$_);
	#$data[1] =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.CDS.mfa/g;
	$tracking{$data[0]} = $data[1];
}

close IN;

# check whether a sufficient number of proteins have been output from each genome

my @file_nums = (1..$count);
my @rerun;

# retrieve the number of CDSs predicted in the reference sequence

my $ref_mfa = $tracking{$refnum};
$ref_mfa =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.CDS.mfa/g;
my $ref_CDS = `grep -c '>' $ref_mfa | tr -d "\n"`;

# check that a reasonable number of CDSs have been predicted for each assembly

foreach my $num (@file_nums) {
	my $CDS_mfa = $tracking{$num};
	$CDS_mfa =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.CDS.mfa/g;
	if (-e "$dir/$CDS_mfa" && -e "$dir/seq.$num.csv" && -e "$dir/$tracking{$num}.report") {
		my $test_CDS = `grep -c '>' $CDS_mfa`;
		chomp $test_CDS;
		if ($test_CDS < (0.75*$ref_CDS)) {		# check for good assembly and gene prediction
			push(@rerun,$num);
		}
	} else {
		push(@rerun,$num);
	}
}

$" = ",";

# rerun jobs where necessary; else run clustering jobs

if (scalar(@rerun) > 0) {
	if ($assembly == 1) {
		system "bsub -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid"."glim[@rerun]".' -q normal_serial -o log.%I -e err.%I ./JobArray.\${LSB_JOBINDEX}';
	} else {
		system "bsub -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' -J $jobid"."glim[@rerun]".' -q normal_serial -o log.%I -e err.%I ./JobArray.\${LSB_JOBINDEX}';
	}
	system "bsub -J \"$jobid\".CHK -w \"ended($jobid"."glim)\" -o all.log -e all.err $script_dir/clustering_checkpoint.pl $count $jobid $assembly $refnum $reference";
} else {
	# concatenate the assembly reports
	
	if (-e "$dir/all.strains.assembly_stats.out") {
		system "rm $dir/all.strains.assembly_stats.out";
	}
	
	system "echo -e \"File\tLength\t#_contigs\tN50\t#_CDS\" > all.strains.assembly_stats.out";
	
	foreach my $num (keys %tracking) {
		unless ($num == $refnum) {
			system "cat $tracking{$num}.report >> all.strains.assembly_stats.out";
			system "rm $tracking{$num}.report";
		}
	}
	
	# submit the concatenation script
	
	system "bsub -J \"$jobid"."CAT\" -o all.log -e all.err -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' ./concatenation_script.sh";
	
	# run BLAST comparisons
	
	system "bsub -w \"ended($jobid"."CAT)\" -M 2000000 -R 'select[mem>2000] rusage[mem=2000] hname!=bc-11-4-10' -J $jobid"."UBLA[1-$count]".' -q normal_serial -o log.%I -e err.%I ./BlastJob.\${LSB_JOBINDEX}';

	system "bsub -w \"ended($jobid"."CAT)\" -M 2000000 -R 'select[mem>2000] rusage[mem=2000] hname!=bc-11-4-10' -J $jobid"."FBLA[1-$count]".' -q normal_serial -o log.%I -e err.%I ./FilterBlast.\${LSB_JOBINDEX}';

	
	# run Cogtriangles on the BLAT comparisons

	# calculate approximate memory requirements - 13 Gb per 100 strains, 20 Gb starting level - this gives the 98 Gb needed to process the SPARC collection, need to fix this for other bugs
	
	my @blat_array = (1..$count);
	my @filter_array = (1..$count);
	$" = ",";
		
	system "bsub -w \"ended($jobid"."UBLA[1-$count]) && ended($jobid"."FBLA[1-$count])\" -o all.log -e all.err $script_dir/blast_checkpoint.pl @blat_array @filter_array $count $jobid $reference $refnum";
	
	#my $bignum = 10000000;
#	my $smallnum = 10000;
#	my $queue = "long";
#
#	if ($count > 100) {
#		$bignum = (int($count/100)*13000000)+20000000;
#		$smallnum = $bignum/1000;
#		if ($smallnum > 20000) {
#			$queue = "hugemem";
#		}
#	}
#	
#	system "bsub -w \"ended($jobid"."UBLA[1-$count]) && ended($jobid"."FBLA[1-$count])\" -o all.log -e all.err -q long -J $jobid.cogA -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' ./cogtriangle_run_A.sh";
#	
#	system "bsub -w \"ended($jobid.cogA)\" -o all.log -e all.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' -J $jobid.cogB ./cogtriangle_run_B.sh";
#	
#	system "bsub -w \"ended($jobid.cogB)\" -o all.log -e all.err -q long -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid.cogC ./cogtriangle_run_C.sh";
#	
#	system "bsub -w \"ended($jobid.cogC)\" -o all.log -e all.err -q long -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid.cluchk $script_dir/post_process_COGs.pl all.strains.cls.out.csv ./blaf/filtered.all.strains.blast.tab all.strains.csv $refnum";
#	
#	my $embl = $reference;
#	$embl =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.embl/g;
#	
#	system "bsub -o all.log -e all.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' -w \"ended($jobid.cluchk)\" $script_dir/extract_annotate_orthologues.pl $count $refnum $embl";
#	
	# tidy up files that are no longer needed
	
	foreach my $num (@file_nums) {
		my $file_dna = $tracking{$num};
		unless ($file_dna =~ /.dna/) {
			$file_dna =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.dna/g;
			system "rm $file_dna";
		}
		my $file_fa = $tracking{$num};
		unless ($file_fa =~ /.fa/) {
			$file_fa =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.fa/g;
			system "rm $file_fa";
		}
	}
	
	my $refstem = $tracking{$refnum};
	$refstem =~ s/_1.fastq|.dna|.fasta|.seq|.fa//g;
	system "rm $tracking{$refnum}.report $refstem.train $refstem.longorfs reference_training.sh ";
}

print STDERR "Completed clustering_checkpoint.pl\n";
