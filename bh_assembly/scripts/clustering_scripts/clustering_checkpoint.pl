#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by run_clustering.pl | Submits BLAT jobs, COGtriangles runs and       #
# extract_annotate_orthologues.pl                                                 #
# Checks whether gene predictions have run on each genome appropriately           #
###################################################################################

# set up script_dir and slurm variables
my $script_dir = $ENV{SCRIPT_DIR};
my $slurmqueue = $ENV{SLURMQUEUE};
my $slurmtime = $ENV{SLURMTIME};
my $slurmmem = $ENV{SLURMMEM};
my $slurmexclude = $ENV{SLURMEXCLUDE};

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

my @file_nums = (1..$count);

# check if any samples need to be rerun
my @rerun;

foreach my $num (@file_nums) {
	unless ($num == $refnum) {
		my $done = $tracking{$num};
		$done =~ s/_1.fastq.gz/.fa.done/g;
		unless (-e "$dir/$done") {
			push(@rerun,$num);
		} 
	} 
}

my $jobarrayid;
my $clustercheckpointid;

if (scalar(@rerun) > 0) {
	if ($assembly == 1) {
		write_suffix_array_slurm_script("ReRunJobArrays.sh", "JobArray"); # (name of slurm job array batch script (must match in sbatch below),  prefix of indexed jobs for job array)
		my $jobmem = $slurmmem*8;
		$jobarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude  --mem=$jobmem -n 1 -t $slurmtime --array=@rerun --job-name=$jobid.REGLIM --wrap=\"./ReRunJobArrays.sh\" | awk ' { print \$4 }'`;
		chomp $jobarrayid;
		print STDERR "Submitted job $jobarrayid - to rerun putative protein extraction jobs @rerun\n";
	} else {
		write_suffix_array_slurm_script("ReRunJobArrays.sh", "JobArray"); # (name of slurm job array batch script (must match in sbatch below),  prefix of indexed jobs for job array)
		$jobarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude  --mem=$slurmmem -n 1 -t $slurmtime --array=[@rerun] --job-name=$jobid.REGLIM --wrap=\"./ReRunJobArrays.sh\" | awk ' { print \$4 }'`;
		chomp $jobarrayid;
		print STDERR "Submitted job $jobarrayid - to rerun putative protein extraction jobs @rerun\n";
	}
	my $clustercheckpointid=`sbatch -d afterok:$jobarrayid --exclude=$slurmexclude --mem=200 -n 1 -t 10 --job-name=$jobid\".\"RECLUSCHK -p $slurmqueue --wrap=\"$script_dir/clustering_checkpoint.pl $count $jobid $assembly $refnum $reference\"| awk ' { print \$4 }'`;
	chomp $clustercheckpointid;
	print STDERR "Resubmitted job $clustercheckpointid - to check clusters\n";

} else {
	# concatenate the assembly reports
	if (-e "$dir/all.strains.assembly_stats.out") {
		system "rm $dir/all.strains.assembly_stats.out";
	}

	system "echo -e \"File\tLength\t#_contigs\tN50\t#_CDS\" > all.strains.assembly_stats.out";

	# JH filenames were wrong, .report files are based on fa filenames, used to use the fastq.gz filenames output at start of master script
	foreach my $num (keys %tracking) {
		unless ($num == $refnum) {
			my $filename = $tracking{$num};
			$filename =~ s/_1.fastq.gz/.fa/g;
			system "cat $filename.report >> all.strains.assembly_stats.out";
			system "rm $filename.report";
		}
	}

	# submit the concatenation script
	my $concatjobid=`sbatch -n 1 --exclude=$slurmexclude --mem=200 -t 10 -p $slurmqueue --job-name=${jobid}.CAT --wrap=\"./concatenation_script.sh\" | awk ' { print \$4 }'`;
	chomp $concatjobid;
	print STDERR "Submitted job $concatjobid - concatenation script\n";
	
	# run BLAST comparisons
	write_suffix_array_slurm_script("RunBlastArrays.sh", "BlastJob"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
	my $blastarrayid=`sbatch -d afterok:$concatjobid -p $slurmqueue --exclude=$slurmexclude --mem=$slurmmem -n 1 -t $slurmtime --array=1-$count --job-name=$jobid.UBLA --wrap=\"./RunBlastArrays.sh" | awk ' { print \$4 }'`;
	chomp $blastarrayid;
	print STDERR "Submitted job $blastarrayid - to run Blast\n";
	
	write_suffix_array_slurm_script("RunFilterBlastArrays.sh", "FilterBlast"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
	my $filterblastarrayid=`sbatch -d afterok:$concatjobid -p $slurmqueue --exclude=$slurmexclude --mem=$slurmmem -n 1 -t $slurmtime --array=1-$count --job-name=$jobid.FBLA --wrap=\"./RunFilterBlastArrays.sh" | awk ' { print \$4 }'`;
	chomp $filterblastarrayid;
	print STDERR "Submitted job $filterblastarrayid - to run Filtered Blast\n";
	
	# run Cogtriangles on the BLAT comparisons

	# calculate approximate memory requirements - 13 Gb per 100 strains, 20 Gb starting level - this gives the 98 Gb needed to process the SPARC collection, need to fix this for other bugs

	my @blat_array = (1..$count);
	my @filter_array = (1..$count);
	$" = ",";


	my $blast_checkpointid=`sbatch -d afterok:$blastarrayid:$filterblastarrayid -n 1 --exclude=$slurmexclude --mem=200 -t 10  -p $slurmqueue --job-name=${jobid}.BLACHK --wrap=\"$script_dir/blast_checkpoint.pl @blat_array @filter_array $count $jobid $reference $refnum\" | awk ' { print \$4 }'`;
	chomp $blast_checkpointid;
		
	foreach my $num (@file_nums) {
		my $file_dna = $tracking{$num};
		unless ($file_dna =~ /.dna/) {
			$file_dna =~ s/_1.fastq$|_1.fastq.gz$|.dna|.fasta|.seq|.fa/.dna/g; #JH mod to remove correct file if input is gzipped
			system "rm $file_dna";
		}
		my $file_fa = $tracking{$num};
		unless ($file_fa =~ /.fa/) {
			$file_fa =~ s/_1.fastq$|_1.fastq.gz$|.dna|.fasta|.seq|.fa/.fa/g; #JH mod to remove correct file if input is gzipped
			system "rm $file_dna";
			system "rm $file_fa";
		}
	}

	my $refstem = $tracking{$refnum};
	$refstem =~ s/_1.fastq$|_1.fastq.gz$|.dna|.fasta|.seq|.fa//g; #JH mod to remove correct file if input is gzipped
	#JH system "rm $file_dna";
	system "rm $tracking{$refnum}.report $refstem.train $refstem.longorfs reference_training.sh ";
}

print STDERR "Completed clustering_checkpoint.pl\n";


##############
# SUBROUTINES #
##############


sub write_suffix_array_slurm_script {
	my $scriptname=$_[0];
	my $stageprefix=$_[1];
	open SLURMSCRIPTFH, "> $scriptname";
	print SLURMSCRIPTFH "#!/bin/bash\n";
	print SLURMSCRIPTFH "./$stageprefix.";
	print SLURMSCRIPTFH '$SLURM_ARRAY_TASK_ID';
	close SLURMSCRIPTFH;
	system "chmod +x $scriptname";
}