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
my $jobmem = $slurmmem*8;
my $concatjobid=`sbatch -n 1 --mem=$jobmem -t $slurmtime --open-mode=append -o all.log -e all.err -p $slurmqueue --job-name=${jobid}.CAT --wrap=\"./concatenation_script.sh\" | awk ' { print \$4 }'`;
chomp $concatjobid;
print STDERR "Submitted job $concatjobid - concatenation script\n";
#system "bsub -J \"$jobid"."CAT\" -o all.log -e all.err -M 8000000 -R 'select[mem>8000] rusage[mem=8000]' ./concatenation_script.sh";

# run BLAST comparisons
write_suffix_array_slurm_script("RunBlastArrays.sh", "BlastJob"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
my $blastarrayid=`sbatch -d afterok:$concatjobid -p $slurmqueue --mem=$slurmmem -n 1 -t $slurmtime --array=1-$count --job-name=$jobid.UBLA --wrap=\"./RunBlastArrays.sh" | awk ' { print \$4 }'`;
chomp $blastarrayid;
print STDERR "Submitted job $blastarrayid - to run Blast\n";
#system "bsub -w \"ended($jobid"."CAT)\" -M 2000000 -R 'select[mem>2000] rusage[mem=2000] hname!=bc-11-4-10' -J $jobid"."UBLA[1-$count]".' -q normal_serial -o log.%I -e err.%I ./BlastJob.\${LSB_JOBINDEX}';

write_suffix_array_slurm_script("RunFilterBlastArrays.sh", "FilterBlast"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
my $filterblastarrayid=`sbatch -d afterok:$concatjobid -p $slurmqueue --mem=$slurmmem -n 1 -t $slurmtime --array=1-$count --job-name=$jobid.FBLA --wrap=\"./RunFilterBlastArrays.sh" | awk ' { print \$4 }'`;
chomp $filterblastarrayid;
print STDERR "Submitted job $filterblastarrayid - to run Filtered Blast\n";
#system "bsub -w \"ended($jobid"."CAT)\" -M 2000000 -R 'select[mem>2000] rusage[mem=2000] hname!=bc-11-4-10' -J $jobid"."FBLA[1-$count]".' -q normal_serial -o log.%I -e err.%I ./FilterBlast.\${LSB_JOBINDEX}';


# run Cogtriangles on the BLAT comparisons

# calculate approximate memory requirements - 13 Gb per 100 strains, 20 Gb starting level - this gives the 98 Gb needed to process the SPARC collection, need to fix this for other bugs

my @blat_array = (1..$count);
my @filter_array = (1..$count);
$" = ",";

my $blast_checkpointid=`sbatch -d afterok:$blastarrayid:$filterblastarrayid -n 1 --mem=$slurmmem -t $slurmtime --open-mode=append -o all.log -e all.err -p $slurmqueue --job-name=${jobid}.BLACHK --wrap=\"$script_dir/blast_checkpoint.pl @blat_array @filter_array $count $jobid $reference $refnum\" | awk ' { print \$4 }'`;
chomp $blast_checkpointid;
	

	
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


print STDERR "Completed clustering_checkpoint.pl\n";


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