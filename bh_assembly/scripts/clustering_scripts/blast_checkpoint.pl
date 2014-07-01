#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by clustering_checkpoint.pl | Submits failed BLAT jobs, COGtriangles  #
# runs and extract_annotate_orthologues.pl                                        #
# Checks whether BLAT jobs have completed successfully                            #
###################################################################################

# set up script_dir and slurm variables
my $script_dir = $ENV{SCRIPT_DIR};
my $slurmqueue = $ENV{SLURMQUEUE};
my $slurmtime = $ENV{SLURMTIME};
my $slurmmem = $ENV{SLURMMEM};


# set up other variables

my @blat_array = split(/\,/,$ARGV[0]);
my @filter_array = split(/\,/,$ARGV[1]);
my $count = $ARGV[2];
my $jobid = $ARGV[3];
my $reference = $ARGV[4];
my $refnum = $ARGV[5];
my @blast_rerun;
my @filter_rerun;
my %filter_output;
my %blast_output;

# read blast output file

open IN, 'blast_outputs.log' or die print STDERR "Could not open blast_outputs.log!\n";

	foreach (<IN>) {
	chomp;
	my @data = split(/\s+/,$_);
	if ($data[0] eq 'FilterBlast') {
		$filter_output{$data[1]} = $data[2];
	} elsif ($data[0] eq 'BlastJob') {
		$blast_output{$data[1]} = $data[2];
	}
}

close IN;

system "mv blast_outputs.log previous.blast_outputs.log";

# summarise hashes into arrays of jobs that need rerunning

foreach my $key (@filter_array) {
	unless ($key == 0) {
		if (!defined $filter_output{$key}){
			push(@filter_rerun,$key);
		} else {
			if ($filter_output{$key} ne "0") {
				push(@filter_rerun,$key);
			}
		}
	}
}

foreach my $key (@blat_array) {
	unless ($key == 0) {
		if (!defined $blast_output{$key}){
			push(@blast_rerun,$key);
		} else {
			if ($blast_output{$key} ne "0") {
				push(@blast_rerun,$key);
			}
		}
	}
}
# if BLAT jobs completed, then submit clustering jobs; else rerun BLATs


if ($#filter_rerun == -1 && $#blast_rerun == -1) {
	
	my $cogAid=`sbatch -n 1 --mem=200 -t 10 -p $slurmqueue --job-name=$jobid.COGA --wrap=\"./cogtriangle_run_A.sh\" | awk ' { print \$4 }'`;
	chomp $cogAid;
	#system "bsub -o all.log -e all.err -q long_serial -J $jobid.cogA -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' ./cogtriangle_run_A.sh";
	print STDERR "Submitted cogtriangle-A job $cogAid\n";

	my $cogBid=`sbatch -d afterok:$cogAid -n 1 --mem=200 -t 10 -p $slurmqueue --job-name=$jobid.COGB --wrap=\"./cogtriangle_run_B.sh\" | awk ' { print \$4 }'`;
	chomp $cogBid;
	#system "bsub -w \"ended($jobid.cogA)\" -o all.log -e all.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' -J $jobid.cogB ./cogtriangle_run_B.sh";
	print STDERR "Submitted cogtriangle-B job $cogBid\n";


	my $cogCid=`sbatch -d afterok:$cogBid -n 1 --mem=$slurmmem -t 10 -p $slurmqueue --job-name=$jobid.COGC --wrap=\"./cogtriangle_run_C.sh\" | awk ' { print \$4 }'`;
	chomp $cogCid;
	#system "bsub -w \"ended($jobid.cogB)\" -o all.log -e all.err -q long_serial -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid.cogC ./cogtriangle_run_C.sh";
	print STDERR "Submitted cogtriangle-C job $cogCid\n";


	my $postprocesscogsid=`sbatch -d afterok:$cogCid -n 1 --mem=200 -t 10 -p $slurmqueue --job-name=$jobid.PPCOGS --wrap=\"$script_dir/post_process_COGs.pl all.strains.cls.out.csv ./blaf/filtered.all.strains.blast.tab all.strains.csv $refnum\" | awk ' { print \$4 }'`;
	chomp $postprocesscogsid;
	#system "bsub -w \"ended($jobid.cogC)\" -o all.log -e all.err -q long_serial -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid.cluchk $script_dir/post_process_COGs.pl all.strains.cls.out.csv ./blaf/filtered.all.strains.blast.tab all.strains.csv $refnum";
	print STDERR "Submitted COG post-processing job $postprocesscogsid\n";


	my $embl = $reference;
	$embl =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.embl/g;
	
	my $extractannotateorthologsid=`sbatch -d afterok:$postprocesscogsid -n 1 --mem=200 -t 10  -p $slurmqueue --job-name=$jobid.EXANORTH --wrap=\"$script_dir/extract_annotate_orthologues.pl $count $refnum $embl\" | awk ' { print \$4 }'`;
	chomp $extractannotateorthologsid;
	#system "bsub -o all.log -e all.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' -w \"ended($jobid.cluchk)\" $script_dir/extract_annotate_orthologues.pl $count $refnum $embl";
	print STDERR "Submitted Ortholog extraction and annotation job $extractannotateorthologsid \n";


} else {
	
	my $dependency_string;
	
	$" = ",";
	
	if ($#filter_rerun > -1) {

		write_suffix_array_slurm_script("ReRunFilterBlastArrays.sh", "FilterBlast"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
		my $filterblastarrayid=`sbatch -p $slurmqueue --mem=$slurmmem -n 1 -t $slurmtime --array=@filter_rerun --job-name=$jobid.REFBLA --wrap=\"./ReRunFilterBlastArrays.sh" | awk ' { print \$4 }'`;
		chomp $filterblastarrayid;
		print STDERR "Submitted job $filterblastarrayid - to rerun Filtered Blasts @filter_rerun\n";
		$dependency_string = $filterblastarrayid;

	} else {
		push(@filter_rerun,"0");
	}
	
	if ($#blast_rerun > -1) {
					
		write_suffix_array_slurm_script("ReRunBlastArrays.sh", "BlastJob"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
		my $blastarrayid=`sbatch -p $slurmqueue --mem=$slurmmem -n 1 -t $slurmtime --array=@blast_rerun --job-name=$jobid.REUBLA --wrap=\"./ReRunBlastArrays.sh" | awk ' { print \$4 }'`;
		chomp $blastarrayid;
		print STDERR "Submitted job $blastarrayid - to rerun Blasts @blast_rerun\n";


		if (defined($filterblastarrayid) && length($filterblastarrayid) > 0) {
			$dependency_string="$filterblastarrayid:$blastarrayid";
		} else {
			$dependency_string = $blastarrayid;
		}
				
	} else {
		push(@blast_rerun,"0");
	}
	
	my $blast_checkpointid=`sbatch -d afterok:$dependency_string -n 1 --mem=200 -t 10  -p $slurmqueue --job-name=${jobid}.REBLACHK --wrap=\"$script_dir/blast_checkpoint.pl @blast_rerun @filter_rerun $count $jobid $reference $refnum\" | awk ' { print \$4 }'`;
	chomp $blast_checkpointid;

}


print STDERR "Completed blast_checkpoint.pl\n";


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
