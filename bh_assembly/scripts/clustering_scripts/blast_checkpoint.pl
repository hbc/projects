#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by clustering_checkpoint.pl | Submits failed BLAT jobs, COGtriangles  #
# runs and extract_annotate_orthologues.pl                                        #
# Checks whether BLAT jobs have completed successfully                            #
###################################################################################

# set up script_dir variable

my $script_dir = $ENV{SCRIPT_DIR};
my $slurmqueue = $ENV{SLURMQUEUE};

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
		if ($filter_output{$key} ne "0") {
			push(@filter_rerun,$key);
		}
	}
}

foreach my $key (@blat_array) {
	unless ($key == 0) {
		if ($blast_output{$key} ne "0") {
			push(@blast_rerun,$key);
		}
	}
}

# if BLAT jobs completed, then submit clustering jobs; else rerun BLATs

if ($#filter_rerun == -1 && $#blast_rerun == -1) {
	
	my $bignum = 10000000;
	my $smallnum = 10000;
	my $queue = "long_serial";

	if ($count > 100) {
		$bignum = (int($count/100)*13000000)+20000000;
		$smallnum = $bignum/1000;
		if ($smallnum > 20000) {
			$queue = "bigmem";
		}
	}
	
	
	my $cogAid=`sbatch -n 1 --mem=16000 -t 1-0 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.cogA --wrap=\"./cogtriangle_run_A.sh\" | awk ' { print \$4 }'`;
	chomp $cogAid;
	#system "bsub -o all.log -e all.err -q long_serial -J $jobid.cogA -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' ./cogtriangle_run_A.sh";

	my $cogBid=`sbatch -d afterok:$cogAid -n 1 --mem=$smallnum -t 1-0 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.cogB --wrap=\"./cogtriangle_run_B.sh\" | awk ' { print \$4 }'`;
	chomp $cogBid;
	#system "bsub -w \"ended($jobid.cogA)\" -o all.log -e all.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' -J $jobid.cogB ./cogtriangle_run_B.sh";
	
	my $cogCid=`sbatch -d afterok:$cogBid -n 1 --mem=16000 -t 1-0 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.cogC --wrap=\"./cogtriangle_run_C.sh\" | awk ' { print \$4 }'`;
	chomp $cogCid;
	#system "bsub -w \"ended($jobid.cogB)\" -o all.log -e all.err -q long_serial -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid.cogC ./cogtriangle_run_C.sh";
	
	my $postprocesscogsid=`sbatch -d afterok:$cogCid -n 1 --mem=16000 -t 1-0 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.cluchk --wrap=\"$script_dir/post_process_COGs.pl all.strains.cls.out.csv ./blaf/filtered.all.strains.blast.tab all.strains.csv $refnum\" | awk ' { print \$4 }'`;
	chomp $postprocesscogsid;
	#system "bsub -w \"ended($jobid.cogC)\" -o all.log -e all.err -q long_serial -M 16000000 -R 'select[mem>16000] rusage[mem=16000]' -J $jobid.cluchk $script_dir/post_process_COGs.pl all.strains.cls.out.csv ./blaf/filtered.all.strains.blast.tab all.strains.csv $refnum";
	
	my $embl = $reference;
	$embl =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.embl/g;
	
	my $extractannotateorthologsid=`sbatch -d afterok:$postprocesscogsid -n 1 --mem=$smallnum -t 1-0 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.EAO --wrap=\"$script_dir/extract_annotate_orthologues.pl $count $refnum $embl\" | awk ' { print \$4 }'`;
	chomp $extractannotateorthologsid;
	#system "bsub -o all.log -e all.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' -w \"ended($jobid.cluchk)\" $script_dir/extract_annotate_orthologues.pl $count $refnum $embl";
	
} else {
	
	my $new_id = int(rand(1000));
	
	my $dependency_string;
	
	$" = ",";
	
	if ($#filter_rerun > -1) {
		
		print STDERR "bsub -M 2000000 -R 'select[mem>2000] rusage[mem=2000]' -J $new_id"."FBLA[@filter_rerun]".' -q normal_serial -o log.%I -e err.%I ./FilterBlast.\${LSB_JOBINDEX}'."\n";
		system "bsub -M 2000000 -R 'select[mem>2000] rusage[mem=2000]' -J $new_id"."FBLA[@filter_rerun]".' -q normal_serial -o log.%I -e err.%I ./FilterBlast.\${LSB_JOBINDEX}';
	
		$dependency_string = "ended($new_id"."FBLA)";
		
	} else {
		push(@filter_rerun,"0");
	}
	
	if ($#blast_rerun > -1) {
		
		system "bsub -M 2000000 -R 'select[mem>2000] rusage[mem=2000]' -J $new_id"."UBLA[@blast_rerun]".' -q normal_serial -o log.%I -e err.%I ./BlastJob.\${LSB_JOBINDEX}';
		
		if (defined($dependency_string) && length($dependency_string) > 0) {
			$dependency_string.=" && ";
		}
		
		$dependency_string.="ended($new_id"."UBLA)";
	} else {
		push(@blast_rerun,"0");
	}
	
	print STDERR "bsub -o all.log -e all.err -w \"$dependency_string\" $script_dir/blast_checkpoint.pl @blast_rerun @filter_rerun $count $jobid $reference $refnum\n";
	system "bsub -o all.log -e all.err -w \"$dependency_string\" $script_dir/blast_checkpoint.pl @blast_rerun @filter_rerun $count $new_id $reference $refnum";
}


print STDERR "Completed blast_checkpoint.pl\n";
