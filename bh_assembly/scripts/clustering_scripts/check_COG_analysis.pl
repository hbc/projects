#! /usr/bin/perl -w
use strict;

###################################################################################
# Submitted by extract_annotate_orthologues.pl | Submits align_clusters.pl or     #
# orthologue_alignment_processing.pl                                              #
# Checks whether the alignment script has successfully run on all COGs            #
###################################################################################

# set up script_dir and slurm variables
my $script_dir = $ENV{SCRIPT_DIR};
my $slurmqueue = $ENV{SLURMQUEUE};
my $slurmtime = $ENV{SLURMTIME};
my $slurmmem = $ENV{SLURMMEM};
my $slurmexclude = $ENV{SLURMEXCLUDE};
my $dir = `pwd | tr -d "\n"`;

# set up variables

my $range = 1000;
my $jobid = int(rand($range));

my @jobarray;
my $max_jobarray = $ARGV[0];
my @smalljobarray;
my $max_smalljobarray = $ARGV[1];
my $count = $ARGV[2];
my $reference = $ARGV[3];
my $refnum = $ARGV[4];
my $n = 1;

# create filname variable 

# check whether all normal queue jobs have completed; rerun if not, delete files if they have

while ($n <= $max_jobarray) {
	if (-e "AlignArray.$n") {
		open IN, "AlignArray.$n" or die print STDERR "Cannot open file AlignArray.$n!\n";
		foreach (<IN>) {
			chomp;
			my @contents = split(/\s+/,$_);
			my @cluster_files = ("$contents[6]/cluster."."$contents[4]".".aln","$contents[6]/cluster."."$contents[4]".".seq","$contents[6]/cluster."."$contents[4]".".DNA.seq");
			if (-e "$contents[6]/cluster.$contents[4].report") {					# check if the report has been generated
				if ($contents[3] ne "0") {
					if (-e "$contents[6]/cluster.$contents[4].DNA.aln" && -e "$contents[6]/partial.$contents[4].core_genome.tab") {								# check if core COGs are aligned
						clean_up(@cluster_files);
					} else {
						push(@jobarray,$_);
					}
				} else {
					push(@cluster_files,"$contents[6]/cluster."."$contents[4]".".DNA.aln");
					clean_up(@cluster_files);
				}
			} else {
				push(@jobarray,$_);
			}
		}
		close IN;
	} else {
		die print STDERR "Cannot open AlignArray.$n!\n";
	}
	$n++;
}

# same process for all small queue jobs

$n = 1;

while ($n <= $max_smalljobarray) {
	if (-e "SmallAlignArray.$n") {
		open IN, "SmallAlignArray.$n" or die print STDERR "Cannot open file SmallAlignArray.$n!\n";
		foreach (<IN>) {
			chomp;
			my @contents = split(/\s+/,$_);
			my @cluster_files = ("$contents[6]/cluster."."$contents[4]".".aln","$contents[6]/cluster."."$contents[4]".".seq","$contents[6]/cluster."."$contents[4]".".DNA.seq");
			if (-e "$contents[6]/cluster.$contents[4].report") {					# check if the report has been generated
				if ($contents[3] ne "0") {
					if (-e "$contents[6]/cluster.$contents[4].DNA.aln") {								# check if core COGs are aligned
						clean_up(@cluster_files);
					} else {
						push(@smalljobarray,$_);
					}
				} else {
					push(@cluster_files,"$contents[6]/cluster."."$contents[4]".".DNA.aln");
					clean_up(@cluster_files);
				}
			} else {
				push(@smalljobarray,$_);
			}
		}
		close IN;
	} else {
		die print STDERR "Cannot open SmallAlignArray.$n!\n";
	}
	$n++;
}

# remove old job array files

system "rm AlignArray.* SmallAlignArray.*";

# rerun alignment jobs if necessary; otherwise submit alignment processing job

$" = ", ";

my $dependency_string;
my $jobarray_n = 0;
my $smalljobarray_n = 0;
my $alignarrayid;
my $smallalignarrayid;

if ($#jobarray > -1) {
	$jobarray_n = 1;
	foreach my $jobtext (@jobarray) {		# sort into ascending order to avoid larger num files overwriting smaller num files
		open ARR, "> AlignArray.$jobarray_n" or die print STDERR "Cannot open file AlignArray.$jobarray_n!\n";
		print ARR "$jobtext\n";
		close ARR;
		system "chmod +x AlignArray.$jobarray_n";
		$jobarray_n++;
	}
	$jobarray_n--;
	if ($jobarray_n > 1) {
		write_suffix_array_slurm_script("RAlignArrays.sh", "AlignArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
		my $alignarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude --mem=$slurmmem -n 1 -t 60 --array=1-$jobarray_n --job-name=$jobid.RALN --wrap=\"./RAlignArrays.sh\" | awk ' { print \$4 }'`;
		chomp $alignarrayid;
	#JH	system "bsub -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -q normal_serial -J $jobid"."Raln[1-$jobarray_n]".' -o log.%I -e err.%I ./AlignArray.\$LSB_JOBINDEX';
	} else {
		write_suffix_array_slurm_script("RAlignArrays.sh", "AlignArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
		my $alignarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude --mem=$slurmmem -n 1 -t 60 --array=1 --job-name=$jobid.RALN --wrap=\"./RAlignArrays.sh\" | awk ' { print \$4 }'`;
		chomp $alignarrayid;
	#JH system "bsub -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -q normal_serial -J $jobid"."Raln[1]".' -o log.%I -e err.%I ./AlignArray.\$LSB_JOBINDEX';
	}
	$dependency_string = $alignarrayid;
}

my $ten_count = 1;

if ($#smalljobarray > -1) {
	$smalljobarray_n = 1;
	foreach my $jobtext (@smalljobarray) {	# sort into ascending order to avoid larger num files overwriting smaller num files
		open ARR, ">> SmallAlignArray.$smalljobarray_n" or die print STDERR "Cannot open file SmallAlignArray.$smalljobarray_n!\n";
		print ARR "$jobtext\n";
		close ARR;
		system "chmod +x SmallAlignArray.$smalljobarray_n";
		$ten_count++;
		if ((($ten_count/10) - (int($ten_count/10))) == 0) {
			$smalljobarray_n++;
			$ten_count = 1;
		}
	}
	if ($ten_count == 1) {
		$smalljobarray_n--;
	}
	#JH  if (defined($dependency_string) && length($dependency_string) > 1) {
	# 	$dependency_string.=" && ";
	# }
	if ($smalljobarray_n > 1) {
		
		write_suffix_array_slurm_script("RSmallAlignArrays.sh", "SmallAlignArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
		my $smallalignarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude --mem=$slurmmem -n 1 -t 60 --array=1-$smalljobarray_n --job-name=$jobid.SRALN --wrap=\"./RSmallAlignArrays.sh\" | awk ' { print \$4 }'`;
		chomp $smallalignarrayid;
	 #JH system "bsub -M 4000000 -R 'select[mem>4000] rusage[mem=4000]' -q short_serial -J $jobid"."RSaln[1-$smalljobarray_n]".' -o Slog.%I -e Serr.%I ./SmallAlignArray.\$LSB_JOBINDEX';
		$dependency_string ="$smallalignarrayid";
	} else {
		write_suffix_array_slurm_script("RSmallAlignArrays.sh", "SmallAlignArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
		my $smallalignarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude --mem=$slurmmem -n 1 -t 60 --array=1 --job-name=$jobid.SRALN --wrap=\"./RSmallAlignArrays.sh\" | awk ' { print \$4 }'`;
		chomp $smallalignarrayid;
		#JH system "bsub -M 4000000 -R 'select[mem>4000] rusage[mem=4000]' -q normal_serial -J $jobid"."RSaln[1]".' -o Slog.%I -e Serr.%I ./SmallAlignArray.\$LSB_JOBINDEX';
		$dependency_string = $smallalignarrayid;
	}
	if (defined($alignarrayid) && length($alignarrayid) > 0) {
		$dependency_string="$alignarrayid:$smallalignarrayid";
	}
} #JH

if ($#jobarray == -1 && $#smalljobarray == -1) {
	my $jobmem = $slurmmem*20;	
	`sbatch -n 1 --exclude=$slurmexclude --mem=$jobmem -t $slurmtime -p $slurmqueue --job-name=${jobid}.OALPROC --wrap=\"$script_dir/orthologue_alignment_processing.pl $count $reference $refnum\"`;
	#JH	system "bsub -o all.log -e all.err -M 18000000 -R 'select[mem>18000] rusage[mem=18000]' $script_dir/orthologue_alignment_processing.pl $count $reference $refnum";
} else {
	`sbatch -d afterok:$dependency_string -n 1 --exclude=$slurmexclude --mem=$slurmmem -t $slurmtime -p $slurmqueue --job-name=${jobid}.COGCHK --wrap=\"$script_dir/check_COG_analysis.pl $jobarray_n $smalljobarray_n $count $reference $refnum\"`;
	#JH system "bsub -o all.log -e all.err -M 1000000 -R 'select[mem>1000] rusage[mem=1000]' -w \"$dependency_string\" -J $jobid"."COGCHK $script_dir/check_COG_analysis.pl $jobarray_n $smalljobarray_n $count $reference $refnum";	
}

print STDERR "Completed check_COG_analysis.pl\n";

##############
# SUBROUTINE #
##############

sub clean_up {
	my @cfiles = shift;
	foreach my $f (@cfiles) {
		if (-e "$f") {
			system "rm $f";
		}
	}
}

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

