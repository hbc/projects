#! /usr/bin/perl -w
use strict;
use Bio::SeqIO;

###################################################################################
# Submitted by clustering_checkpoint.pl | Submits check_COG_analysis.pl and       #
# align_clusters.pl                                                               #
# Processes clusters before analysis, identifies the core genome                  #
###################################################################################

# set up script_dir and slurm variables
my $script_dir = $ENV{SCRIPT_DIR};
my $slurmqueue = $ENV{SLURMQUEUE};
my $slurmtime = $ENV{SLURMTIME};
my $slurmmem = $ENV{SLURMMEM};
my $slurmexclude = $ENV{SLURMEXCLUDE};



# random number for job submission/dependencies

my $range = 1000;
my $jobid = int(rand($range));

my $count = $ARGV[0];
my $refnum = $ARGV[1];
my $refstring = "seq"."$refnum";
my $reference = $ARGV[2];

# process COGtriangles output to identify single copy core CDSs

my %clusters;
my %clusters_strains;
my %cds_cluster;
my $clucount;

my %SCcore;
my $SCcount = 1;
my %refvers;

# read in clustering files for reference and query strains

print STDERR "Reading in clusters\n";

my @cluster_files = ("all.strains.filtered.cls.csv","reference.strain.filtered.cls.csv");

foreach my $file (@cluster_files) {
	open CLU, "$file" or die print STDERR "Cannot open cluster file $file!\n";
	while (my $line = <CLU>) {
		chomp $line;
		my @data = split(/\,/,$line);
		if (defined($data[1]) && $data[1] =~ /^CLS/) {
			my @strain = split(/orf/,$data[0]);
			unless ((defined($cds_cluster{$data[0]}{$data[1]}) && $cds_cluster{$data[0]}{$data[1]} == 1)) {
				if ($strain[0] eq $refstring) {
					$refvers{$data[1]} = $data[0];
				} else {
					push(@{$clusters{$data[1]}},$data[0]);
					if (defined($clusters_strains{$data[1]}{$strain[0]})) {
						$clusters_strains{$data[1]}{$strain[0]}++;
					} else {
						$clusters_strains{$data[1]}{$strain[0]} = 1;
					}
				}
			}
			$cds_cluster{$data[0]}{$data[1]} = 1;
		}
	}
	close CLU;
}

print STDERR "Clusters read; now printing presence/absence matrix\n";

open PAM, "> presence_absence_matrix.out" or die print STDERR "Could not open presence absence matrix for writing!\n";

print PAM "Taxon";

foreach my $cluster_name (sort keys %clusters_strains) {
	print PAM "\t$cluster_name";
}

print PAM "\n";

my $num = 1;

while ($num <= $count) {
	unless ($num == $refnum) {
		my $seqfull = "seq".$num;
		print PAM "$seqfull";
		#JH foreach my $cluster (sort keys %{%clusters_strains}) { #JH remove extra curly braces
		foreach my $cluster (sort keys %clusters_strains) { 
			if (defined($clusters_strains{$cluster}{$seqfull}) && $clusters_strains{$cluster}{$seqfull} >= 1) {
				print PAM "\t1";
			} else {
				print PAM "\t0";
			}
		}
		print PAM "\n";
	}
	$num++;
}

close PAM;

print STDERR "Printed; now submitting pangenome calculation job\n";

my $running_total = 0;

# run pangenome analysis

my $bignum = 10000000;
my $smallnum = 10000;
my $queue = "long_serial";

if ($count > 100) {
	$bignum = int($count/100)*10000000;
	$smallnum = $bignum/1000;
}

my $ccpgid=`sbatch -n 1 --exclude=$slurmexclude  --mem=$smallnum -t 60 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.CALCCOGPAN --wrap=\"$script_dir/calculate_COG_pangenome.pl strain.info all.strains.cls.out.csv $refnum\" | awk ' { print \$4 }'`;
chomp $ccpgid;
#JH system "bsub -o all.log -e log.err -q $queue -M $bignum -R 'select[mem>$smallnum] rusage[mem=$smallnum]' $script_dir/calculate_COG_pangenome.pl strain.info all.strains.cls.out.csv $refnum";

print STDERR "Submitted; now identifying single copy orthologues\n";

# identify single copy orthologue groups

my %core;
my %SClist;

foreach $clucount (sort keys %clusters_strains) {
	if (scalar(keys %{$clusters_strains{$clucount}}) == ($count-1)) {		# checks there is at least one entry per strain
		foreach my $strain (sort keys %{$clusters_strains{$clucount}}) {	# checks there is only one entry per strain
			$running_total+=$clusters_strains{$clucount}{$strain};
		}
		if ($running_total == ($count-1)) {					# total count excluding the reference
			if (defined($refvers{$clucount}) && length($refvers{$clucount}) > 0) {
				$core{$SCcount} = $refvers{$clucount};
			} else {
				$core{$SCcount} = 0;
			}
		} else {
			$core{$SCcount} = 0;
		}
	} else {
		$core{$SCcount} = 0;
	}
	@{$SCcore{$SCcount}} = @{$clusters{$clucount}};
	$SClist{$SCcount} = $clucount;
	$SCcount++;
	$running_total = 0;
}

my $max_SCcount = $SCcount-1;

# identify corresponding protein sequences and align with MUSCLE

$num = 1;

print STDERR "Done; now submitting alignment jobs\n";

my $normal_count = 1;
my $small_count = 1;
$refstring.="orf";

# subdividing folders
my $ten_count = 1;
my $thousand_count = 1;
my $base_dir  = `pwd | tr -d "\n"`;
my $working_dir = $base_dir."/aligndir_".$thousand_count;
system "mkdir $working_dir; cd $working_dir; ln -s ../protein_sequences.db .; ln -s ../DNA_sequences.db .; ln -s ../$reference .; cd ..;";

open SMA, "> SmallAlignArray.$small_count";

foreach $SCcount (sort keys %SCcore) {						# submit analysis jobs for clusters to different queues depending on the number of members they contain
	if (($thousand_count/1000)-(int($thousand_count/1000)) == 0) {
		my $work_dir_num = ($thousand_count/1000)+1;
		$thousand_count = 1;
		$working_dir = $base_dir."/aligndir_".$work_dir_num;
		system "mkdir $working_dir; cd $working_dir; ln -s ../protein_sequences.db .; ln -s ../DNA_sequences.db .; ln -s ../$reference .; cd ..;";
	}
	if (scalar(@{$SCcore{$SCcount}}) >= 75) {
		open OUT, "> AlignArray.$normal_count";
		print OUT "$script_dir/align_cluster.pl $reference $refstring $core{$SCcount} $SCcount $SClist{$SCcount} $working_dir @{$SCcore{$SCcount}}\n";
		close OUT;
		system "chmod +x AlignArray.$normal_count";
		$normal_count++;
		print $normal_count;
	} elsif (scalar(@{$SCcore{$SCcount}}) > 0) {
		if (($ten_count/10)-(int($ten_count/10)) == 0) {
			close SMA;
			system "chmod +x SmallAlignArray.$small_count";
			$small_count++;
			$ten_count = 1;
			open SMA, "> SmallAlignArray.$small_count";
		}
		print SMA "$script_dir/align_cluster.pl $reference $refstring $core{$SCcount} $SCcount $SClist{$SCcount} $working_dir @{$SCcore{$SCcount}}\n";
		$ten_count++;
	}
	$thousand_count++;
}

close SMA;
system "chmod +x SmallAlignArray.$small_count";

if ($ten_count != 1) {
	$small_count++;
}

$normal_count--;
$small_count--;
my $dependency_string;
my $alignarrayid;
my $smallalignarrayid;

# run analysis of clusters containing more than 74 members

if ($normal_count >= 1) { #JH changed from >1 to >=1
	write_suffix_array_slurm_script("AlignArrays.sh", "AlignArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
	my $alignarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude  --mem=2500 -n 1 -t 60 --array=1-$normal_count --job-name=$jobid.ALN --wrap=\"./AlignArrays.sh\" | awk ' { print \$4 }'`;
	chomp $alignarrayid;
	#JH system "bsub -M 2500000 -R 'select[mem>2500] rusage[mem=2500]' -q normal_serial -J $jobid"."aln[1-$normal_count]".'%100 -o log.%I -e err.%I ./AlignArray.\$LSB_JOBINDEX';
	#JH $dependency_string = "ended($jobid"."aln[1-$normal_count])";
	$dependency_string="$alignarrayid"; #JH
}


# run analysis of clusters containing fewer than 75 members

if ($small_count >= 1) { #JH changed from >1 to >=1
	write_suffix_array_slurm_script("SmallAlignArrays.sh", "SmallAlignArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
	my $smallalignarrayid=`sbatch -p $slurmqueue --exclude=$slurmexclude  --mem=2500 -n 1 -t 60 --array=1-$small_count --job-name=$jobid.SMALN --wrap=\"./SmallAlignArrays.sh\" | awk ' { print \$4 }'`;
	chomp $smallalignarrayid;
	$dependency_string="$smallalignarrayid"; #JH

	#JHsystem "bsub -q short_serial -M 2500000 -R 'select[mem>2500] rusage[mem=2500]' -J $jobid"."Saln[1-$small_count]".'%100 -o Slog.%I -e Serr.%I ./SmallAlignArray.\$LSB_JOBINDEX';
	if (defined($alignarrayid) && length($alignarrayid) > 0) {
		$dependency_string="$alignarrayid:$smallalignarrayid";
	
	#JH if (defined($dependency_string) && length($dependency_string) > 0) {
	#JH	$dependency_string.=" && ";
	}
	#JH $dependency_string.="ended($jobid"."Saln[1-$small_count])";
}

# submit checkpointing and alignment processing jobs

my $ccaid=`sbatch -d afterok:$dependency_string -n 1 --exclude=$slurmexclude  --mem=2500 -t 60 -o all.log -e all.err -p $slurmqueue --job-name=$jobid.CHKCOG --wrap=\"$script_dir/check_COG_analysis.pl $normal_count $small_count $count $reference $refnum\" | awk ' { print \$4 }'`;
chomp $ccaid;
#JH system "bsub -o all.log -e all.err -M 2500000 -R 'select[mem>2500] rusage[mem=2500]' -w \"$dependency_string\" -J $jobid"."COGCHK $script_dir/check_COG_analysis.pl $normal_count $small_count $count $reference $refnum;

print STDERR "Done\nCompleted extract_annotate_orthologues\n";


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
