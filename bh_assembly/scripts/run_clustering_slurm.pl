#! /usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
use DBI;

# set up scripts directory

my $script_dir = "/n/home08/jhutchin/consults/bh_assembly/scripts/clustering_scripts";
$ENV{SCRIPT_DIR} = $script_dir;

# random number for job submission/dependencies

my $range = 1000;
my $jobid = int(rand($range));
print STDERR "$jobid\n";

# deal with command line options

my $help = 0;
my $reference;

GetOptions (
	"r=s" => \$reference,
	"h" => \$help
);

if ($help == 1 || $#ARGV == -1) {
	croak();
}

unless (grep(/$reference/,@ARGV)) {
	push(@ARGV,$reference);
}

my $file;
my $count = 1;
my $countb = 1;
my %tracking;
my @csv_files;
my @aa_files;
my @dna_files;
my @filtered;
my $refnum;
my $refstem;

##########################################################################################
# CLEANUP

# remove any pre-existing analyses
print STDERR "\n\nRemoving previous analysis files\n";
my @previous = ("all.log","all.err","protein_sequences.db","DNA_sequences.db",
	"core_genome.tab","all.strains.csv","all.strains.aa*","filtered*aa","all.strains.blast",
	"seq*protseq","JobArray.*","log.*","err.*","seq.*.csv","RAxML*all.strains","shuffled*fastq",
	"all.strains.icm","all.strains.COG.report","reference.strain.filtered.cls.csv");
foreach $file (@previous) {
	if (-e $file) {
		system "rm $file";
	}
}
my @arrays = ("Output.*","SelfBlast.1","JobArray.1","log.1","err.1","seq.1.csv","GlimArray.1","AlignArray.1",
	"SmallAlignArray.1");
foreach $file (@arrays) {
	if (-e $file) {
		$file =~ s/1/\*/g;
		system "rm $file";
	}
}
my @folders = ('blan','blaf','tmp','BLASTconv','blastDB','pangenome','aligndir_*');
foreach my $folder (@folders) {
	if (-e $folder) {
		system "rm -rf $folder";
	}
}

my @folder_arrays = ('aligndir_1');
foreach my $folder (@folder_arrays) {
	if (-e $folder) {
		$folder =~ s/1/\*/g;
		system "rm -rf $folder";
	}
}

##########################################################################################
# SETUP

# Setup SQL databases
# create databse for protein sequences and connnect, fail safely if not
my $prot_dbfile = "protein_sequences.db";
my $prot_db = DBI->connect("dbi:SQLite:dbname=$prot_dbfile","","");
my $create = $prot_db->prepare("create table Proteins (Label text, Sequence text)") or die print STDERR "Could not prepare statement for ProtDB\n";
$create->execute or die print STDERR "Could not execute statement on ProtDB";
undef($create);
$prot_db->disconnect;
# create database for DNA sequences and connect, fail safely if not
my $DNA_dbfile = "DNA_sequences.db";
my $DNA_db = DBI->connect("dbi:SQLite:dbname=$DNA_dbfile","","");
$create = $DNA_db->prepare("create table DNA (Label text, Sequence text)");
$create->execute;
undef($create);
$DNA_db->disconnect;

# make EMBL of predicted genes, 
# extract protein sequences, 
# make CSV files for COGtriangles, 
# run BLAT comparisons for clustering

print STDERR "\nPreparing genomes for analysis\n";

open LSE, "> LSE.csv";
foreach $file (@ARGV) {
	# find the reverse paired read file if input forward paried read file (defined by presence of "_R1_" in filename)
	if ($file =~ /_R1_/) {
		print "Finding other paired read file\n";
		my @b = split(/_R1_/,$file);
		my $forward = $b[0]."_1.fastq";
		if ($file =~ /.gz/) {
			$forward.=".gz";
		}
		system "cp -s $file $forward";
		my $original_reverse = $file;
		$original_reverse =~ s/_R1_/_R2_/g;
		my $new_reverse = $forward;
		$new_reverse =~ s/_1.fastq/_2.fastq/g;
		if (-e "$original_reverse") {
			system "cp -s $original_reverse $new_reverse";
		} else {
			print STDERR "Cannot find reverse fastq $original_reverse\n";
			exit(0);
		}
		$file = $forward;
	}
	
	# printing LSE file for performing all pairwise comparisons in clustering
	foreach my $fileb (@ARGV) {
		unless ($file eq $fileb) {
			print LSE "seq"."$count,seq"."$countb\n";
		}
		$countb++;
	}
	$countb = 1;
	
	# make script for running gene prediction script
	open JA, "> JobArray.$count" or croak();
	print JA "$script_dir/extract_putative_proteins.pl $file $count $reference";
	close JA;
	my $outname = $file;
	$outname =~ s/_1.fastq|_1.fastq.gz|.dna|.fasta|.seq|.fa/.aa/g; #JH mod to give correct name with fastq.gz inputs
	my $DNA_outname = $file;
	$DNA_outname =~ s/_1.fastq|_1.fastq.gz|.dna|.fasta|.seq|.fa/.CDS.mfa/g; #JH mod to give correct name with fastq.gz inputs
	push(@aa_files,$outname);
	push(@csv_files,"seq.$count.csv");
	push(@dna_files,$DNA_outname);
	push(@filtered,"filtered.$outname");
	system "chmod +x JobArray.$count";
	
	# makeprint script for running filtered and unfiltered BLAT comparisons
	open BL, "> BlastJob.$count" or croak();
	print BL "blat blastDB/all.strains.aa $outname -prot -out=blast8 unfilteredblast.$count.out; echo -e \"BlastJob\t$count\t\$?\" >> blast_outputs.log";
	close BL;
	system "chmod +x BlastJob.$count";
	open BL, "> FilterBlast.$count" or croak();
	print BL "blat blastDB/filtered.all.strains.aa filtered.$outname -prot -out=blast8 filteredblast.$count.out; echo -e \"FilterBlast\t$count\t\$?\" >> blast_outputs.log";
	close BL;
	system "chmod +x FilterBlast.$count";
	
	# keep track of file numbering
		$tracking{$count} = $file;
	if ($file eq $reference) {
		$refnum = $count;
		$refstem = $outname;
		$refstem =~ s/.aa//g;
	}
	$count++;
}
close LSE;


# store information on strains in a file
open DATA, "> strain.info" or croak();
foreach my $num (sort keys %tracking) {
	print DATA "$num\t$tracking{$num}\n";
}
close DATA;


# write script for concatenating data
open CS, '> concatenation_script.sh' or croak();
#print CS "mkdir blastDB; lfs setstripe blastDB -c -1\n";
print CS "mkdir blastDB\n"; #JHremove failing lustre command
print CS "cat @aa_files > blastDB/all.strains.aa\ncat @csv_files > all.strains.csv\ncat @dna_files > blastDB/all.strains.dna\n";
#print CS "seg blastDB/all.strains.aa -n -x > blastDB/filtered.all.strains.aa\n";
print CS "cat @filtered > blastDB/filtered.all.strains.aa\n";
print CS "cd blastDB\n";
print CS "formatdb -i all.strains.aa -p T\n";
print CS "formatdb -i all.strains.dna -p F\n";
print CS "cd ..\n";
print CS "cat seq.*.protseq > seq.all.protseq\n";
print CS "sqlite3 protein_sequences.db \".import seq.all.protseq Proteins\"\n";
print CS "cat seq.*.DNAseq > seq.all.DNAseq\n";
print CS "sqlite3 DNA_sequences.db \".import seq.all.DNAseq DNA\"\n";
print CS "rm *.DNAseq *.protseq\n";
print CS "echo 'Completed concatenation_script.sh' >> all.err;";
close CS;
system "chmod +x concatenation_script.sh";


$count--;

# write script to train Glimmer and Prodigal on the reference sequence
print STDERR "\nSubmitting genome analysis jobs\n";

open REF, "> reference_training.sh";
print REF "
long-orfs -n -t 1.15 $reference $refstem.longorfs
extract -t $reference $refstem.longorfs > $refstem.train
build-icm -r all.strains.icm < $refstem.train
prodigal -t all.strains.prod.train < $reference
";
close REF;
system "chmod +x ./reference_training.sh";

# submit job and track
my $reftrainjobid=`sbatch -n 1 --mem=1000 -t 10 -o all.log -e all.err -p serial_requeue --job-name=${jobid}.reftrain --wrap=\"./reference_training.sh\" | awk ' { print \$4 }'`;
chomp $reftrainjobid;
print STDERR "Submitted job $reftrainjobid - to train Glimmer and Prodigal on reference\n";


# run JobArray scripts to assemble and predict putative proteins
# only submit with assembly memory requirements if fastqs are present in the input list
my $assembly = 0;
my $jobarrayid="";
if (grep(/fastq/,@ARGV)) {
	write_suffix_array_slurm_script("RunJobArrays.sh", "8000", "1", "60", "serial_requeue", "JobArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
	$jobarrayid=`sbatch -d afterok:$reftrainjobid --array=1-$count --job-name=$jobid.glim --wrap=\"./RunJobArrays.sh" | awk ' { print \$4 }'`;
	chomp $jobarrayid;
	print STDERR "Submitted job $jobarrayid - to extract putative proteins\n";
	$assembly = 1;
} else {
	write_suffix_array_slurm_script("RunJobArrays.sh", "2000", "1", "60", "serial_requeue", "JobArray"); # (name of slurm job array batch script (must match in sbatch below), memmory, nodes, time in minutes, queue, prefix of indexed jobs for job array)
	$jobarrayid=`sbatch -d afterok:$reftrainjobid --array=1-$count --job-name=$jobid.glim --wrap=\"./RunJobArrays.sh" | awk ' { print \$4 }'`;
	chomp $jobarrayid;
	print STDERR "Submitted job $jobarrayid - to extract putative proteins\n";
}


# run checkpointing following gene prediction
my $clustercheckpointid=`sbatch -e clustercheck.err -o clustercheck.out -d afterok:$jobarrayid --mem=2000 -n 1 -t 10 --job-name=$jobid\".\"CHK -p serial_requeue -o all.log -e all.err --wrap=\"$script_dir/clustering_checkpoint.pl $count $jobid $assembly $refnum $reference\"| awk ' { print \$4 }'`;
chomp $clustercheckpointid;
print STDERR "Submitted job $clustercheckpointid - to check clusters\n";
#system "bsub -J \"$jobid\".CHK -w \"ended($jobid"."glim[1-$count])\" -o all.log -e all.err $script_dir/clustering_checkpoint.pl $count $jobid $assembly $refnum $reference";                   to see exactly what this module does.	


# print files for COGtriangles clustering (three stages have different memory requirements)
open COGA, "> cogtriangle_run_A.sh" or croak();
open COGB, "> cogtriangle_run_B.sh" or croak();
open COGC, "> cogtriangle_run_C.sh" or croak();

# part A - low memory, processes BLAT results
print COGA "rm *.aa *.CDS.mfa;
mkdir blan;
mkdir blaf;
cat unfilteredblast.*.out > blan/raw.all.strains.blast.tab;
rm unfilteredblast.*.out;
cat filteredblast.*.out > blaf/filtered.all.strains.blast.tab;
rm filteredblast.*.out;
mkdir BLASTconv;
COGmakehash -i=all.strains.csv -o=./BLASTconv -s=\",\" -n=1;
COGreadblast -r -d=./BLASTconv/ -e=0.1 -q=1 -t=2;
rm *mod;
mkdir tmp;
COGlse -d=./BLASTconv/ -j=LSE.csv -p=all.strains.csv -o=all.strains.lse.csv -t=./tmp;
echo 'Completed run_cogtriangles_A.sh!' >> all.err";
close COGA;

# part B - high memory, runs clustering
print COGB "
COGtriangles -i=./BLASTconv -q=all.strains.csv -l=all.strains.lse.csv -o=all.strains.cls.csv -t=0.5 -e=0.01 -n=\"CLS\" -s=1;
echo 'Completed run_cogtriangles_B.sh!' >> all.err";
close COGB;

# part C - low memory processes the clustering output into a unique cluster for each CDS
print COGC "$script_dir/make_COG_addendum.pl all.strains.csv all.strains.cls.csv;
cat all.strains.cls.csv addendum.out > all.strains.cls.chk.csv;
rm all.strains.cls.csv;
COGcognitor -i=./BLASTconv -t=all.strains.cls.chk.csv -q=all.strains.csv -o=all.strains.cls.out.csv -c=1;
rm all.strains.cls.chk.csv;
echo 'Completed run_cogtriangles_C.sh!' >> all.err";
close COGC;

system "chmod +x cogtriangle_run_A.sh";
system "chmod +x cogtriangle_run_B.sh";
system "chmod +x cogtriangle_run_C.sh";


##############
# SUBROUTINES #
##############

sub croak {
	die print STDERR "
Error: $!\nUsage: -r [reference fasta] [list of fastqs/fastas/multifastas]

The reference genome will not feature in the final analysis unless it also features in the list

";
}

sub write_suffix_array_slurm_script {
	my $scriptname=$_[0];
	my $memory=$_[1];
	my $nodes=$_[2];
	my $minutes=$_[3];
	my $queue=$_[4];
	my $stageprefix=$_[5];
	open SLURMSCRIPTFH, "> $scriptname";
	print SLURMSCRIPTFH "#!/bin/bash\n";
  	print SLURMSCRIPTFH "#SBATCH --mem=$memory\n";
  	print SLURMSCRIPTFH "#SBATCH -n $nodes\n";
  	print SLURMSCRIPTFH "#SBATCH -t $minutes\n";
	print SLURMSCRIPTFH "#SBATCH -p $queue\n";
	print SLURMSCRIPTFH "#SBATCH -o log_%A_%a\n";
	print SLURMSCRIPTFH "#SBATCH -e err_%A_%a\n";
	print SLURMSCRIPTFH "./$stageprefix.";
	print SLURMSCRIPTFH '$SLURM_ARRAY_TASK_ID';
	close SLURMSCRIPTFH;
	system "chmod +x $scriptname";
}
