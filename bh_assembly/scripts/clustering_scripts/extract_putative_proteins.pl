#! /usr/bin/perl -w
use strict;
use Bio::SeqIO;

###################################################################################
# Submitted by run_clustering.pl | Submits nothing                                #
# Runs gene predictions, writes out EMBL files and protein sequences              #
###################################################################################

# set up script_dir variable, slurm variables
my $script_dir = $ENV{SCRIPT_DIR};

# the scripts called by this script are on directory higher in the tree than the clustering scripts, chop off last directory
$script_dir =~ s/\/[^\/]+$//;
my $slurmqueue = $ENV{SLURMQUEUE};
my $slurmtime  = $ENV{SLURMTIME};
my $slurmmem   = $ENV{SLURMMEM};

# set up other variables

my $tab;
my $embl;
my $emblseq;
my $header;
my $embltemp;
my $dna;
my $count = $ARGV[1];
my %tracking;
my $file      = $ARGV[0];    # original fastq file
my $mfa       = $file;
my $reference = $ARGV[2];    # original reference fa file
my $mod;

########################################################################
# FILE SETUP

# set up file names for each strain
# unzipping files and finding
if ( $file =~ /_1.fastq.gz/ ) {
    my $unzipped = $file;
    $unzipped =~ s/.gz//g;
    system "gunzip -c $file > $unzipped";
    $file = $unzipped;
    if ( $file =~ /_1.fastq/ ) {
        my $reverse = $file;
        $reverse =~ s/_1.fastq/_2.fastq/g;
        if ( -e "$reverse.gz" ) {
            system "gunzip -c $reverse.gz > $reverse";
        }
        else {
            print STDERR "Cannot find reverse fastq to unzip!\n";
            exit(1);
        }
    }
    $mfa =~ s/.gz//g;
    $mfa =~ s/_1.fastq/.fa/g;
}

$dna = $embl = $emblseq = $tab = $header = $embltemp = $mod = $file;

unless ( $file =~ /.dna/ ) {
    $dna =~ s/_1.fastq|.fasta|.seq|.fa/.dna/g;
}

if ( $file =~ /.dna|.fasta|.seq/ ) {
    $mfa =~ s/.dna|.fasta|.seq/.fa/g;
    system "cp $file $mfa";
}

$tab =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.tab/g;
$emblseq =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.emblseq/g;
$embl =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.embl/g;
$header =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.header/g;
$embltemp =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.temp/g;
$mod =~ s/_1.fastq|.dna|.fasta|.seq|.fa/.mod/g;

########################################################################
# VELVET and ABACAS
# shuffle data if provided as fastq and run abacas
if ( $file =~ /_1.fastq/ ) {
    my $working_dir = `pwd`;
    chomp $working_dir;
    my $forward = $file;
    my $reverse = $file;
    $reverse =~ s/_1.fastq/_2.fastq/g;
    system
"$script_dir/velvet_assembly.sh -f $forward -r $reverse -p -s shuffled.$file";
    $file =~ s/.fastq//g;
    system "cp $reference shuffled.$file" . "_velvet";
    chdir "$working_dir/shuffled.$file" . "_velvet";
    system
"$script_dir/abacas.1.3.1.pl -abcNm -p nucmer -q contigs.fa -r $reference";
    system "cat contigs.fa_"
      . "$reference.MULTIFASTA.fa contigs.fa_"
      . "$reference.contigsInbin.fas > $mfa";
    system "cat $mfa | grep -v '>' > joined.seq";
    $file =~ s/_1/.fa/g;
    system "echo \">$file\" | cat - joined.seq > ../$file";
    chdir "$working_dir";
}

########################################################################

# add the three frame stop to the ends of the assembly

open FILE, $mfa
  or die print STDERR
  "Could not open file $mfa corresponding to input $file!\n";

open MFA, "> $mod" or croak();
open DNA, "> $dna" or croak();

my $file_header = 0;

print DNA ">$file\n";

my ( $bit, $end );

foreach (<FILE>) {
    chomp;
    if (/^>/) {
        print MFA "$_\n";
    }
    elsif ( $file_header == 0 ) {
        print MFA "TTATTTATTTA" . "$_\n";
        $file_header = 1;
        $bit         = substr( $_, 0, 49 );
        $end         = substr( $_, 49, 11 );
        print DNA "TTATTTATTTA" . "$bit";
    }
    else {
        print MFA "$_\n";
        $bit = substr( $_, 0, 49 );
        print DNA "\n$end" . "$bit";
        if ( length($_) > 49 ) {
            $end = substr( $_, 49, 11 );
        }
        else {
            $end = "";
        }
    }
}

print MFA "TAAATAAATAA\n";

my $line = $bit . $end . "TAAATAAATAA";

if ( length($line) > 60 ) {
    my $p        = $bit . $end;
    my $prev     = length($p);
    my $residual = 60 - $prev;
    $bit = substr( $line, $prev, $residual );
    $end = substr( $line, 60 );
    print DNA "$bit\n$end\n";
}
else {
    print DNA "$line\n";
}

close FILE;
close MFA;

# identify the positions of contig breaks

my @breaks;
my $pos_count = 1;

open IN, $mod or die print STDERR "Could not open file $file; died\n";

foreach (<IN>) {
    chomp;
    if (/>/) {
        push( @breaks, $pos_count );
    }
    else {
        my @data = split( //, $_ );
        $pos_count += scalar(@data);
    }
}

close IN;

# predict CDSs using Glimmer and Prodigal and create EMBL

unless ( $file =~ /.dna/ ) {

#	system "grep -v '>' $mod | seqret -filter | sed 's/EMBOSS_001/$file/g' > $dna";
}

print STDERR "Predicting genes for $file using Prodigal...";
if ( -e "all.strains.prod.train" ) {
    system "prodigal -t all.strains.prod.train < $dna > $file.prod.tab";
}
else {
    die print STDERR
      "Could not use Prodigal training file all.strains.prod.train\n";
}

open PROD, "$file.prod.tab"
  or die print STDERR "Prodigal gene prediction on $file failed!\n";

my $p = 0;
my %prodigal_start;
my %prodigal_end;
my %prodigal_strand;

foreach (<PROD>) {
    chomp;
    $_ =~ s/\.|\(|\)/ /g;
    $_ =~ s/>|<//g;
    my @data = split( /\s+/, $_ );
    if (/^     CDS/) {
        if (/complement/) {
            $prodigal_start{$p}  = $data[3];
            $prodigal_end{$p}    = $data[4];
            $prodigal_strand{$p} = -1;
        }
        else {
            $prodigal_start{$p}  = $data[2];
            $prodigal_end{$p}    = $data[3];
            $prodigal_strand{$p} = 1;
        }
    }
    $p++;
}

close PROD;

#system "rm $file.prod.tab";

print STDERR "Predicting genes for $file using Glimmer3...";
system "glimmer3 $dna all.strains.icm $file.glim";

# remove predicted CDSs that span contig breaks

open PRED, "$file.glim.predict"
  or die print STDERR "Glimmer3 gene prediction on $file failed!\n";
open RED, "> filtered.$file.glim.predict"
  or die print STDERR
  "Could not write data to file filtered.$file.glim.predict\n";

my $broken         = 0;
my $prodigal_match = 0;
my $break;
my $comp_start;
my $comp_end;

foreach (<PRED>) {
    chomp;
    my @data = split( /\s+/, $_ );
    $broken         = 0;
    $prodigal_match = 0;
    if (/^>/) {
        print RED "$_\n";
    }
    else {
        if ( $data[3] > 0 ) {
            my $gob1 = $data[1] + 0.25 * ( $data[2] - $data[1] );
            my $gob2 = $data[2] - 0.25 * ( $data[2] - $data[1] );
            foreach $p ( keys %prodigal_start ) {
                if ( $prodigal_strand{$p} == 1 ) {
                    my $pob1 =
                      $prodigal_start{$p} +
                      0.25 * ( $prodigal_end{$p} - $prodigal_start{$p} )
                      ; # look for overlap in the central 50%, ignore overlaps that only involve the ends
                    my $pob2 = $prodigal_end{$p} -
                      0.25 * ( $prodigal_end{$p} - $prodigal_start{$p} );
                    if (   ( $pob1 >= $gob1 && $pob1 <= $gob2 )
                        || ( $pob2 >= $gob1 && $pob2 <= $gob2 ) )
                    {
                        if ( ( $prodigal_end{$p} - $prodigal_start{$p} ) <
                            ( $data[2] - $data[1] ) )
                        {
                            $comp_start = $data[1];
                            $comp_end   = $data[2];
                        }
                        else {
                            $comp_start = $prodigal_start{$p};
                            $comp_end   = $prodigal_end{$p};
                        }
                        $prodigal_match = 1;
                        foreach $break (@breaks) {
                            if ( $break >= $comp_start && $break <= $comp_end )
                            {
                                if ( $broken == 0 || $break < $broken ) {
                                    $broken = $break;
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            my $gob1 = $data[2] + 0.25 * ( $data[1] - $data[2] );
            my $gob2 = $data[1] - 0.25 * ( $data[1] - $data[2] );
            foreach $p ( keys %prodigal_start ) {
                if ( $prodigal_strand{$p} == -1 ) {
                    my $pob1 =
                      $prodigal_start{$p} +
                      0.25 * ( $prodigal_end{$p} - $prodigal_start{$p} )
                      ; # look for overlap in the central 50%, ignore overlaps that only involve the ends
                    my $pob2 = $prodigal_end{$p} -
                      0.25 * ( $prodigal_end{$p} - $prodigal_start{$p} );
                    if (   ( $pob1 >= $gob1 && $pob1 <= $gob2 )
                        || ( $pob2 >= $gob1 && $pob2 <= $gob2 ) )
                    {
                        if ( ( $prodigal_end{$p} - $prodigal_start{$p} ) <
                            ( $data[2] - $data[1] ) )
                        {
                            $comp_start = $data[1];
                            $comp_end   = $data[2];
                        }
                        else {
                            $comp_start = $prodigal_start{$p};
                            $comp_end   = $prodigal_end{$p};
                        }
                        $prodigal_match = 1;
                        foreach $break (@breaks) {
                            if ( $break >= $comp_start && $break <= $comp_end )
                            {
                                if ( $broken == 0 || $break > $broken ) {
                                    $broken = $break;
                                }
                            }
                        }
                    }
                }
            }
        }
        if ( $prodigal_match == 1 ) {
            if ( $broken == 0 ) {
                if ( $data[3] > 0 ) {
                    print RED
                      "$data[0]\t$comp_start\t$comp_end\t$data[3]\t$data[4]\n";
                }
                else {
                    print RED
                      "$data[0]\t$comp_end\t$comp_start\t$data[3]\t$data[4]\n";
                }
            }
            else {
                if ( $data[3] > 0 ) {
                    my $new_end = $broken - 1;
                    if ( abs( 1 + $new_end - $comp_start ) >= 99 ) {
                        print RED
"$data[0]\t$comp_start\t$new_end\t$data[3]\t$data[4]\n";
                    }
                }
                else {
                    my $new_end = $broken;
                    if ( abs( 1 + $new_end - $comp_end ) >= 99 ) {
                        print RED
                          "$data[0]\t$comp_end\t$new_end\t$data[3]\t$data[4]\n";
                    }
                }
            }
        }
    }
}

close PRED;
close RED;

system "$script_dir/glimmer3totab.pl filtered.$file.glim.predict > $tab";
system "cat $dna | seqret -filter -osformat EMBL > $embltemp";

open TMP,  "$embltemp";
open HEAD, "> $header";
open SEQ,  "> $emblseq";

foreach (<TMP>) {
    if (/^ID|^DE|^XX/) {
        print HEAD "$_";
    }
    else {
        print SEQ "$_";
    }
}

print HEAD
  "XX\nDE   Draft genome\nXX\nFH   Key             Location/Qualifiers\nFH\n";

close TMP;
close SEQ;
close HEAD;
system "grep -v 'ID' $tab | sed s/label=\"orf/label=\"seq" . "$count"
  . "orf/g | cat $header - $emblseq > $embl";

#system "rm $file.glim*";
#system "rm $tab $emblseq $header $embltemp filtered.$file.glim.predict";

# get assembly stats for each of the strains

my $stats     = `$script_dir/fac.pl $file | head -n 2 | tail -n 1 | tr -d "\n"`;
my @stats     = split( /\s+/, $stats );
my $cds_count = 0;

# extract protein sequences from EMBL and write CVS file for orthologue analysis

print STDERR "done\nExtracting protein sequences...";
open CSV,  "> seq.$count.csv";
open PSEQ, "> seq.$count.protseq";
open DSEQ, "> seq.$count.DNAseq";

my $infile     = Bio::SeqIO->new( -file => $embl, -format => "EMBL" );
my $annotation = $infile->next_seq;
my $outname    = $embl;
$outname =~ s/embl/aa/g;
my $DNA_outname = $embl;
$DNA_outname =~ s/embl/CDS.mfa/g;
my $outfile = Bio::SeqIO->new( -format => 'Fasta', -file => "> $outname" );
my $DNA_outfile =
  Bio::SeqIO->new( -format => 'Fasta', -file => "> $DNA_outname" );

for my $gene ( $annotation->get_SeqFeatures ) {
    if ( $gene->primary_tag eq "CDS" ) {
        $cds_count++;
        for my $name ( $gene->get_tag_values('label') ) {
            print CSV "$name,seq" . "$count\n";
            my $DNA_seq    = $gene->seq;
            my $DNA_string = $DNA_seq->seq();
            $DNA_string =~
              s/R|Y|W|S|M|K|H|B|D|V|X/N/g;    # deal with ambiguity codes
            my $prot_seq  = $gene->seq->translate();
            my $aa_string = $prot_seq->seq;
            $prot_seq->display_name($name);
            $prot_seq->desc( "seq" . "$count" );
            $outfile->write_seq($prot_seq);
            $DNA_seq->display_name($name);
            $DNA_seq->desc( "seq" . "$count" );
            $DNA_outfile->write_seq($DNA_seq);
            print PSEQ "$name|$aa_string\n";
            print DSEQ "$name|$DNA_string\n";
        }
    }
}
close CSV;
close DSEQ;
close PSEQ;

# report assembly statistics

open REP, "> $file.report"
  or die print STDERR "Cannot open output file $file.report\n";
print REP "$file\t$stats[8]\t$stats[0]\t$stats[6]\t$cds_count\n";
close REP;

# filter extracted proteins with seg

system
"segmasker -in $outname -outfmt fasta | perl -lane 'if (/>/) {print \$_;} else {\$_ =~ tr/a-z/X/; print \$_;}' > filtered.$outname";

open TRACK, "> $file.done"
  or die print STDERR "Cannot open tracking file $file.done\n";
print TRACK "$file";
close TRACK;

print STDERR "done\nCompleted extract_putative_proteins.pl\n";
