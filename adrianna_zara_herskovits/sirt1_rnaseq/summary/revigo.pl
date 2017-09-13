use strict;
use WWW::Mechanize;
use WWW::Mechanize::FormFiller;
use URI::URL;
use Getopt::Long;


#import data from R
my $goterms;
my $gopvals;
my $cutoff;
my $organism;
my $ispvalue;
my $whatisbetter;
my $measure;
my $domain;

GetOptions ("goterms=s" => \$goterms,
            "gopvals=s" => \$gopvals, 
            "cutoff=s"    => \$cutoff,
			"organism=s" => \$organism,
			"ispvalue=s" => \$ispvalue,
			"whatisbetter=s" => \$whatisbetter,
			"measure=s" => \$measure,
         "domain=s"=> \$domain);
			 
# prep go terms and pvals to be put into hash and joined
## parse options strings into arrays
my @go_terms=split(/,/, $goterms);
my @go_pvals=split(/,/, $gopvals);
## put arrays into hash as keys and values
my %go_hash;
@go_hash{@go_terms} = @go_pvals;
## join keys and values of hash into long string
my $go_str = join("\n", map { "$_  $go_hash{$_}" } keys %go_hash);
$go_str = $go_str."\n";


# setup mechanize
my $agent = WWW::Mechanize->new( autocheck => 1 );
my $formfiller = WWW::Mechanize::FormFiller->new();
$agent->env_proxy();

# go to page and fill out forms with options
$agent->get('http://revigo.irb.hr/');
$agent->form_number(1) if $agent->forms and scalar @{$agent->forms};
$formfiller->add_filler( 'goList' => Fixed => $go_str);
$formfiller->add_filler( 'cutoff' => Fixed => $cutoff );
$formfiller->add_filler( 'isPValue' => Fixed => $ispvalue );
$formfiller->add_filler( 'whatIsBetter' => Fixed => $whatisbetter );
$formfiller->add_filler( 'goSizes' => Fixed => $organism );
$formfiller->add_filler( 'measure' => Fixed => $measure );
$formfiller->fill_form($agent->current_form);

# run analysis
$agent->click("startRevigo");

if($domain eq "MF"){
 $agent->follow_link( url => "toR_treemap.jsp?table=3");
} elsif($domain eq "CC"){
    $agent->follow_link( url => "toR_treemap.jsp?table=2");
    } elsif($domain eq "BP"){
          $agent->follow_link( url => "toR_treemap.jsp?table=1");
    }


# remove header from results
$agent->delete_header;

# dump results
my $results=$agent->content;

$results =~ s/tmPlot/treemap/g;
$results =~ s/;//g;
$results =~ s/CCCCCCAA/CCCCCC00/;
$results =~ s/position.legend = "none"/position.legend = "none",/;

my $addstring = "border.col=c(\"black\",\"gray35\"),\npalette.HCL.options=list(hue_start=80,hue_perm=TRUE),\nborder.lwds=c(5,2),\nfontsize.labels=c(36,18),\nfontcolor.labels=c(\"black\",\"gray34\"),\nfontsize.title=16";
my @add = split(/\n/ ,$addstring);


my @splitresults = split("\n", $results);

my @finalscript;

foreach (@splitresults) {
   if ($_ =~ /^pdf/){
      next;
     } elsif ($_ =~ /^dev.off/) {
      next;
      } elsif ($_ =~ /position.legend/)  {
      push @finalscript, $_;
      push @finalscript, @add;
   } else {
      push @finalscript, $_;
   }
}


print join("\n", @finalscript);
