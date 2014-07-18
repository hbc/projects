#! /usr/bin/perl
use strict;
use warnings;

my $get;

my $first = 1;

while (<>)
{

    if (/^>(\S+)/){
	
	if ($first){
	    $first = 0;
	}
	else{
		print "//\n";
	}
		
	print "ID   $1             \n";
    }
    #if (/[\s+]?(\d+)\s+(\d+)\s+(\d+)\s+\[([+-]\d+)\s+L=\s*\d+\]\s*(.*)/)
    elsif (/^(\w+)\s+(\d+)\s+(\d+)\s+([+-])\d\s+(\d+\.*\d*)/)
    {
        #my ($id, $start, $end, $frame, $message) = ($1, $2, $3, $4, $5);
        my ($id, $start, $end, $direction, $score) = ($1, $2, $3, $4, $5);
	my $location;
	

        if ($direction eq '+' )
        {
            $location = "$start..$end";
        }
        else
        {
            $location = "complement($end..$start)";
        }

	unless ($direction  eq '+' && $start > $end || $direction  eq '-' && $start < $end){

       	 print "FT   CDS             $location\n";
       	 print "FT                   /note=\"Raw score $score\"\n";
       	 print "FT                   /label=$id\n";
       	 print "FT                   /colour=4\n";
       	 print "FT                   /method=\"GLIMMER\"\n";
	}
    }
}

