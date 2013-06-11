#!/usr/bin/perl

 
use strict; 
use warnings;

open (INPUT, $ARGV[0]) or die;
while (<INPUT>) {
	my @array = split /\t/;
	print "@".$array[0]."_00".$array[1]."_FC182:".$array[2].":".$array[3].":".$array[4].":".$array[5]."#".$array[6]."/".$array[7]."\n".$array[8]."\n"."+".$array[0]."_00".$array[1]."_FC182:".$array[2].":".$array[3].":".$array[4].":".$array[5]."#".$array[6]."/".$array[7]."\n".$array[9]."\n";    
}	
close INPUT;


