#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

statsfile=sys.argv[1]
if len(sys.argv)>2:
        mincov=float(sys.argv[2])
else:
        mincov=10

lines=open(statsfile,'rU').readlines()

bins={}

for line in lines[1:]:
        words=line.strip().split()
        bin=round(float(words[5]),0)
        if bins.has_key(bin):
                bins[bin]=bins[bin]+int(words[1])
        else:
                bins[bin]=int(words[1])
                
maxcov=0
best_cov=0


for bin in bins.keys():
        if bins[bin]>maxcov and bin > mincov-1:#the > mincov  is to stop problems with lots of reads with coverage of 1 spoiling the results - should change for different exp coverages... maybe can be handled better
                maxcov=bins[bin]
                best_cov=bin
                
print int(best_cov)
