#!/bin/bash
intersectBed -a ref/Homo_sapiens.GRCh37.68.gtf -b $1 | cut -f9 | cut -f1 -d";" | cut -f3 -d" " | awk '{count[$1]++} END {for (j in count) print j, count[j]}' | tr -d \" > $1.counts
