#!/bin/bash
intersectBed -a ref/human_refgene.bed -b $1 | cut -f4 | cut -f1,2 -d"_" | awk '{print $1$2}' | awk '{count[$1]++} END {for (j in count) print j, count[j]}' > $1.counts

