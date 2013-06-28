## run in directory with teh data
BEDFILES="open_adenoma1_top1000.bed
open_adenoma2_top1000.bed
open_adenoma3_top1000.bed
open_adenoma4_top1000.bed
open_adenoma5_top1000.bed
open_carcinoma1_top1000.bed
open_carcinoma2_top1000.bed
open_carcinoma3_top1000.bed
open_carcinoma4_top1000.bed
open_carcinoma5_top1000.bed
"

DOMAIN=$1

subdomain1=`basename "$DOMAIN" ".genome.bed"`
subdomain2=`echo ${subdomain1#gold.}`
echo $subdomain2

for BEDFILEA in $BEDFILES
do
	subnameA1=`basename "$BEDFILEA" "_top1000.bed"`
	subnameA2=`echo ${subnameA1#open_}`
	echo $subnameA2
		
	for BEDFILEB in $BEDFILES
		do
		subnameB1=`basename "$BEDFILEB" "_top1000.bed"`
		subnameB2=`echo ${subnameB1#open_}`
	  	echo " " $subnameB2
  		IntervalStats -q $BEDFILEA -r $BEDFILEB -d $DOMAIN -o ${subnameA2}.${subnameB2}.${subdomain2}.out
  done
done
