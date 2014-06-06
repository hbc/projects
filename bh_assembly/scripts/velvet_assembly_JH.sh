#!/bin/bash

usage()
{
cat << EOF

USAGE: $0 [options] 

Optimises parameters in Velvet.

OPTIONS:
   -e [int]		Expected coverage			[Default = 20]
   -d			Dirty mode. Does not remove Graph files etc. from velvet folder
   -f [fastq file]	Fastq file containing (forward/shuffled) reads
   -i [int]		Mean insert length 			[Default = 180]
   -k [int]		Starting kmer value 			[Odd number from 21 to 3/4 of read length. Default = 2/3 of read length]
   -n			Turn off scaffolding (no Ns in assembly)
   -p			Paired-end
   -w			WGA data
   -r [fastq file]	Fastq file containing reverse reads	[Optional]
   -s [file name]	Name for shuffled fastq output file 	[Optional]
   -o ["string"]	Any other velvetg commands
   -h			Show this message

EXAMPLES:
   Single end reads starting from a kmer of 23:
   $0 -f bob.fastq -k 23
   Paired end reads in two fastq files with scaffolding turned off:
   $0 -p -n -i 180 -f bob_1.fastq -r bob_2.fastq -s bob_shuffled.fastq 
   Paired end reads in a single shuffled fastq file with minimum contig length set to 100 and read tracking on:
   $0 -p -i 180 -f bob_shuffled.fastq -o "-min_contig_lgth 100 -read_trkg yes"

   
EOF
}


script_dir="/net/hsphfs1/srv/export/hsphfs1/share_root/chb/projects/bh_assembly/scripts/"
velvet_dir=""

PAIRED=0
EXPCOV=20
INSERT=400
KMER=0
SCAFFOLD=yes
FORWARD=
REVERSE=
SHUFFLED=
OTHERS=
DIRTY=0
WGA=0
while getopts “hi:e:f:k:npwr:s:o:d” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             INSERT=$OPTARG
             ;;
         d)
             DIRTY=1
             ;;
         e)
             EXPCOV=$OPTARG
             ;;
         f)
             FORWARD=$OPTARG
             ;;
         k)
             KMER=$OPTARG
             ;;
         p)
             PAIRED=1
             ;;
         w)
             WGA=1
             ;;
         n)
             SCAFFOLD=no
             ;;
         r)
             REVERSE=$OPTARG
             ;;
         s)
             SHUFFLED=$OPTARG
             ;;
         o)
             OTHERS=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


if [ "$WGA" -eq "1" ]
	then velvet_dir="/nfs/users/nfs_s/sh16/velvet-sc/"
fi




if [ -z "$FORWARD" ]
then
	echo
	echo ERROR: No fastq file selected
	echo
    usage
    exit 1
fi
if [ -n "$REVERSE" ] && [ -z "$SHUFFLED" ]
then
	echo
	echo ERROR: You need to give a file name for the shuffled sequences
	echo
	usage
	exit 1
fi

if [ "${FORWARD##*.}" == "gz" ]
    #then READLENGTH=$(zcat $FORWARD | head -n 2 | tail -n 1 |wc -L)
	then NFORWARD=${FORWARD%.gz}
	gunzip -c $FORWARD > $NFORWARD
	FORWARD=$NFORWARD
fi

READLENGTH=$(head -n 2 $FORWARD | tail -n 1 |wc -L)

if [ "${REVERSE##*.}" == "gz" ]
    #then READLENGTH=$(zcat $FORWARD | head -n 2 | tail -n 1 |wc -L)
	then NREVERSE=${REVERSE%.gz}
	gunzip -c $REVERSE > $NREVERSE
	REVERSE=$NREVERSE
fi

if [ $READLENGTH -lt "54" ]
	then READLENGTH=54
fi

MAX_KMER=$[$READLENGTH*3/4]

if [ $MAX_KMER -gt "127" ]
	then MAX_KMER=127
fi

MIN_KMER=21


if [ $KMER -eq "0" ]
	then KMER=$[$READLENGTH*2/3]
fi

if [ $KMER -gt "127" ]
	then KMER=127
fi

rem=$(( $KMER % 2 ))
 
if [ $rem -eq 0 ]
	then KMER=$[$KMER+1]
fi



lastcov=$EXPCOV
mincov=$[$EXPCOV/2] #change this to change the minimum coverage that can be accepted
covcutoff=$mincov
if [ "$mincov" -lt "10" ]
	then
	if [ "$EXPCOV" -ge "10" ]
		then
		mincov=10
		fi
	fi

readtype="short"





if [ "$PAIRED" -eq "1" ]
	then readtype="shortPaired"
fi
if [ -n "$REVERSE" ]
	then echo Shuffling sequences...
	shuffleSequences_fastq.pl $FORWARD $REVERSE $SHUFFLED
	fastq=$SHUFFLED
else fastq=$FORWARD
fi

logfile=${fastq%.fastq}_velvet_assembly.log

echo velvet_assembly.sh $@ > $logfile

echo Running velveth with kmer of $KMER
if [ "${fastq##*.}" == "gz" ]
    then ${velvet_dir}velveth ${fastq%.fastq}_velvet $KMER -fastq.gz -$readtype $fastq > ${fastq%.fastq}_velvet.log
else ${velvet_dir}velveth ${fastq%.fastq}_velvet $KMER -fastq -$readtype $fastq > ${fastq%.fastq}_velvet.log
fi
echo Running velvetg

if [ "$PAIRED" -eq "1" ]
	then ${velvet_dir}velvetg ${fastq%.fastq}_velvet -ins_length $INSERT > ${fastq%.fastq}_velvet.log

else ${velvet_dir}velvetg ${fastq%.fastq}_velvet > ${fastq%.fastq}_velvet.log

fi

expcov=$(${script_dir}velvet_stats_2_av_cov.py ${fastq%.fastq}_velvet/stats.txt $mincov)

echo Found expected coverage of $expcov
echo kmer=$KMER expected coverage=$expcov  >> $logfile
grep "n50" ${fastq%.fastq}_velvet.log | head -n 1
grep "n50" ${fastq%.fastq}_velvet.log | head -n 1 >> $logfile

if [ "$expcov" -gt $[$EXPCOV] ]

	then KMER=$[$KMER+2]
	while [ "$KMER" -lt $MAX_KMER ]
	do
		echo Expected coverage above $EXPCOV, increasing kmer
		echo Running ${velvet_dir}velveth with kmer of $KMER
		
		if [ "${fastq##*.}" == "gz" ]
		then ${velvet_dir}velveth ${fastq%.fastq}_velvet $KMER -fastq.gz -$readtype $fastq > ${fastq%.fastq}_velvet.log
		else ${velvet_dir}velveth ${fastq%.fastq}_velvet $KMER -fastq -$readtype $fastq > ${fastq%.fastq}_velvet.log
		    fi
		
		echo Running ${velvet_dir}velvetg
		if [ "$PAIRED" -eq "1" ]
			then ${velvet_dir}velvetg ${fastq%.fastq}_velvet -ins_length $INSERT > ${fastq%.fastq}_velvet.log
		else ${velvet_dir}velvetg ${fastq%.fastq}_velvet > ${fastq%.fastq}_velvet.log
		fi
		expcov=$(${script_dir}velvet_stats_2_av_cov.py ${fastq%.fastq}_velvet/stats.txt $mincov)
		echo Found expected coverage of $expcov
		echo kmer=$KMER expected coverage=$expcov  >> $logfile
		grep "n50" ${fastq%.fastq}_velvet.log | head -n 1
		grep "n50" ${fastq%.fastq}_velvet.log | head -n 1 >> $logfile
		if  [ "$expcov" -lt $[$EXPCOV] ]
			then echo Found best combination:
			if  [ "$KMER" -gt $[$MIN_KMER] ]
				then KMER=$[$KMER-2]
				expcov=$lastcov
			fi
			echo kmer=$KMER expected coverage=$expcov
			break
		fi
		lastcov=$expcov
		KMER=$[$KMER+2]
	done

elif [ "$expcov"  -lt $[$EXPCOV] ]
	then KMER=$[$KMER-2]
	while [ "$KMER" -ge $MIN_KMER ]
	do
		echo Expected coverage below $EXPCOV, decreasing kmer
		echo Running velveth with kmer of $KMER
		
		if [ "${fastq##*.}" == "gz" ]
		then ${velvet_dir}velveth ${fastq%.fastq}_velvet $KMER -fastq.gz -$readtype $fastq > ${fastq%.fastq}_velvet.log
		else ${velvet_dir}velveth ${fastq%.fastq}_velvet $KMER -fastq -$readtype $fastq > ${fastq%.fastq}_velvet.log
		    fi
		echo Running velvetg
		if [ "$PAIRED" -eq "1" ]
			then ${velvet_dir}velvetg ${fastq%.fastq}_velvet -ins_length $INSERT > ${fastq%.fastq}_velvet.log
		else ${velvet_dir}velvetg ${fastq%.fastq}_velvet > ${fastq%.fastq}_velvet.log
		fi
		expcov=$(${script_dir}velvet_stats_2_av_cov.py ${fastq%.fastq}_velvet/stats.txt $mincov)
		echo Found expected coverage of $expcov
		echo kmer=$KMER expected coverage=$expcov  >> $logfile
		grep "n50" ${fastq%.fastq}_velvet.log | head -n 1
		grep "n50" ${fastq%.fastq}_velvet.log | head -n 1 >> $logfile
		if  [ "$expcov" -gt $[$EXPCOV-1] ]
			then echo Found best combination:
			echo kmer=$KMER expected coverage=$expcov
			break
		fi
		KMER=$[$KMER-2]
	done
else
	echo Found best combination:
	echo kmer=$KMER expected coverage=$expcov
fi

if [ "$KMER" -lt $MIN_KMER ]
	then KMER=$MIN_KMER
fi
if [ "$KMER" -gt $MAX_KMER ]
        then KMER=$MAX_KMER
fi

covcutoff=$[$expcov/3]
if [ "$covcutoff" -gt "10" ]
		then covcutoff=10
fi

echo Final run parameters: >> $logfile
echo kmer=$KMER expected coverage=$expcov coverage cutoff=$covcutoff  >> $logfile 

echo Running final run with kmer of $KMER, expected coverage of $expcov and coverage cutoff of $covcutoff

if [ "$PAIRED" -eq "1" ]
	echo ${velvet_dir}velvetg ${fastq%.fastq}_velvet -ins_length $INSERT -scaffolding $SCAFFOLD -exp_cov $expcov -cov_cutoff $covcutoff $OTHERS
	then ${velvet_dir}velvetg ${fastq%.fastq}_velvet -ins_length $INSERT -scaffolding $SCAFFOLD -exp_cov $expcov -cov_cutoff $covcutoff $OTHERS > ${fastq%.fastq}_velvet.log
else echo ${velvet_dir}velvetg ${fastq%.fastq}_velvet -exp_cov $expcov -cov_cutoff $covcutoff $OTHERS
${velvet_dir}velvetg ${fastq%.fastq}_velvet -exp_cov $expcov -cov_cutoff $covcutoff $OTHERS > ${fastq%.fastq}_velvet.log
fi

grep "n50" ${fastq%.fastq}_velvet.log | head -n 1
grep "n50" ${fastq%.fastq}_velvet.log | head -n 1 >> $logfile

rm ${fastq%.fastq}_velvet.log
if [ "$DIRTY" -eq "0" ]
	then rm -rf ${fastq%.fastq}_velvet/*Graph* ${fastq%.fastq}_velvet/Sequences ${fastq%.fastq}_velvet/Roadmaps ${fastq%.fastq}_velvet/Log
fi

newname=${fastq%.fastq}

sed 's/length_[1234567890]*_//g' ${fastq%.fastq}_velvet/contigs.fa | sed  s/NODE/"${newname##*/}"/g | sed 's/\.[1234567890]*//g' > ${fastq%.fastq}_velvet/contigs.renamed
mv ${fastq%.fastq}_velvet/contigs.renamed ${fastq%.fastq}_velvet/contigs.fa

echo "Assembly finished"
