# Generation of potential proteins using Pacbio data from the Ocwieja paper

## Downloading Pacbio Consensus Reads

Script to download data from SRA website uses a text file with the FTP site for each sample as a line in the file. Use script to download the data for each line in the file (sample)

```bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    wget $line
done < "$1"
```

From the SRA files, extract the FASTA sequences of the consensus reads. More information about the code can be found at http://seqanswers.com/forums/archive/index.php/t-58456.html

```bash
for sra in *.sra; do ~/tools/sratoolkit.2.7.0-centos_linux64/bin/fastq-dump.2.7.0 --fasta $sra; done
```

Combine all fasta files together into a single file for alignment since we do not care about where the sequence came from.

```bash
cat "names of all fasta files"
```

## Filtering the Pacbio sequences for length and contamination with RT primer

Remove sequences less than 40 nucleotides and sequences matching the RT primer similar to the Ocwieja paper (Ocwieja protocol removed sequences with subreads less than 100nt - these should have been removed prior to generation of the consensus reads, so we do not need to perform this step)

```bash
fasta_formatter -i SRR528772.fasta -w 10000 | fastx_clipper -l 40

fastx_clipper -a CTCCACACTAACACTTGTCTCTCCG #RTPrime
```

## Align Pacbio reads

Align consensus reads against HIV strain 89.6 reference and human reference using STAR. STAR is the recommended aligner by PacBio as detailed [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study:-Optimizing-STAR-aligner-for-Iso-Seq-data). However, GMAP is the aligner used in the Pacbio pipeline.

During alignment, we mapped to human and HIV at the same time to identify human contamination for removal. To perform alignment, gene annotation files were not used since it was not important that we receive the best human alignment, we only needed to know if the reads align to human or not. HIV strain 89.6 doesn't have a well-annotated genome, and our goal was to identify non-annotated junctions, so we did not include the annotation file for human or HIV. 

We used parameters suggested for aligning Pacbio reads using STAR from the [Pacbio recommendations](https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study:-Optimizing-STAR-aligner-for-Iso-Seq-data). In addition, to the suggested parameters, the thresholds were relaxed for identifying splice junctions, including the number of unique reads per splice junction (a single read for any junction type) and minimum overhang (nine nucleotides for all junctions). By relaxing these parameters, it is understood that many junctions may be false positives.

```bash
## Generating index for combined human and HIV genomes

module load seq/STAR/2.5.2b

STARlong \
--runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /n/data1/cores/bcbio/legall_hiv_pacbio/STAR/index \
--genomeFastaFiles /groups/shared_databases/igenome/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa /n/data1/cores/bcbio/legall_hiv_pacbio/references/hiv.U39362.fa

## Performing alignment 

module load seq/STAR/2.5.2b

STARlong \
--runThreadN 6 \
--runMode alignReads \
--genomeDir /n/data1/cores/bcbio/legall_hiv_pacbio/STAR/index \
--readFilesIn /n/data1/cores/bcbio/legall_hiv_pacbio/QC/filtered_ccs.fasta \
--outSAMattributes NH HI NM MD \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--outSJfilterOverhangMin 9 9 9 9 \
--outSJfilterCountUniqueMin 1 1 1 1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000
```
Approximately, 57.45% of reads were identified as too short.

To better identify any novel junctions, STARlong was used to perform 2-pass mapping, which uses the list of junctions from the 1st pass as gene annotations for the 2-pass with the --sjdbFileChrStartEnd option.

```bash
module load seq/STAR/2.5.2b

STARlong \
--runThreadN 6 \
--runMode alignReads \
--genomeDir /n/data1/cores/bcbio/legall_hiv_pacbio/STAR/index \
--readFilesIn /n/data1/cores/bcbio/legall_hiv_pacbio/QC/filtered_ccs.fasta \
--outSAMattributes NH HI NM MD \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--outSJfilterOverhangMin 9 9 9 9 \
--outSJfilterCountUniqueMin 1 1 1 1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000 \
--sjdbFileChrStartEnd /n/data1/cores/bcbio/legall_hiv_pacbio/STAR/SJ.out.tab

# Total number of junctions (HIV or human)
wc -l SJ.out.tab: 681 

# Total number of HIV junctions
grep "U39362.2" SJ.out.tab | wc -l: 617
```

The total number of input reads to STAR was 512,181, with an average input read length of 789. The number of reads that were uniquely mapping was 230,631 (45.03%). The number of identified splice junctions had the following breakdown: 219425 GT/AG, 6807 GC/AG, 8882 AT/AC, and 7236 non-canonical. The mismatch rate per base was 2.53%. 673 reads mapped to multiple loci (0.13%), and 842 mapped to too many loci (0.16%) and were dropped. **The main concern with the alignment is the loss of over half of the reads identified as 'too short' (53.37%). To note, the STAR default requires mapped length to be > 2/3 of the total read length.**

## Remove human contamination from the Pacbio reads

Currently, the aligned reads have both HIV reads and human contamination. To remove human contamination by extracting reads mapping to HIV. Reads mapping only to HIV 89.6 genome, U39362.2, were extracted.

```bash

module load seq/samtools/1.3 

## Convert SAM with header to BAM
samtools view -bS Aligned.out.sam > human_hiv_aligned.bam

## Sort and index bam to use ‘view’ command
samtools sort human_hiv_aligned.bam -o human_hiv_aligned_sorted.bam
samtools index human_hiv_aligned_sorted.bam

## Extract only those reads aligning to HIV (U39362.2)
samtools view human_hiv_aligned_sorted.bam "U39362.2" -o hiv_aligned.bam

## Check extraction
samtools view human_hiv_aligned_sorted.bam | grep U39362.2 | wc -l: 230016
samtools view human_hiv_aligned_sorted.bam | wc -l: 233262
samtools view hiv_aligned.bam | wc -l: 230016
```

## Nucleotide sequence extraction

After removing human contamination, 230,016 reads remained. To generate potential open reading frames, the nucleotide sequences were extracted. A BED file was required to generate the FASTA sequences, so the coordinates were first converted to BED format prior to extraction.

```bash
# Extracting HIV sequences from BAM file

module load seq/BEDtools/2.26.0

bedtools bamtobed -i ../STAR/pass2/hiv_aligned.bam -split > hiv_aligned.bed

# The "-split" argument reports each portion of a “split” BAM (i.e., having an “N” CIGAR operation) alignment as distinct BED intervals. Therefore, the BED file does not output the intronic sequences, only the exons if the read is split.

bedtools getfasta -fi ../references/hiv.U39362.fa -bed hiv_aligned.bed -name -fo hiv_aligned_reads_split.fa

# Combine exons for each read

## Create a list of read header names to automate merging exons together

grep ">" hiv_aligned_reads_split.fa | cut -c 2- | sort -u > header_list

## Create a script, transcript_sequences_extraction.sh, to combine exons:

# Echo name of transcript, then merge sequences of exons together

>for name in $(cat header_list.txt)
>do
>echo ">$name" >> hiv_aligned_reads_merged.fa
>grep -A1 $name hiv_aligned_reads_split.fa | grep -v $name >> hiv_aligned_reads_merged.fa
>done

# Exons for each read are on separate lines, so need to merge the lines together for each read

fasta_formatter -i hiv_aligned_reads_merged.fa -w 0 > hiv_aligned_reads_final.fa

grep ">" hiv_aligned_reads_final.fa | wc -l: 230016

```

## Identification of ORFs and potential proteins

To identify the potential open reading frames (ORFs) using the getorf tool from the Emboss suite of tools, any ORFs were identified at any location in the read sequences using standard code and alternative initiation codons. The lowest minimum nucleotide size (30) was used, and ORFs were defined as a region that began with a START codon and ended with a STOP codon. We only found ORFs on the forward sequence, as no known transcripts are known to be encoded on the reverse strand for HIV. The identified ORFs were output as potential proteins.

```bash
# Get ORFs using standard code with alternative initiation codons, min nucleotide size of 30, with ORF defined as a region that begins with a START and ends with a STOP codon, and only finding ORFs in the forward sequence (not including reverse complement - using emboss suite getorf command available on Orchestra version 6.6.0

module load seq/emboss/6.6.0

getorf -sequence hiv_aligned_reads_merged.fa -outseq ./pacbio_potential_orfs_split.fa -table 1 -find 1 -reverse No

# Does not wrap lines of peptides, so need to merge the lines together for each read

fasta_formatter -i pacbio_potential_orfs_split.fa -w 0 > pacbio_potential_orfs_merged.fa

cp pacbio_potential_orfs_merged.fa pacbio_potential_orfs_merged.fa_copy
```

463,629 potential proteins were identified using the getorf program. To determine whether the potential proteins from the Pacbio analysis corresponded to the potential proteins identified from the transcript coordinates given in the Ocwieja paper, the Pacbio proteins were compared to those in the paper using BLAT with a filter for minimum identity equal to 100%. 

```bash
## BLAT the sequences against the ocwieja transcripts to determine known transcripts for the pacbio data 

blat ../ocwieja_analysis/unique_896_potential_proteins.fa pacbio_potential_orfs_merged.fa -prot -minIdentity=100 -out=pslx pacbio_orfs_in_ocwieja_data_100ident.pslx
```
# Statistical analysis in R

Analyzing the correspondence between the Pacbio analysis proteins and the Ocwieja paper proteins by importing the output of BLAT into R.

```r
library(ggplot2)

# Looking out BLAT results
pacbio_to_ocwieja <- read.table("pacbio_orfs_in_ocwieja_data_100ident.pslx", header = F, skip = 2, sep = "\t")

#List of all ocwieja proteins identified
all_ocwieja_proteins <- scan("../../ocwieja_analysis/all_ocwieja_protein_names.txt", what=character())

all_ocwieja_proteins <- all_ocwieja_proteins[grep(">", all_ocwieja_proteins)]

all_ocwieja_proteins <- sapply(strsplit(as.character(all_ocwieja_proteins), ">"), "[", 2)
```

Exploring the output of BLAT, which aligned only those pacbio proteins with 100% minimum identity to the ocwieja paper-derived proteins.

```r
summary_all_proteins <- summary(pacbio_to_ocwieja$V14)
capture.output(summary_all_proteins, file = "summary_all_proteins.txt")
levels(pacbio_to_ocwieja$V14)
```
108 theoretical proteins were identified by the predicting ORFs from the transcripts in the Ocwieja paper

```r
length(unique(pacbio_to_ocwieja$V10))
```
149,239 of the 463,629 HIV Pacbio proteins were returned from the BLAT analysis as aligning to the proteins from the Ocwieja paper

To determine which proteins identified by Ocwieja are or are not present in the Pacbio proteins identified by BLAT

```r
proteins_present <- all_ocwieja_proteins[which(all_ocwieja_proteins %in% levels(pacbio_to_ocwieja$V14))]
```
94 proteins of the 108 are present in the Pacbio output

```r
proteins_not_present <- all_ocwieja_proteins[which(!(all_ocwieja_proteins %in% levels(pacbio_to_ocwieja$V14)))]
```
14 proteins are not present in the Pacbio output

```r
## Proteins identified during the Ocwieja analysis that did not have Pacbio protein align had the following headings: 

## [1] "vif2__5"  "vif2__9"  "vif2__10" "vif2__12" "vif2__17"
## [6] "vif2__27" "vif2__30" "vif2__32" "vif2__34" "vif2__35"
## [11] "vif2__49" "vpr3__11" "tat5__11" "env4__11"
```

Identification of which Ocwieja proteins had full pacbio potential proteins that aligned identically to them.

```r
## Indicating full length read by specifying the Query start = 0 (.psl is zero-based) and the Query end = the length of the Query

fl_matching <- subset(pacbio_to_ocwieja, V12 == 0 & V13 == V11)
levels(fl_matching$V14)
```

69 proteins from the Ocwieja paper had full-length Pacbio proteins that aligned to them; however, some of the full-length Pacbio proteins aligned uniquely to proteins, while other full-length reads aligned to multiple proteins

```r
summary_nonunique_proteins <- summary(fl_matching$V14)
capture.output(summary_nonunique_proteins, file = "summary_full_nonunique_proteins.txt")
```

Identification of the number of Pacbio proteins that had full-length reads aligning uniquely to proteins in the Ocwieja paper

```r
uniq_matching <- which(!(fl_matching$V10 %in% fl_matching$V10[duplicated(fl_matching$V10)]))

full_length_uniq_matching <- fl_matching[uniq_matching, ]
```

49,797 of the 149,239 of the Pacbio potential proteins were full-length and aligned uniquely to a single Ocwieja protein

```r
length(which(levels(full_length_uniq_matching$V14) %in% full_length_uniq_matching$V14))
```

51 proteins from the Ocwieja paper have full-length Pacbio proteins uniquely aligning to them

```r
summary_unique_proteins <- summary(full_length_uniq_matching$V14)
capture.output(summary_unique_proteins, file = "summary_full_unique_proteins.txt")
```

Identification of the number of proteins that had partial-length reads aligning uniquely to proteins in the Ocwieja paper

```r
partial_uniq_matching <- which(!(pacbio_to_ocwieja$V10 %in% pacbio_to_ocwieja$V10[duplicated(pacbio_to_ocwieja$V10)]))

partial_length_uniq_matching <- pacbio_to_ocwieja[partial_uniq_matching, ]
```

34,642 of the 149,239 pacbio potential proteins were partial-length and aligned uniquely to a single Ocwieja protein

```r
length(which(levels(partial_length_uniq_matching$V14) %in% partial_length_uniq_matching$V14))
```

44 proteins have partial-length pacbio reads uniquely aligning to them

```r
summary_partial_unique_proteins <- summary(partial_length_uniq_matching$V14)
capture.output(summary_partial_unique_proteins, file = "summary_partial_unique_proteins.txt")
```

How many of the partial-length pacbio reads align to proteins not identified with full-length uniquely aligning reads?

```r
partial_uniq_proteins <- levels(partial_length_uniq_matching$V14)[levels(partial_length_uniq_matching$V14) %in% partial_length_uniq_matching$V14]

full_length_uniq_proteins <- levels(full_length_uniq_matching$V14)[levels(full_length_uniq_matching$V14) %in% full_length_uniq_matching$V14]

length(which(!(partial_uniq_proteins %in% full_length_uniq_proteins)))
```

5 proteins with partial-length pacbio reads uniquely aligning to them did not have full-length pacbio reads uniquely aligning to them. 

Therefore, a total of 51 + 5 = 56 proteins of the 108 (51.9%)identified in the Ocwieja paper were supported by full-length or partial-length Pacbio potential proteins that were uniquely aligning. 

The proteins identified from the Pacbio analysis included novel proteins identified in the Ocwieja paper.

```r
# Plotting

summary_pb_to_oc <- pacbio_to_ocwieja %>% 
        group_by(V14) %>%
        summarise(no_rows = length(V14))

ggplot(summary_fl_matching) +
        geom_histogram(aes(x = no_rows))

summary_fl_matching <- fl_matching %>% 
        group_by(V14) %>%
        summarise(no_rows = length(V14))

summary_uniq_fl_matching <- full_length_uniq_matching %>% 
        group_by(V14) %>%
        summarise(no_rows = length(V14))

summary_uniq_pl_matching <- partial_length_uniq_matching %>% 
        group_by(V14) %>%
        summarise(no_rows = length(V14))
```
# Liftover from HIV strain 89.6 to NL4-3

## Change FASTA files to 2bit to use kenttools

```bash
~/tools/kentUtils/bin/faToTwoBit ../../references/hiv.U39362.fa U39362_hiv_896.2bit
~/tools/kentUtils/bin/faToTwoBit ../../references/AF324493_hiv_nl43_ref_seq.fasta AF324493_hiv_nl43.2bit
~/tools/kentUtils/bin/twoBitInfo U39362_hiv_896.2bit U39362_hiv_896.chromInfo
~/tools/kentUtils/bin/twoBitInfo AF324493_hiv_nl43.2bit AF324493_hiv_nl43.chromInfo
```

## Use BLAT to create PSL file aligning 89.6 to NL4-3

```bash
module load seq/blat/35

blat chain/AF324493_hiv_nl43.2bit chain/U39362_hiv_896.2bit psl/89_to_nl.psl -tileSize=12 -noHead -minScore=100

~/tools/kentUtils/bin/axtChain -psl ../psl/89_to_nl.psl AF324493_hiv_nl43.2bit U39362_hiv_896.2bit  89_to_nl.chain -linearGap=loose

```

## Change coordinate system by creating a LFT file

```bash
~/tools/kentUtils/bin/liftUp -pslQ ../psl/89_to_nl.psl 89_to_nl.lft warn 89_to_nl_old.psl
```
## Chain together the coordinates from the LFT file to create a CHAIN file

```bash
~/tools/kentUtils/bin/axtChain -psl ../psl/89_to_nl.psl AF324493_hiv_nl43.2bit U39362_hiv_896.2bit  89_to_nl.chain -linearGap=loose
```

## Make alignment nets from chains

```bash
~/tools/kentUtils/bin/chainNet 89_to_nl.chain AF324493_chrom.sizes U39362_chrom.sizes ../net/89_to_nl.net /dev/null
```

## Create liftover chain file

```bash
 ~/tools/kentUtils/bin/netChainSubset ../net/89_to_nl.net 89_to_nl.chain ../89_to_nl_chain_net.chain 
 ```
 
 ## Liftover of 89.6 coordinates to NL4-3 using net chain files
 
 ```bash
 ~/tools/kentUtils/bin/liftOver ../../getORFs/hiv_aligned.bed ../89_to_nl_chain_net.chain hiv_converted_to_nl.bed ../unMapped/unmapped_89_to_nl
 ```
# Ocwieja analysis

# Extracting sequences using a GTF-like file of transcripts

vim ocwieja_transcripts_bed.txt ":%s/^M/\r/g"
remove spaces and "" 
remove any old .fai files

#Using HIV 89.6 strain - accession # U39362.2
bedtools getfasta -fi U39362.2_hiv_sequence.fasta -bed ocwieja_transcripts_bed.txt -name -fo ocwieja_transcript_sequences.fa

# Merge sequences for each individual transcript
Use transcript_sequences_extraction.sh
Merge consequtive sequences in vim by "gJ" in command mode

# Get ORFs using standard code with alternative initiation codons, min nucleotide size of 30, with ORF defined as a region that begins with a START and ends with a STOP codon, and only finding ORFs in the forward sequence (not including reverse complement - using emboss suite getorf command available on Orchestra version 6.6.0

getorf -sequence merged_transcript_sequences.fa -outseq ./potential_orfs.fa -table 1 -find 1 -reverse No 

cp potential_orfs.fa potential_orfs.fa_copy

# Collapse redundant protein fasta sequences using awk
Using the potential_orfs.fa_copy file do the following:
in vim, remove all new lines ':% s/\n/'
in vim, add new lines before > ':%s/>/\r>/g'
in vim, add new line after ] ':%s/]\s/]\r/g'
in command line, collapse duplicates
awk '!x[$0]++' potential_orfs.fa_copy > unique_potential_orfs.fa

# Remove header lines for the duplicated sequences that were already removed in text wrangler by finding and deleting

# Align sequences using blat to HIV strain NL4-3 accession number AF324493.2
blat AF324493_hiv_nl43_ref_seq.fa ../ocwieja_transcript_sequences.fa -t=dna -q=dna -out=blast aligned_transcripts_NL43.blast
