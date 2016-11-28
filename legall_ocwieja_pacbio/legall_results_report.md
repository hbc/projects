# Identification of unknown HIV proteins through re-analysis of Pacbio data

The goal of this project is to provide analysis of Pacbio RNA-Seq reads to generate a more complete database of possible HIV proteins for analysis of Mass Spectrometry data for the Le Gall group. The reads were obtained from the [SRA database](http://www.ncbi.nlm.nih.gov/sra/), which had been deposited for the publication of the paper, [Dynamic regulation of HIV-1 mRNA populations analyzed by single-molecule enrichment and long-read sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3488221/) by Ocwieja et. al (Nucleic Acids Res. 2012 Nov; 40(20): 10345–10355). 

The Ocwieja paper had identified novel splice sites, leading to identification of a new class of HIV transcripts and two confirmed new proteins for the HIV 89.6 strain (Tat8c and Ref). Within the Ocwieja paper, the genomic coordinates for splice sites identified from the Pacbio reads were given along with how the splice sites were combined to form transcripts in tables within the supplement of the paper. We performed two different analyses to help identify any possible proteins generated from HIV strain 89.6:

1. Use the splice sites given in the Ocwieja paper to reconstruct transcripts and identify open-reading frames and any possible proteins.
2. Download the Pacbio data used in the Ocwieja paper from the SRA database, then align the transcripts to the HIV genome, allowing for non-canonical splice sites. Use the aligned reads to identify open-reading frames (ORFs) and any possible proteins. We compared the proteins identified with the possible proteins from the Ocwieja paper to demonstrate the method of analysis was sufficient for identifying known and novel proteins.

Since the Pacbio data was generated from the HIV strain 89.6, and Mass Spec facility uses the NL4-3 strain, we converted the read alignment coordinates from the 89.6 strain to the NL4-3 strain, then determined ORFs and potential proteins for that particular strain.

## Analysis of splice sites identified in the Ocwieja paper

For the analysis of reconstructing transcripts from the Ocwieja paper, we manually entered the coordinates for each exon, then generated the HIV 89.6 nucleotide sequence for those genomic coordinates using BEDTools (version 2.26.0) and the reference sequence for HIV strain 89.6 (GenBank: U39362.2). The start of the first exon and stop of the last exon for each transcript were extended to begin at the beginning and end at the end of the HIV genome, respectively.  

>Quinlan, AR, Hall, IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15; 26(6): 841–842.

The sequences were output for each exon, and for processing the exons into transcripts, we combined exons with the same header information, and removed lines separating the combined sequences using the [FASTX Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) (version 0.0.13).

To identify potential ORFs, the Emboss suite of tools was used (version 6.6.0). ORFs were identified at any location in the read sequences using standard code and alternative initiation codons. The lowest minimum nucleotide size (30) was used, and ORFs were defined as a region that began with a START codon and ended with a STOP codon. We only found ORFs on the forward sequence, as no known transcripts are known to be encoded on the reverse strand for HIV. The identified ORFs were output as potential proteins.

>Rice,P. Longden,I. and Bleasby,A. EMBOSS: The European Molecular Biology Open Software Suite. Trends Genet. 2000 Jun; 16(6):276-7.

Finally, redundant proteins were removed and a total of 108 potential proteins were identified. The potential proteins were examined to ensure the analysis was able to identify all known HIV proteins and the two novel proteins identified in the Ocwieja paper. The proteins Gag and Pol are not generated from a spliced transcript; therefore, they were not listed within the splice tables identified from the Ocwieja paper and not pursued for further investigation for novel splice sites. 

All known HIV proteins were identified among the list of 108 potential proteins generated from the Ocwieja tables, including Env, Vpu, Vif, Tat, Vpr, Rev, and Nef. The two novel proteins identified in the Ocwieja paper were also present, Tat8c and Ref. 

To identify the corresponding proteins for the HIV NL4-3 strain, the genomic coordinates for the each exon identified in the Ocwieja paper was converted from the genomic coordinates for HIV 89.6 strain (GenBank: U39362.2) to the NL4-3 strain (GenBank: AF324493.2) using [KentUtils](https://github.com/ENCODE-DCC/kentUtils) liftover (version 302). Using the new coordinates the same workflow was completed as for the 89.6 strain, which included extracting the nucleotide sequence for each exon, combining exons to form transcripts, and identifying potential ORFs and proteins. The potential proteins were examined for the novel proteins identified by Ocwieja for Tat8c and Ref. A very similar protein was identified for Ref and a shorter protein was identified for Tat8c.

## Analysis of the Pacbio reads generated for the Ocwieja paper

After downloading the SRA files, the consensus Pacbio reads were extracted using the [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software) (version 2.7.0). Sequences from all samples were combined, and filtering was performed similar to the Ocwieja paper. Sequences less than 40 nucleotides and sequences matching the RT primer were removed using the FASTX Toolkit. 

After filtering, the consensus Pacbio reads were aligned to the HIV genome using [STAR](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635), a splice-aware aligner with options for aligning long reads generated by Pacbio. The human and HIV genomes were both provided during the alignment to identify and remove human contamination. To improve the likelihood of identification of novel splice junctions, no annotations were provided for human or HIV. 

Using STAR (version 2.5.2b), the STARlong mode was used and specific STAR parameters were included based on [Pacbio recommendations](https://github.com/PacificBiosciences/cDNA_primer/wiki/Bioinfx-study:-Optimizing-STAR-aligner-for-Iso-Seq-data). We also altered parameters to increase the likelihood of identifying more non-canonical splice junctions, by including splice junctions that are only supported by a single uniquely mapping read and reducing the required overhang for non-canonical splice junctions.

To also better identify any novel junctions, STARlong was used to perform 2-pass mapping, where the splice junctions identified during the alignment were used as annotations when performing an additional round of alignment. The 2-pass mapping increased the number of mapped reads; however, over half of the reads were discarded due to poor alignment of the reads. 

To summarize the alignment statistics, the total number of input reads to STAR was 512,181, with an average input read length of 789. The number of reads that were uniquely mapping was 230,631 (45.03%). The number of identified splice junctions had the following breakdown: 219425 GT/AG, 6807 GC/AG, 8882 AT/AC, and 7236 non-canonical. The mismatch rate per base was 2.53%. 673 reads mapped to multiple loci (0.13%), and 842 mapped to too many loci (0.16%) and were dropped. Most of these statistics are quite good, and the main concern with the alignment is the loss of over half of the reads identified as 'too short' (53.37%) based on the STAR default requiring mapped length to be > 2/3 of the total read length. This step may require optimization in the future if this procedure for identifying potential proteins is continued.

>Dobin, A., et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1; 29(1): 15-21.

The reads aligning to the human genome were removed using [SAMtools](http://samtools.sourceforge.net) (version 1.3). After removing human contamination, 230,016 reads of the 233,262 total aligned reads remained. 

>Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009). The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
>
>Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. Epub 2011 Sep 8. [PMID: 21903627]

To generate potential open reading frames, the nucleotide sequences were extracted. A BED file was required to generate the FASTA sequences, so the coordinates were first converted to BED format using BEDtools prior to extraction. To combine the exons into transcripts, exons with the same header information were combined and lines separating the combined sequences were removed using the FASTX Toolkit.

ORFs and potential proteins were identified similar to the procedure used in the previous analysis. A total of 463,629 potential proteins were identified (many redundant due to the same protein being generated by multiple reads), which corresponded to a total of 7,797 unique protein sequences. To determine whether the potential proteins from the Pacbio analysis corresponded to the potential proteins identified from the transcript coordinates given in the Ocwieja paper, the Pacbio proteins were compared to those in the paper using BLAT (version 35) with a filter for minimum identity equal to 100%.

>Kent WJ. BLAT - the BLAST-like alignment tool. Genome Res. 2002 Apr;12(4):656-64.

Of the total 7,797 unique potential proteins identified by the Pacbio reads, 1,089 of the proteins aligned completely to the Ocwieja potential proteins. Therefore, the remaining 6,708 potential proteins include potential novel HIV proteins, in addition to false positives (proteins with incorrectly identified splice junctions and ORFs not actually translated in HIV). 

Of the 108 potential proteins identified from the Ocwieja paper coordinates, 94 were identified within the Pacbio potential proteins. The 14 potential proteins not identified by the Pacbio data corresponded to very short proteins between 10-15 amino acids in length. Of the 463,629 potential proteins identified by the Pacbio data, 149,239 were identified by BLAT as aligning to at least some portion of the 108 Ocwieja potential proteins. Many of the 149,239 potential proteins aligned to more than a single protein in the 108, with only 54,411 potential proteins aligning uniquely to a single Ocwieja protein. 

A Pacbio protein that aligned in it's entirety (full length) to an Ocwieja protein is more likely to have originated from the designated HIV Ocwieja protein than a Pacbio protein that did not completely align. Uniquely aligning full-length Pacbio proteins identified 51 of the 108 Ocwieja proteins. Descriptive statistics describing the Pacbio proteins can be found in the table below:


|               | # Pacbio proteins Aligning to Ocwieja Proteins | # Ocwieja Proteins Identified |
| ----------------------------- |:------------------------------------:|:-----------------------------:|
| Total Pacbio proteins identified as aligning to Ocwieja proteins | 149,239 / 463,629 | 94 / 108 |
| Uniquely aligning full-length Pacbio proteins | 49,797 / 463,629 | 51 / 108 |
| Total uniquely aligning Pacbio proteins | 54,411 / 463,629 | 56 / 108 |


The Pacbio potential proteins were also examined to ensure the analysis was able to identify all known HIV proteins and the two novel proteins identified in the Ocwieja paper. All known HIV proteins were identified by uniquely aligning full-length Pacbio potential proteins, including Env, Vpu, Vif, Tat, Vpr, Rev, and Nef. The two novel proteins identified in the Ocwieja paper were also present, Tat8c and Ref. 

To be conservative, we estimated only full-length, uniquely aligning Pacbio proteins to determine the number of Pacbio reads aligning to each HIV protein. The numbers of unique full-length Pacbio proteins aligning to each of the known HIV proteins are displayed in the table below:

| # HIV Protein | # Reads Resulting in Uniquely Aligning Pacbio Full-Length Proteins |
|:-------------:|:-----------------------------------:|
| Vif | 126 |
| Vpr | 116 |
| Tat | 651 |
| Rev | 4 * |
| Vpu | 2,545 |
| Env | 3,206 |
| Nef | 1,263 |
| Tat8c | 12 |
| Ref | 295 * |

_* The Rev and Ref proteins do not have uniquely mapping reads because both full-length proteins are entirely contained within different full-length potential proteins identified during the Ocwieja analysis. The numbers for the Rev and Ref proteins reflect the number of full-length Pacbio potential proteins that align to the entire Rev/Ref protein, without containing any additional sequence._

The low number of full-length uniquely mapping reads for Rev are likely due to the fact that the majority of the protein overlaps the Tat protein. Therefore, it is likely that many reads may have originated from Rev or Tat, but were not counted since they were not unique to either protein.

Finally, to identify the corresponding proteins for the HIV NL4-3 strain, the genomic coordinates for the each exon identified during the Pacbio analysis was converted from the genomic coordinates of the HIV 89.6 strain (GenBank: U39362.2) to the NL4-3 strain (GenBank: AF324493.2) using KentUtils liftover similar to the Ocwieja analysis. Using the new coordinates the same workflow was completed as for the 89.6 strain to identify potential ORFs and proteins. 
