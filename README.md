# DeltaMP
 **A flexible, reproducible and resource efficient metabarcoding amplicon pipeline for HPC**

DeltaMP is a command line tool for high performance computer taking advantage of queueing systems to parallelize and standardize bioinformatic processing of metabarcoding raw read libraries.

DeltaMP is initially developed to process **16S, 18S, 28S, ITS, COI and rbcl** raw read libraries with the most up-to-date bioinformatic workflows, but can also handle any other barcoding targets (e.g. 23/28S, rbcL).

DeltaMP intend to be accessible for non-bioinformatician users with its fully tunable workflows based on a TAB-separated configuration file.

DeltaMP integrate a checkpointing feature which enable easy and efficient comparisons between different workflows applied on the same set of read libraries.

Last but not least, DeltaMP produces version controlled, reproducible and fully documented OTU tables in TAB-separated and BIOM formats readily usable for downstream taxonomic and OTU diversity analyses.

## Table of Contents
 - [Installation](#installation)
 - [Dependencies](#dependencies)
 - [Quick start](#quick-start)
 - [Usage instructions](#usage-instructions)
 - [Configuration file](#configuration-file)
   * [PROJECT section](#project-section)
   * [LIBRARIES section](#libraries-section)
   * [TARGET section](#target-section)
   * [PAIR-END section](#pair-end-section)
   * [TRIMMING section](#trimming-section)
   * [PIPELINE section](#pipeline-section)
   * [SAMPLES section](#samples-section)
 - [Pipeline execution](#pipeline-execution)
   * [Quick start](#quick-start)
   * [best practices](#best-practices)
   * [pipeline execution options](#pipeline-execution- options)
   * [Checkpointing](#checkpointing)
 - [Pipeline analysis steps](#pipeline-analysis-steps)
 - [Outputs](#outputs)
 - [Troubleshooting](#troubleshooting)
 - [References](#references)


## Installation

Source code of DeltaMP version 0.5 is available at <https://github.com/lentendu/DeltaMP/releases/tag/v0.5>

The repository could also be cloned using for example:
```
git clone https://github.com/lentendu/DeltaMP.git
```

Then install following the [installation instructions](INSTALL).


## Dependencies

DeltaMP is intend to be used on a HPC with a job scheduler (i.e. batch-queuing system).

### supported job schedulers:

**[SLURM](https://slurm.schedmd.com)**

**Grid Engine**
+ [Oracle Grid Engine](http://www.oracle.com/technetwork/oem/grid-engine-166852.html)
+ [Open Grid Scheduler](http://gridscheduler.sourceforge.net)
+ [Univa Grid Engine](http://www.univa.com/products)


### compulsory softwares:

+ [vsearch v2](https://github.com/torognes/vsearch) (Rognes et al., 2016)
+ [MOTHUR v1.44+](http://www.mothur.org) (Schloss et al., 2009)
+ [SeqKit v0.15.0+](https://github.com/shenwei356/seqkit)
+ [cutadapt v1.10+](https://cutadapt.readthedocs.io/en/stable/) (Martin, 2011)
+ [PANDAseq v2.10+](https://github.com/neufeld/pandaseq) (Masella et al., 2012)
+ [biom-format](http://biom-format.org) (McDonald et al., 2012)
+ [GNU parallel](http://savannah.gnu.org/projects/parallel)
+ [R v3](https://cran.r-project.org/) with the plyr, tidyverse, cowplot, foreach, doParallel and seqinr packages
+ [Ghostscript](https://www.ghostscript.com/)
+ [WebLogo](https://github.com/WebLogo/weblogo)
+ [archive-sum](https://github.com/idiv-biodiversity/archive-sum)


### optional softwares:

+ [QIIME v1.8+](http://qiime.org) for 454 demultiplexing ([minimal install](http://qiime.org/1.9.0/install/install.html#native-base) is enough; Caporaso et al., 2010)
+ [FlowClus](https://github.com/jsh58/FlowClus) for denoising of 454 flows (Gaspar and Thomas, 2015)
+ [dada2](https://github.com/benjjneb/dada2) R package for error model based correction of Illumina reads (Callahan et al., 2016)
+ [NGmerge](https://github.com/harvardinformatics/NGmerge) for paired-end merging of Illumina reads via novel empirically-derived models of sequencing errors (Gaspar, 2018)
+ [swarm v2](https://github.com/torognes/swarm) (Mahé et al., 2015)
+ [sumatra](https://git.metabarcoding.org/obitools/sumatra) and [sumaclust](https://git.metabarcoding.org/obitools/sumaclust) (Mercier et al., 2013)
+ [MCL](http://micans.org/mcl) (van Dongen, 2000)
+ [CD-HIT & cd-hit-454](http://weizhongli-lab.org/cd-hit/) for (pre-)clustering (Fu et al., 2012)
+ [ITSx](http://microbiology.se/software/itsx/) (Bengtsson‐Palme Johan et al. 2013)
+ [IBM Aspera Connect](https://downloads.asperasoft.com/connect2/) for faster download of raw reads libraries from the ENA SRA public database

All this dependencies need to be available through the $PATH environmental variable or need to be loaded by the DeltaMP module file.


## Quick start




## Usage instructions

Usages could be accessed via the command line:

```
deltamp -h
```

```
NAME
	DeltaMP version 0.5 - a flexible, reproducible and resource efficient metabarcoding amplicon pipeline for HPC

SYNOPSIS
	Usage: deltamp [-a account_or_project_name] [-cdfhnqtx] [-m max_running_tasks] [-p reference_subproject] [-r step] configuration_file

DESCRIPTION
	-h	display this help and exit

	-a ACCOUNT
		account or project name for the job queuing system

	-c	check dry run (-d option) by preserving the SUBPROJECT and PROJECT directories in the output directory

	-d	dry run avoid any submission to the queuing system and only output submission informations

	-f	display default values for optional configuration parameters

	-n	avoid checkpointing (i.e. searching for previous SUBPROJECT) and run the pipeline from the beginning

	-m max_running_tasks
		fix a maximum number of concurrently running tasks per array jobs, default is 400

	-p SUBRPOJECT
		only execute additional jobs or with different input variables compared to the reference SUBPROJECT. This will only work if both configuration files have the same input libraries and
		output directory.

	-q	proceed until quality step only

	-r STEP
		restart pipeline computation from STEP. Replace STEP by 'list' to list all available steps of the subproject associated with the provided configuration file.

	-t	proceed until demultiplexing step only

	-x	delete the subproject associated with the provided configuration file

AUTHOR
	Guillaume Lentendu, Christina Weißbecker, Anna Heintz-Buschart and Tesfaye Wubet

REPORTING BUGS
	Submit suggestions and bug-reports at <https://github.com/lentendu/DeltaMP/issues>, send a pull request on <https://github.com/lentendu/DeltaMP>, or compose an e-mail to Guillaume Lentendu <guilaume.lentendu@unine.ch>.

COPYRIGHT
	Copyright (C) 2018 Guillaume Lentendu, Christina Weißbecker, Anna Heintz-Buschart and Tesfaye Wubet

	This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
```

## Configuration file

To run, DeltaMP always required a configuration file and, optionaly, local copies of raw sequence libraries.

The configuration is encoded in UTF8 and the fields are TAB-separated.

If the file is edited under Windows, pay attention to use Unix compliant end of line (\n).

The following sections describe the parameters to provide in the configuration file.

Only the project name parameter and a sample list (see [SAMPLE section](#sample-section)) are compulsory.

A parameter is recognize by its exact definition in the first column of the configuration file.

If an optional parameter is let empty or if its definition is not found in the first column, its default value will be used instead.

The default values for the optional parameters could be displayed by using the -f option (see the [usage instructions](#usage-instructions)).

### PROJECT section

+ __Project name__: name of your project. This can be shared among multiple datasets and targets.

+ __User name__: full name of the person in charge of the analysis. The default is $USER.

+ __Path to output location__: full path to the output location. All final outputs will be copied into sub-directories of that directory. The default is /home/$USER .

+ __Path to execution location__: full path of the working directory. All computations will take place in that directory. The default is /work/$USER . 


### LIBRARIES section

+ __Sequencing technology__: “Illumina” or “454”. The default is Illumina.

+ __Directory path to archives or libraries OR BioProject accession__: full path of the directory where to find the sequence libraries or archives. If the raw sequences are publicly available at the European Nucleotide Archive (https://www.ebi.ac.uk/ena), only provide the BioProject accession number. The libraries will then automatically be retrieved from the provided URLs in the BARCODE section. The default is /home/$USER.

+ __Tab separated list of archives__: filename of archives containing raw libraries, one per column. The default is “no” for no archive.

+ __Libraries already demultiplexed__: “yes” or “no”. If the libraries are not demultiplexed, one fastq/sff file per sample in each library will be produced. An archive containing all demultiplexed libraries and labeled with "raw_fastq"/“raw_sff” will be created in the output directory. The default is “yes” for libraries already demultiplexed.

+ __Strict barcode anchoring__: “yes” or “no”. Wether the barcodes are strictly anchored at the 5'-end or partial match are also allowed (non-internal barcodes) in raw undemultiplexed libraries. The default is "yes".

+ __Single barcode side__: *Illumina specific parameter*, "both", "forward" or "reverse". For single-indexing (one barcode per library), provide the position of the barcode, "both" for both forward and reverse primers, "forward" for barcode only associated with the forward primer, "reverse"  for barcode only associated with the reverse primer. The default is “both”.

+ __Bind barcode to primer for demultiplexing__: *Illumina specific parameter*, "yes" or "no". For dual-indexing, control binding of barcodes with the forward/reverse primer to increase demultiplexing accuracy. The first barcode in a pair is only search associated with the forward primer at 5'-end and the second barcode is only search associated with the reverse primer at 5'-end, in both paired libraries. The default is “no”.

+ __Orientation threshold__: *Illumina specific parameter*, a percentage of raw demultiplexed reads above which to keep using reads in a certain orientation (forward primer in R1 and reverse primer in R2, and/or reverse primer in R1 and forward primer in R2). The default is 1.

### TARGET section

+ __Target organisms__: Latin name of the targeted group as labelled in the reference sequence database (eg. Fungi, Glomeromycota, Bacteria, Archea, Eukaryota), or protists (Eukaryota excluding Fungi, Metazoa and Streptophyta). The OTUs assigned to the target group will be extracted from the main OTU table and copied into an additional OTU table labeled with the provided name (see [Outputs](#outputs)). The default is “protist”.

+ __Target region__: any name is accepted here, but specific workflows will be allowed with “16S”, “18S”, "28S", “ITS”, “COI” or "rbcl" respectively for prokaryotes, eukaryotes (18S and 28S), Fungi, Metazoa or plants workflows. The default is “18S”.

+ __Forward primer name__: the biological forward primer. The 5'-3' orientation of the forward primer must be the same than the reference database sequence orientation. The default is “TAReuk454FWD1”.

+ __Forward primer sequence (5' to 3')__: the sequence of the biological forward primer. The default is “CCAGCASCYGCGGTAATTCC”.

+ __Reverse primer name__: the biological reverse primer. The default is “TAReukREV3”.

+ __Reverse primer sequence (5' to 3')__: the sequence of the biological reverse primer. The default is “ACTTTCGTTCTTGATYRA”.

+ __Sequencing direction__: *454 specific parameter*, “forward”, “reverse” or "both" for libraries (column 3 of [SAMPLES section](#samples-section)) predominently orientated in the forward or reverse primer direction or for bi-directional sequencing, respectively. The default is “forward”.

+ __Remove primers at reads ends__: "no", "3prime", "5prime" or "both". Using cutadapt, "3prime" remove non-anchored reverse-complemented reverse and forward primers from 3'-end of forward and reverse libraries, respectively; "5prime" remove non-anchored forward and reverse primers from 5'-end of forward and reverse libraries, respectively; "both" combined both previous removal using the [linked adapter strategy](https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters-combined-5-and-3-adapter) with anchored 5' adapter; value "no" avoid primer clipping. For "3prime", "5prime" and "both" values, unmatched reads are checked against the opposite primer/linked-adapter to allow for biderectionnal sequencing strategy detection, and the reads matched in this second step are reverse-complemented and append to the reads match at the first step. The default is "5prime".

+ __Anchored primers__: “yes” or “no”. Wether the primers are anchored to the 5'-end or not in demultiplexed libraries. The default is "yes".


### PAIR-END section
*all Illumina-specific parameters*

+ __Pair-end algorithm__: "vsearch", "ngmerge", "simple_bayesian", "ea_util", "flash", "pear", "rdp_mle", "stitch" or "uparse" . Either using the vsearch "--fastq_mergepairs" option, the NGmerge tool or any algorithm of PandaSeq. See the [PandaSeq manual](https://storage.googleapis.com/pandaseq/pandaseq.html) (or `man pandaseq`) for algorithm descriptions. The default is “simple_bayesian”.

+ __Pair-end similarity threshold__: a number between 0 and 1 setting a similarity threshold ratio on the alligned region from which a pair is conserved. This value is provided to the -t option of pandaseq. This value is substract to 1 and then multiply by 100 to be provided to the --fastq_maxdiffpct option of vsearch. The default is 0.6 (i.e. vsearch default is 40).

+ __Minimum overlap__: a number between 0 and ∞ provided to the --fastq_minovlen option of vsearch or the -o option of pandaseq. The default is 10.

+ __Keep unmerged__: “yes” or “no”. Wether to keep high-quality reads which cannot be pair-end assembled alongside with paired reads for further processing or not. A padding sequence of "NNNNN" is added between the unpaired reads. The default is "no".

+ __Avoid pair-end__: “yes” or “no”. In case of low quality reads or too long target fragment, the pair-end assembly procedure is skipped and paired sequences are simply joined with a padding sequence of "NNNNN". Use this option as a last resort in case all other parameter tunning did not provide suitable results. The default is "no".


### TRIMMING section
+ __Denoising__: *454 specific parameter*, "yes" or "no", whether to perform 454 flowgram denoising using FlowClus or nothing. The default is "no".

+ __Minimum number of flows__: *454 specific parameter*, a number between 100 and 600, minimum length of flowgrams to be kept for denoising. The default is 360.

+ __Number of mismatches allowed on the barcode sequence__: a number between 0 and 2, or 'a' for automatic detection of maximum allowed mismatch on a barcode to avoid mislabeling of sequence inside one library. This value is used for demultiplexing only. The default is 1.

+ __Number of mismatches allowed on the primer sequence__: a number between 0 and 10. For Illumina libraries, this number will only be used to calculate allowed dissimilarity (number of mismatch / primer length) for primer detection and removal by cutadapt. The default is 6.

+ __Maximum number of ambiguities allowed in the sequence__: a number between 0 and ∞. The default is 0.

+ __Maximum homopolymer length allowed__: a number between 0 and 20. The default is 10.

+ __Minimum sequence length__: a number between 50 and the set Maximum sequence length. This parameter is used for filtering pair-end assembled reads in the case of Illumina libraries. The default is 50.

+ __Maximum sequence length__: a number between the set Minimum sequence length and 1000. This parameter is used for filtering pair-end assembled reads only in the case of Illumina libraries. The differrence between the __Maximum sequence length__ and the __Minimum sequence length__ added to the __Minimum overlap__ defines the maximum overlap length allowed for PandaSeq based pair-end assembly. The default is 600.

+ __Type of quality filtering__: "average" or "maxee". Filter sequences either based on the average quality or based on the maximum expected error over the sequence, after any length truncation. Filtering based on maximum expected error is only possible for Illumina reads. The default is "average".

+ __Truncation before pair-end__: *Illumina specific parameter*, "yes" or "no", only for available for _Type of quality filtering__ set to "maxee". Perform length truncation on unpairred reads using a maximum expected error, and eventually a tuncation length in order to improve pair-end assembly efficiency. The default is "no".

+ __Minimum average quality on the trimmed sequence length__: a number between 20 and 30 (Phred score). Only consider when the parameter __Type of quality filtering__ is set to "average". The read quality average is calculated after optimal sequence length trimming. The default is 20.

+ __Maximum expected error on the trimmed length__: *Illumina specific parameter*, a number between 0.5 and ∞, by step of 0.5. Only consider when the parameter __Type of quality filtering__ is set to "maxee". The default is 2.

+ __Minimum value for the maximum expected error on the trimmed length__: *Illumina specific parameter*, a number between 0.5 and ∞, by step of 0.5. Only consider when the parameter __Type of quality filtering__ is set to "maxee". Only use for trimming length and quality optimization The default is 0.5.

+ __Increment for maximum expected error optimization__: *Illumina specific parameter*, a number between 0.1 and 1, which is a value by which the __Maximum expected error on the trimmed length__ is incremented for the optimization procedure. The default is 0.5.

+ __Minimum length truncation of unpaired reads__: *Illumina specific parameter*, "no" or a number between __Minimum sequence length__ / 2 and ( __Maximum sequence length__ + __Minimum overlap__ ) / 2. When a length is provided, a length truncation is performed before pair-end join, using the length maximizing the length of the read while minimizing the maximum expected error on the truncated length, which allows to keep the read counts fixed by the __Minimum number of trimmed reads per sample__ parameter. The truncation length and the maximum expected error are optimized for all combinations of libraries (R1 and R2) and directions (forward and reverse primers) separately, but are fixed for all reads in a combination. To allow for read truncation only based on maximum expected error (thus producing reads of variable sizes), set to "no". The default is "no".

+ __Maximum length truncation of unpaired reads__: *Illumina specific parameter*, "no" or a read length __Minimum sequence length__ and __Maximum sequence length__ , and higher than __Minimum length truncation of unpaired reads__ + __Increment for length truncation optimization__ .

+ __Increment for length truncation optimization__: *Illumina specific parameter*, a number between 1 and 10, which is the number of nucleotide by which the __Minimum length truncation of unpaired reads__ is incremented for the optimization procedure. The default is 5.

+ __Expected mean length of amplified barcode gene__: *Illumina specific parameter*, "no" or a read length between __Minimum sequence length__ and __Maximum sequence length__ use to optimize truncation length and maxEE parameters, so that assembled pair-end reads reach at least this length

+ __Reuse previous subproject optimized quality parameters__: *Illumina specific parameter*, "no" or a path to a previous suproject execution directory. This will reuse the same optimize quality parameters for length and quality fitering truncation if the primers and the options __Type of quality filtering__, __Truncation before pair-end__ and __Minimum length truncation of unpaired reads__ are identical. The default is no.


### PIPELINE section

+ __Minimum number of trimmed reads per sample__: a number between 0 and ∞. If the minimum of raw, raw with primers or trimmed reads is not reached in any of the samples, the pipeline will stop after the quality step and report these reads counts in the quality step .out logfile. If the value is between 0 and 1, no threshold is apply on raw reads and the value is used as ratio of reads in the previous filtration step for the tresholding of raw with primers, pair-end and trimmed reads (e.g. read count of trimmed reads have to be > to the provided ratio times read count of pair-end reads). The default is 2000.

+ __Minimum number of pair-end reads per sample__: *Illumina specific parameter*, a number between 0 and ∞. If the minimum number of pair-end reads is not reached in any of the samples, the pipeline will stop after the quality step and report these reads counts in the quality step .out logfile. If the value is between 0 and 1, no threshold is apply on raw reads and the value is used as ratio of reads in the previous filtration step for the tresholding. The default is 0.8.

+ __Samples percentage allowed below the read count threshold__: percentage of samples which are allowed to have their read count below the threshold set by the __Minimum number of trimmed reads per sample__ parameter. These samples below the threshold are not considered for quality, maxee and length truncation optimization. A value of 100 percent disable thresholding and quality parameter optimizations and use the minimum (average quality, length truncation) or maximum (expected error) values provided for quality filtering. The default is 0.

+ __ITSx region to extract__: "no", “ITS1” or “ITS2” to skip or to extract either regions from fungal ITS reads using ITSx. The default is “no”.

+ __Pre-clustering__: "no", "cdhit454" or "mothur" to skip or to use [cd-hit-454](http://weizhong-lab.ucsd.edu/public/?q=softwares/cd-hit-454) or [mothur pre.cluster](https://mothur.org/wiki/pre.cluster/) algorithms to pre-cluster reads after trimming and before chimera removal. For the mothur based pre-clustering, and aligned version of a reference database is required (see below at __Database prefix name__). The default is "no".

+ __Chimera removal__: "before", "after" or "both", to check for chimera before OTU clustering only (in each sample separetedly), after OTU clustering only (among OTU representative sequences), or at both moments. *De-novo* chimera are detected with UCHIME and removed. For dada2 clustering, bimera are always removed using _removeBimera_. The default is "after".

+ __Remove bimera after ASV calling and joining__: *Illumina specific parameter*, "no" or "yes" to remove bimera when clustering with dada2. Independent from the __Chimera removal__ parameter. The default is "no".

+ __Subsampling__: “yes” or “no”. Activating subsampling will randomly select the same number of reads in all samples according to the read count in the sample with the lowest number of reads, after the trimming step. The default is “no”.

+ __Clustering algorithm__: "cd-hit-est”, “vsearch”, “mcl”, “sumaclust”, “swarm” or “dada2“. The default is “mcl”.

+ __Clustering similarity threshold__: a number between 80 and 100 (percent of similarity). Disregarded if the “swarm” or "dada2" clustering algorithms are used. The default is 97.

+ __Single run__: "yes" or "no", consider all libraries coming from the same run in order to build a single error model for dada2. The default is "no".

+ __Cluster with previous subproject reference sequences__: the accepted values are a full path or “no”. If the full path to a previous subproject output directory is provided, all amplicons of this previous subproject will be search for exact match with amplicons from the current subproject. If amplicons of one or multiple OTUs from the current subproject have match in a single OTUs of the previous subproject, these amplicons will be assign to the OTU index name used in the previous subproject. Remaining OTUs from the current subproject without match in the reference suproject have their OTU index names starting after the last OTU of the reference subproject. This option is experimental, use with caution. The default is “no”.

+ __Cluster with previous subproject percent overlap__: a percentage of length overlap between current and previous subproject amplicons. This is usefull to merge OTUs of two subproject when different primers were used to amplify the same region. The default is 95.

+ __Minimum abundance of previous subproject amplicons to cluster with__: the minimum abundance of any dereplicated and non-chimeric amplicon in the previous and the current subprojects to be included in the clustering with previous subproject approach. The default is 4.

+ __Remove singletons before chimera re-check__: the accepted values are “yes” or “no”. This option allow to remove singleton OTUs from sub-sequent analyses in order to reduce computation time and memory footprints, as high number of singletons could have critical effects on grid job with limited resources for both following chimera re-check and taxonomic identification steps. The default is “no”.

+ __Taxonomic classifier__: "bayesian" or "vsearch", to use the RDP naive bayesian classifier (Wang et al., 2007) as implemented in mothur or the usearch global alignment method (Edgard, 2010) as implemented in vsearch. The default is "bayesian".

+ __Directory path to database__: the full path to the directory containing reference sequence database files. The default is /home/$USER .

+ __Database prefix name__: filename prefix identifying database files. A database have to be composed from a TAB-separated "prefix".txt file in which the second column describe the version, citation and reference of the database with the first column containing VERSION, CITATION and FULLCITATION, respectively (see e.g. <test/pr2_4.10.0.txt>) and either from a pair of "prefix".fasta and "prefix".taxonomy mothur's style database format for "bayesian" classifier, or from a "prefix".udb binary file produced by the command `vsearch --makeudb_usearch` for "vsearch" classifier. Optionnaly, a "prefix".align.fasta, containing an aligned version of the "prefix".fasta sequences, is needed for pre-clustering with "mothur". See the <test/> directory for example files in each format. You can use the "update_xxx.sh" scripts in <auxillary_scripts/> to automatically download the PR2, SILVA and UNITE databases and to format them into the previously described formats. The default is "pr2_4.10.0".

+ __Reduce database to amplified fragment__: "yes" or "no", to cut the database reference sequences using either the provided primers or ITSx with the provided 'ITSx region to extract' or do nothing, respectively. If set to "yes", this will use the cut database for taxonomic assignment, while if set to "no" this will use the full database sequences for taxonomic assignment. For SILVA and UNITE databases, if a unique cut read is produced from multiple accessions with different taxonomic path, the cut read will be annotated with the least common ancestor. For PR2 database, cut reads are only dereplicated for each taxonomic path separately. If set to "yes" and the cut database is already present in the "Directory path to database", the database cutting is skipped. The default is "no".

+ __Consensus assignment threshold__: a number between 50 and 100. Consensus threshold percent to assign a taxonomic rank among the matches. For "bayesian" taxonomic calssifier the matches are the 100 bootstrap matches, for "vsearch" the matches are the best match(es) plus the match(es) in a 1 % similarity range below the best match(es) similarity. The default is 60.

+ __Assign all reads__: the accepted values are “yes” or “no”. Activating this option will assign all dereplicated or pre-clustered reads to a taxonomy. A consensus assignment is then determined for each OTU at a threshold of 60 %. If set to “no”, only the most abundant read per OTU will be assigned to taxonomy. The default is “no”.

+ __Full path to additional database prefix name__: the full path to the filename prefix of additional reference sequences to add to the reference database (useful to integrate yet unpublished sequences). The fasta and taxonomy formats need to match with the main database. dub format for vsearch will be generated automatically if needed. The default is "no".

+ __Assign putative function__: "yes" or "no" to assign putative function to Fungi using the FUNGuild database (Nguyen et al., 2016).

+ __Minimum number of sample for abundant OTUs__: a number between 0 and ∞ (see [Outputs](#outputs)). The default is 1.

+ __Minimum number of reads for abundant OTUs__: a number between 0 and ∞ (see [Outputs](#outputs)). The default is 4.


### SAMPLES section

+ __454 libraries__: the columns 1 to 3 have to be filled with barcode sequences, sample names and library filenames or URLs, respectively. If the same barcode was used to sequence in both forward and reverse directions, the library containing the reverse primer at its 5'-end have to be provided in the fourth column.

+ __Illumina libraries__: the columns 1 to 4 have to be filled with barcodes, sample names, the filename or URL of the forward libraries and the filename or URL of the reverse libraries, respectively. The column 1 can remain empty when libraries are already delmultiplexed. Dual index demultiplexing is turn on when a comma separated pair of barcodes is provided for every sample in column 1.

Columns 3 and 4 accept libraies in fastq or sff format with most kinds of compression (.gz, .tar, .tar.gz, .tgz, .bz2, .tar.bz2, .tbz2, .zip).

For ENA libraries, the column 2 have to match ENA "Submitter's sample name" field and column 3 to 4 have to match ENA full ftp or fasp (for aspera connect download) URLs of a run accession as listed in ENA fields “Submitted files (FTP)” or "FASTQ files (FTP)" or “Submitted files (Aspera)” or "FASTQ files (Aspera)".

Example configuration files 'configuration_xxx.tsv' are available in the test/ directory after installation with `make`.


## Pipeline execution

### Quick start

To execute the full pipeline, run the following command in a terminal:
```
deltamp [path/]configuration_file
```
replacing [path/] by the path to your configuration file, and “configuration_file” by its filename. Avoid “path/” if you are already in the right directory.

For each DeltaMP execution, a "PROJECT” directory labeled with the project name will be created into the "Path to output location" (nothing done if this directory already exist). A “SUBPROJECT” directory will be created inside the PROJECT directory and will be labeled as follow: “Project name”\_”Sequencing technology”\_”Target region”\_”unique identifier”. The unique identifier (cksum of date, time and configuration file) enables differentiation between each instance of the pipeline execution for the same project. 

The SUBPROJECT directory is copied into the "Path to execution location" and all necessary jobs are submited to the queue.

Once all step's jobs are completed, the outputs are copied back into the "Path to output location"/PROJECT/SUBPROJECT directory.

### Directories structure

+ output SUBPROJECT directory:
   * from stratup: directory config/ containing the initial configuration file renamed configuration.SUBPROJECT.tsv and internal configuration text files
   * at execution end: files SUBPROJECT.documentation.txt and archives SUBPROJECT.outputs.tar.gz and SUBPROJECT.processing.files.tar.gz (and respective md5sum)
   
+ execution SUBPROJECT directory:
   * archives/: store outputs before copying them back to output directory
   * config/: same as in output directory with additionnaly "jobid" list of jobids and jobnames of queued jobs and lists of files ("xxx.files") and directories ("xxx.dir") in the execution directory after each job step "xxx"
   * libraries/: for raw libraries and reads filtering until quality step
   * log/: output "xxx.out" and error "xxx.err" logfiles of each job step "xxx", split among tasks for array job (e.g. "xxx.1.out")
   * processing/: contains files and directories produced from trim step to the end
   * quality_check/: used for the quality step

### best practices

+ Begin with a dry run to check for correct configuration
```
deltamp -d [path/]configuration_file
```

+ start analyses with a "quality_only" workflow in order to first control read counts and quality:
```
deltamp -q [path/]configuration_file
```
The output statistics will be then useful to check the proper handling of raw reads until trimming step (e.g. demultiplexing of 454 libraries or pair-end assembly of Illumina libraries), before starting more computational intensive OTU clustering and taxonomic identification steps.

+ test DeltaMP behaviors for checkpointing with multiple dry run on a continuously modified configuration files, avoiding the cleaning of the dry run directories:
```
deltamp -cd [path/]configuration_file
```
To run the pipeline normally after such a dry checkpointing, first delete all directories created for dry runs from the output directory, or change the output directory in the configuration file. 

### Checkpointing

In order to avoid to run the pipeline from the beginning if only one parameter is modified in the configuration file (e.g. clustering algorithm), DeltaMP allow re-use of previously produced files in the execution directory by [checkpointing](https://en.wikipedia.org/wiki/Application_checkpointing).

Checkpointing is turn on by default, which means that if any SUBPROJECT with the same target, same sequencing technology and the same list of samples is found in the "Path to output location/PROJECT" directory, all configuration parameters and option values will be compared.

Then come different situations:
+ if a previous SUBPROJECT have the same set of configuration parameters and option values than for the current deltamp call, and error message is issued and nothing appends to avoid any useless re-calculation. If the users really want to re-compute, he/she will have to first delete the previous SUBPROJECT, to change the "Path to output location" or to disable checkpointing with option -n.

+ if a previous SUBPROJECT have the same set of configuration parameters and option values than for the current deltamp call, but the SAMPLES section differs, and error message is issued and nothing appends to avoid checkpointing for different sets of samples/libraries. If the users really want to re-compute, he/she will have to first delete the previous SUBPROJECT, to change the "Path to output location" or to disable checkpointing with option -n.

+ if at least one configuration parameter or option value differs with any previous SUBPROJECT, a new SUBPROJECT will be created and the computation will start directly after the last common step, while all files and directories produced during the common steps of the previous SUBPROJECT are symlinked to the execution directory or hard copied to the output directory of the new SUBPROJECT. If multiple previous SUBPROJECT exist, the previous SUBPROJECT is the one with the highest amount of common steps (then the oldest if still multiple previous SUBPROJECT). 

Relationship tree among successive checkpointed SUBPROJECTs are documented in the output directory of SUBPROJECTs with at least one previous SUBPROJECT at config/tree.summary.

Checkpointing could be repeated multiple times, so long at least one parameter differ with any of the previous SUBPROJECTs.

The previous SUBPROJECT could be hard set with the -p option, otherwise the previous SUBPROJECT with highest number of common steps will be taken as reference.

Checkpointing can be turn off by using the -n option.


## Pipeline analysis steps

Running the DeltaMP command will generate the directories and configuration files necessary to conduct the pipeline analysis as well as submitting the required steps to the queueing system.

The following steps are bash scripts available in the bin directory after build and are named following the scheme “xxx.sh”, xxx being the name of the step.

For 454 or Illumina specific steps, step script filenames follow “454_xxx.sh” and “Illumina_xxx.sh” schemes, respectively.

+ init: create symlink to files of previous SUBPROJECT if checkpointing; one serial job
+ get: raw data copy/download and oligo file(s) creation; one parallel job
+ 454 libraries specific steps
   * demulti: demultiplex sff based on barcodes; target raw sequences in libraries containing multiple target primers are identified when holding the expected forward sequencing primer with a maximum number of mismatches equal to one third of the primer length; one parallel array job per run/lane library
   * sff: sff to fasta + qual and raw reads statistics; one parallel array job per group of 10 libraries
   * opt: sequence count for incremented length and quality values; one parallel array job per group of 10 libraries
   * quality: control sequencing depth for raw and trimmed reads; optimize the length and quality trimming parameters; one serial job
+ Illumina libraries specific steps
   * demulti: demultiplex fastq based on barcodes using cutadapt; one parallel array job per run/lane library
   * fastq: cut primers; truncate; raw reads quality and length statistics; one serial array job per library
   * pair_end: pair-end assembly; convert to fasta + qual; cut and pair-end reads counts and pair-end reads quality and length statistics; one serial array job per library
   * opt: optimization of reads average quality for trimming or reads length and maxee for truncation; one serial array job per library
   * quality: control sequencing depth for raw, cut, pair-end and trimmed reads; optimize the quality trimming parameters; one serial job
+ cut_db: reduce database reference sequences to the amplified region or to the ITS-covered region; format to mothur and vsearch udb; only one time per database version and pair of primer names; one parallel job
+ trim: group each sample reads if multiple libraries for one sample; denoise; trim or truncate; per library ASV inference; dereplicate; subsample; align and remove badly aligned sequences; precluster; remove chimeras; one parallel array job per sample
+ asv: for dada2 denoising of Illumina reads, when the samples could be group by run or library, this step perform the first step of the ASV inference by building one error model per run/library; one parallel array job per library/run
+ OTU: cluster sequences into OTUs or infer ASVs; pick representative sequences; remove singletons (optional); chimera check; one parallel job
+ id: classify all or only OTU representative sequences against a reference database; create OTU consensus assignment if needed; one parallel job
+ end: create OTU tables, extract representative sequences and count reads; combine each logs from all steps; add sample's metadata from ENA BioProject to BIOM file; add functional annotation; one serial job
+ archiving: compress raw demultiplexed reads (if demultiplexing), pipeline outputs and processing file in separated tar.gz archives; copy to the output directory
+ doc: create the complete documentation of the pipeline processes; one serial job

Each step job is waiting in the queue so long its preceeding step job is not completed, except for the cut_db step which do not have to wait on any job completion. For checkpointing, the init step job will wait on the completion of the last common step job from the previous SUBPROJECT.


## Outputs

### documentation:

SUBPROJECT.documentation.txt: The file describes the whole processing of the libraries through the pipeline in a human readable format. This file is not archived and will be directly outputted in the SUBPROJECT output directory.

### files in SUBPROJECT.outputs.tar.gz at the end of a full analyses :

+ TAB-separated OTU matrices ("_OTUS.tsv") and representative sequences ("_repseq.fasta") in four different flavors:
   * “SUBPROJECT.all” for all OTUs
   * “SUBPROJECT.abundant” dominant OTUs as defined by the parameters “Minimum number of reads for dominant OTUs” and "Minimum number of sample for dominant OTUs"
   * “SUBPROJECT.all_TARGET” for all OTUs assigned to the "TARGET" target group as defined by the parameter “Target organisms”
   * “SUBPROJECT.abundant_TARGET” for the dominant OTUs assigned to the "TARGET" target group.

The TAB-separated OTU matrices each contains a dense matrix of read counts per OTUs in each sample, each row corresponding to an OTU and each column to a sample, the second to last column containing the consensus taxonomic assignment (labeled as “taxonomy”), and the last column listing the representative sequence identifiers (labeled as “repseq”). Taxonomic rank with no consensus assignment are labeled with “unidentified”.

+ SUBPROJECT.json.biom: biom sparse OTU matrix for all OTUs, JSON encoded, with BioProject metadata annotation for samples and taxonomy, assignment bootstraps and function annotations for observations
+ SUBPROJECT.read_counts.tsv: tab separated table of read counts in each sample for each step of the pipeline
+ SUBPROJECT.log: contains the compilation of all logs (standard output and standard error normally printed in the terminal) outputted by the different steps in their called order
+ SUBPROJECT.raw_and_pair-end_reads_statistics.pdf: contains graphical representations of length and quality distributions as well as the average quality over the raw sequences for each sequence library, for both the full library and for each sample in the library.
+ configuration.SUBPROJECT.tsv: the original configuration file
+ demultiplexing_check.csv: Informations on demultiplexing efficiency and potential barcode mismatch for each run/lane library (only for originally non-demultiplexed libraries)


## Troubleshooting

Job execution could be visualized with the grid engine native tool `qstat` (grid Engine) or `squeue` (SLURM).

If a job is waiting in error state or if it dependencies is never satisfied, it means that the job or its preceding job exit due to an error detected inside of the job.

If it concerns the quality step, the most common error is a number of reads below the set limit.

For any other jobs, it means that the previous step unexpectedly terminate before its end or that an error was issued by mothur in the previous step.

In all those cases, check the standard output (.out) and standard error (.err) log files of the respective step(s) and array(s), which are situated in the "Path to execution location"/SUBPROJECT/log directory.

To detect job terminated due to overpassing requested memory and/or time, compare the requested memory/time in the problematic step's script with the maximum used memory/runtime during job execution. To print job's record after execution, use `qacct` (grid Engine) or `sacct` (SLURM).

For unsolved issues, send an email to guillaume.lentendu@unine.ch, including the .out and .err log files as well as the output of the qacct command for the problematic job.


## References
+ Bengtsson‐Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A., Wit, P., Sánchez‐García, M., Ebersberger, I., Sousa, F., Amend, A., Jumpponen, A., Unterseher, M., Kristiansson, E., Abarenkov, K., Bertrand, Y. J. K., Sanli, K., Eriksson, K. M., Vik, U., Veldre, V., Nilsson, R. H., Bunce, M., 2013. Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods in Ecology and Evolution 4, 914–919. doi:[10.1111/2041-210X.12073](http://doi.org/10.1111/2041-210X.12073)
+ Callahan, B.J., McMurdie, P.J., Rosen, M.J., Han, A.W., Johnson, A.J.A., Holmes, S.P., 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods 13, 581–583. doi:[10.1038/nmeth.3869](http://doi.org/10.1038/nmeth.3869)
+ Edgar, R.C., 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26, 2460–2461. doi:[10.1093/bioinformatics/btq461](http://doi.org/10.1093/bioinformatics/btq461)
+ Fu, L., Niu, B., Zhu, Z., Wu, S., Li, W., 2012. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics 28, 3150–3152. doi:[10.1093/bioinformatics/bts565](http://doi.org/10.1093/bioinformatics/bts565)
+ Gaspar, J.M., 2018. NGmerge: merging paired-end reads via novel empirically-derived models of sequencing errors. BMC Bioinformatics 19, 536. doi:[10.1186/s12859-018-2579-2](http://doi.org/10.1186/s12859-018-2579-2)
+ Gaspar, J.M., Thomas, W.K., 2015. FlowClus: efficiently filtering and denoising pyrosequenced amplicons. BMC Bioinformatics 16. doi:[10.1186/s12859-015-0532-1](http://doi.org/10.1186/s12859-015-0532-1)
+ Mahé, F., Rognes, T., Quince, C., de Vargas, C., Dunthorn, M., 2015. Swarm v2: highly-scalable and high-resolution amplicon clustering. PeerJ 3, e1420. doi:[10.7717/peerj.1420](http://doi.org/10.7717/peerj.1420)
+ Martin, M., 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal 17, 10–12. doi:[10.14806/ej.17.1.200](http://doi.org/10.14806/ej.17.1.200)
+ Masella, A.P., Bartram, A.K., Truszkowski, J.M., Brown, D.G., Neufeld, J.D., 2012. PANDAseq: paired-end assembler for illumina sequences. BMC Bioinformatics 13, 31. doi:[10.1186/1471-2105-13-31](http://doi.org/10.1186/1471-2105-13-31)
+ McDonald, D., Clemente, J.C., Kuczynski, J., Rideout, J.R., Stombaugh, J., Wendel, D., Wilke, A., Huse, S., Hufnagle, J., Meyer, F., Knight, R., Caporaso, J.G., 2012. The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. GigaScience 1. doi:[10.1186/2047-217X-1-7](http://doi.org/10.1186/2047-217X-1-7)
+ Mercier, C., Boyer, F., Bonin, A., Coissac, É., 2013. SUMATRA and SUMACLUST: fast and exact comparison and clustering of sequences.
+ Nguyen, N.H., Song, Z., Bates, S.T., Branco, S., Tedersoo, L., Menke, J., Schilling, J.S., Kennedy, P.G., 2016. FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild. Fungal Ecology 20, 241–248. doi:[10.1016/j.funeco.2015.06.006](http://doi.org/10.1016/j.funeco.2015.06.006)
+ Rognes, T., Flouri, T., Nichols, B., Quince, C., Mahé, F., 2016. VSEARCH: a versatile open source tool for metagenomics. PeerJ. doi:[10.7717/peerj.2584](http://doi.org/10.7717/peerj.2584)
+ Shen, W., Le, S., Li, Y., Hu, F., 2016. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11, e0163962. [doi:10.1371/journal.pone.0163962](http://doi.org/10.1371/journal.pone.0163962)
+ Schloss, P.D., Westcott, S.L., Ryabin, T., Hall, J.R., Hartmann, M., Hollister, E.B., Lesniewski, R.A., Oakley, B.B., Parks, D.H., Robinson, C.J., Sahl, J.W., Stres, B., Thallinger, G.G., Van Horn, D.J., Weber, C.F., 2009. Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. Applied and Environmental Microbiology 75, 7537–7541. doi:[10.1128/AEM.01541-09](http://doi.org/10.1128/AEM.01541-09)
+ Wang, Q., Garrity, G., Tiedje, J., Cole, J., 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Applied and Environmental Microbiology 73, 5261–5267. doi:[10.1128/AEM.00062-07](http://doi.org/10.1128/AEM.00062-07)
