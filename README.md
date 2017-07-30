# PARCEL

A computational pipeline for analyzing sequencing reads generated from PARCEL experiment to identify genomic regions with 
RNA strutual changes in transcripts.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

## Prerequisites

### Operating Systems

#### Supported Unix distributions

- Ubuntu
- CentOS
- Red Hat Enterprise Linux (please use the CentOS packages and instructions)

#### Job scheduler

- Univa Grid Engine
- TORQUE Resource Manager

#### Tools or packages
- perl >= 5.10
- python >= 3.5.1 (for snakemake)
- R >= 3.1.0
- [GNU parallel](https://www.gnu.org/software/parallel/) >= 20150222
- [GNU sort](https://www.gnu.org/software/coreutils/coreutils.html) >= (GNU coreutils) 8.23
- [pigz](https://zlib.net/pigz/) >= 2.3.1
- [mawk](http://invisible-island.net/mawk/) >= 1.3.4
- [bedtools](https://github.com/arq5x/bedtools2) >= 2.25.0
- [snakemake](http://snakemake.readthedocs.io/en/stable/index.html) >= 3.12.0
- [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html) >= 1.8.1
- [samtools](http://www.htslib.org/doc/samtools.html) >= 1.3.1
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) >= 2.2.4

#### Libraries or modules

##### Perl
- IO::File
- IO::Handle
- List::Util
- Math::Random

can be installed by following commands:

```
perl -MCPAN -e "install App::Cpan"
cpan -i IO::Handle IO::File Math::Random List::Util
```

##### R
- [argparse](https://cran.r-project.org/web/packages/argparse/index.html)
- [adagio](https://cran.r-project.org/package=adagio)
- [data.table](https://cran.r-project.org/web/packages/data.table/) >= 1.10.0
- [edgeR](https://bioconductor.org/packages/edgeR)
- [bedr](https://cran.r-project.org/web/packages/bedr) >= 1.0.2

can be installed by following commands in R:

```{r message = FALSE}
install.packages(c("argparse","adagio","bedr","data.table"));
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
```
## Installing

### Install snakemake

Install snakemake into a virtual environment

```
git clone https://bitbucket.org/snakemake/snakemake.git
cd snakemake
virtualenv -p python3 snakemake
source snakemake/bin/activate
python setup.py install
```
### Download scripts and configuration files from github and add directory of scripts into PATH variable 

```
git clone https://github.com/shenyang1981/PARCEL.git
cd PARCEL/; export PARCELSCRIPTS="${PWD}/scripts"; export PATH="${PARCELSCRIPTS}:$PATH"
```
You may consider put 'export PATH=${PARCELSCRIPTS}:$PATH' <replace PARCELSCRIPTS with real path to PARCEL scripts> into your .bashrc file.
  
### Prepare transcriptome and annotation file

Transcriptome file is in FASTA format and is indexed for Bowtie2.

* transcriptome.fas -- transcriptome file
* transcriptome.size -- Length of each transcript in format: transcriptID{tab}Length
* cdsinfo.txt -- The start and end position of CDS in transcript: transcriptID{tab}start{tab}end{tab}Length

put all files into a folder, "database/C.albican/" for example. Build bowtie2 index with transcriptome file.

```
cd database/C.albican/
bowtie2-build transcriptome.fas transcriptome
```
### Prepare input files and sample information

* sampleList.txt -- information of each sequenced library, including library ID(**LibID**), condition or treatment(**Condition**), replicates(**Replicates**), sequencing batch(**SeqBatch**), experimental batch(**ExperiementalBatch**), comparison batch(**ComparisonBatch**). Samples belonged to **same comparison batch** would be selected for pairwised comparison. 

The format of sampleList.txt is like:

|Species|LibID|Condition|Replicates|SeqBatch|ExperiementalBatch|ComparisonBatch|
|-------|:-----:|:---------:|:----------:|:--------:|:------------------:|:---------------:|
|Candida|V1_1 |control  |rep1      |seq1       |1                 |batch1         |
|Candida|V1_2 |control  |rep2      |seq1       |1                 |batch1         |
|Candida|V1_met_1|met  |rep1      |seq1       |1                 |batch1         |
|Candida|V1_met_2|met  |rep2      |seq1       |1                 |batch1         |

** Note: LibID should be unique as the corresponding sequence should be named as {LibID}.fastq.gz.

* input reads files -- Reads are single-end reads. Name of each file should be {LibID}.fastq.gz (LibID should be as same as in sampleList.txt). All of reads files from same sequencing batch should be put into one folder named by {SeqBatch} as indiciated in sampleList.txt. For example, reads files "V1_1.fastq.gz", "V1_2.fastq.gz", "V1_met_1.fastq.gz" and "V1_met_2.fastq.gz" can be put into folder "input/seq1/"

```
ls input/*
input/sampleList.txt

input/seq1:
V1_1.fastq.gz  V1_2.fastq.gz  V1_met_1.fastq.gz  V1_met_2.fastq.gz
```

### generate config file 

To generate a configuration file for snakemake, several variables need to be defined:
- PARCELSCRIPTS: path to scripts used in pipeline
- PARCELDB: path to folder where transcriptome files are
- PARCELREADSROOT: path to root folder of sequenced reads
- PARCELSAMPLEINFO: path to sampleList.txt file
- PARCELRESULTROOT: path to root folder of results
- PARCELBATCH: batchID indicating which libraries should be selected
- PARCELCONTROL: which condition should be used as control
Configuration file can be generated using script **generateConfigureFile.sh**

```
PARCELDB=database/C.albican/ PARCELREADSROOT=input/ PARCELSAMPLEINFO=input/sampleList.txt PARCELRESULTROOT=result/ PARCELBATCH=batch1 PARCELCONTROL=control generateConfigureFile.sh pipeline/config/conf.template.json > pipeline/config/conf.batch1.json
```
**conf.batch1.json**

## Running Pipeline

The pipeline can be simply run in local mode with configuration file.

```
source {$pathtosnakemake}/snakemake/bin/activate
snakemake -s pipeline/parcek.sk --configfile conf.batch1.json -j 32
```

Or run it by submitting to job scheduler

```
runsnake.sh pipeline/parcek.sk conf.batch1.json testjob 24 24
```

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Miao Sun** 
* **Yang Shen**

## Contact

Please contact us if you find bugs, have suggestions, need help etc. You can either use our mailing list or send us an email:

* [Yang Shen](mailto:sheny@gis.a-star.edu.sg)
* [Niranjan Nagarajan](mailto:nagarajann@gis.a-star.edu.sg)

PARCEL is developed in the Genome Institute of Singapore

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
