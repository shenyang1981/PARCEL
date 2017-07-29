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
- python >= 3.15 (for snakemake)
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

### download scripts and configuration files from github

```
git clone https://github.com/shenyang1981/PARCEL.git
```
### Prepare transcriptome and annotation file

Transcriptome file is in FASTA format and is indexed for Bowtie2.

* transcriptome.fas -- transcriptome file
* transcriptome.size -- Length of each transcript in format: transcriptID{tab}Length
* cdsinfo.txt -- The start and end position of CDS in transcript: transcriptID{tab}start{tab}end{tab}Length

put all files into database
```
cd database/C.albican/
bowtie2-build transcriptome.fas transcriptome
```
### Prepare input files and sample information

* sampleList.txt -- information of each sequenced library, including library ID(**LibID**), condition or treatment(**Condition**), replicates(**Replicates**), sequencing batch(**SeqBatch**), experimental batch(**ExperiementalBatch**), comparison batch(**ComparisonBatch**). Samples belonged to **same comparison batch** would be selected for pairwised comparison. 

The format of sampleList.txt is like:

|Species|LibID|Condition|Replicates|SeqBatch|ExperiementalBatch|ComparisonBatch|
|-------|:-----:|:---------:|:----------:|:--------:|:------------------:|:---------------:|
|Candida|V1_1 |control  |rep1      |1       |1                 |batch1         |
|Candida|V1_2 |control  |rep2      |1       |1                 |batch1         |
|Candida|V1_met_1|met  |rep1      |1       |1                 |batch1         |
|Candida|V1_met_2|met  |rep2      |1       |1                 |batch1         |

** Note: LibID should be unique as the corresponding sequence should be named as {LibID}.fastq.gz.

### generate config file 

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
