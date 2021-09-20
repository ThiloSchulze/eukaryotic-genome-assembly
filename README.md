Eukaryotic Genome Assembly Pipeline
-----------------------------------

> **Note:** This pipeline is currently under construction and contain unfinished parts 🚧.

This is a pipeline for performing _de novo_ whole genome assembly on short read sequences
of eukaryotic origin. Here is a brief summary of what is being executed:

1. Raw reads quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Trimmed reads quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. _De novo_ whole genome assembly ([SPAdes](https://cab.spbu.ru/software/spades/))
5. Assembly quality assessment ([QUAST](http://bioinf.spbau.ru/quast))

### Installing Nextflow

> **Note:** If you already have Nextflow installed, then you may skip to the
> next section.

The workflow manager [Nextflow](https://www.nextflow.io/) is required to run
this pipeline. It requires Bash 3.2 (or later) and Java 8 (or later, up to 15)
to be installed. To download and install Nextflow into the current directory,
run the following in your terminal.

```bash
$ curl -s https://get.nextflow.io | bash
```

Then, make the `nextflow` binary executable by running `chmod +x`:

```bash
$ chmod +x nextflow
```

### Installing Anaconda

> **Note:** If you already have Anaconda installed, then you may skip to the
> next section.

Download Anaconda [from here](https://www.anaconda.com/products/individual)
(I recommend the command-line installer) and run the installer.

### Downloading this Pipeline

To download this pipeline, go to the folder where you wish to store this
pipeline in and then run the following:

```
$ git clone https://github.com/Animal-Evolution-and-Biodiversity/trimming
```

### Formatting the Input Data

The input data may be single-end or paired-end reads. Paired-end reads is the
default. These reads should be in the FASTQ format, either compressed or
uncompressed. The resulting reads are always in a compressed format (`.gzip`).

If your reads are single-ended, use the `--single-end` option. If your reads
are paired-ended, then you don't need to supply any additional options.

### Running this Pipeline

Go to the folder where you downloaded this workflow and then use
`cd eukaryotic-genome-assembly` to enter the base directory. From there you
launch the pipeline like so:

```bash
$ nextflow run . --help
```

### Cleaning Up

Nextflow produces a lot of additional files which takes up space on your drive.
These files are useful when running and troubleshooting the pipeline but may
then safely be removed in our case, since everything you need is saved to the
output directory. For this, we provide a GNU Makefile. Go to the workflow's
base directory and enter the following command to remove unnecessary output
files after a successful run:

```bash
$ make clean
```

© Animal Evolution and Biodiversity 2021
