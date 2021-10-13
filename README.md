Eukaryotic Genome Assembly Pipeline
-----------------------------------

This pipeline performs raw- and trimmed read quality control, adapter trimming,
assembly, and assembly quality assessment on short read sequence data of
eukaryotic origin. Here is a brief list of what is being executed:

1. Raw reads quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Trimmed reads quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. _De novo_ whole genome assembly ([SPAdes](https://cab.spbu.ru/software/spades/))
5. Assembly quality assessment ([QUAST](http://bioinf.spbau.ru/quast))

### Installing Nextflow

[Nextflow](https://www.nextflow.io/) is required to run this pipeline. For
installing Nextflow, simply follow the [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)
in their documentation.

### Downloading this Pipeline

To download this pipeline, go to the folder where you wish to store this
pipeline in and then run the following:

```
$ git clone https://github.com/ThiloSchulze/eukaryotic-genome-assembly
```

### Formatting the Input Data

The input data may be single-end or paired-end reads. Paired-end reads is the
default. If your reads are single-ended, use the `--single-end` option. If your
reads are paired-ended, then you don't need to supply any additional options.
These reads should be in the FASTQ format, either compressed or uncompressed.
The resulting reads are always in a compressed format (`.gzip`).

### Running this Pipeline

Go to the folder where you downloaded this workflow and then use
`cd eukaryotic-genome-assembly` to enter the base directory. From there you
launch the pipeline like so:

```bash
$ nextflow run . --help
```

### Choosing a profile

All programming versioning and installations are handled by the pipeline itself.
We provide the following profiles for managing these installations:

1. `conda` ([Anaconda](https://www.anaconda.com))
2. `docker` ([Docker](https://www.docker.com/))
3. `singularity` ([Singularity](https://sylabs.io/))
4. `cluster` ([Slurm](https://slurm.schedmd.com/documentation.html))

For example, to run this pipeline using Anaconda, launch the pipeline with
the following flag:

```
nextflow run . -profile conda
```

You can provide multiple profiles by separating the profiles with a comma (`,`).
To run the pipeline with the profiles `cluster` and `singularity` enabled,
launch it like so:

```
nextflow run . -profile cluster,singularity
```

### Using a custom configuration

The default configuration is stored in the `nextflow.config`. We encourage you
to replace the default values with whatever is appropriate for your needs. The
easiest way to do this is by creating a custom file. We have provided an
example of such a file in the `configs` directory. To use this configuration,
simply run this pipeline with the following options:

```
nextflow run . -config 'configs/example.config'
```

You could make a copy of this configuration and then make your changes from
there. Then just replace `example.config` with whatever name you chose for your
custom file.

#### In a Linux cluster environment

If you intend to run this pipeline on your local cluster, you may consider using
our wrapper script (`ega_wrapper.sh`) for automatically generating the initial
Nextflow script. This currently only works for paired-end reads.

Your raw reads should be contained within a directory called `raw_reads`. The
program will automatically detect the difference between the forward and reverse
file.

For example, in a directory structure like this:

```
DC-1_Perinereis_nuntia
└── raw_reads
    ├── DC-1_S0_L000_R1_001.fastq.gz
    └── DC-1_S0_L000_R2_001.fastq.gz
```

The following files are created:

```
DC-1_Perinereis_nuntia
├── ega_assembly_out
│   └── batch_job.sh
└── raw_reads
    ├── DC-1_S0_L000_R1_001.fastq.gz
    └── DC-1_S0_L000_R2_001.fastq.gz
```

To launch this job, you run `batch_job.sh`. Note that this file does not launch
the assembly jobs and such, these are controlled by the Nextflow configuration
file and are generated automatically. This is just to generate the batch job
that launches the Nextflow pipeline itself.

Here is an example of what such file might look like:

```
#!/bin/bash
#SBATCH --job-name=DC-1_Per
#SBATCH --partition=medium
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --constraint=scratch2
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --output=/scratch2/thalen/DC-1_Perinereis_nuntia/ega_assembly_out/DC-1_Per_stdout.txt
#SBATCH --error=/scratch2/thalen/DC-1_Perinereis_nuntia/ega_assembly_out/DC-1_Per_stderr.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.waehling@gmx.de

module purge
module load nextflow
module load singularity

export NXF_SINGULARITY_CACHEDIR="/scratch2/thalen/DC-1_Perinereis_nuntia/ega_assembly_out/singularity_cachedir"
mkdir -p "/scratch2/thalen/DC-1_Perinereis_nuntia/ega_assembly_out/singularity_cachedir"

nextflow run "/usr/users/thalen/pipelines/eukaryotic-genome-assembly"\
  -profile cluster,singularity\
  -resume\
  --reads "/scratch2/thalen/DC-1_Perinereis_nuntia/raw_reads/DC-1_S0_L000_R{1,2}_001\.fastq\.gz"
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
