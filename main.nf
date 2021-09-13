#!/usr/bin/env nextflow

nextflow.enable.dsl=2

ch_rawReads = Channel
  .fromFilePairs( params.reads, size: params.single_end ? 1 : 2, type: 'file' )
  .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
  .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.reads}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'\n" +
             "For single-end reads, specify '--single-end'" }

def helpMessage() {
log.info"""
  Usage:
    (1) Single-end reads:
      nextflow run . --single-end --reads '*\\.fastq\\.gz'
    (2) Paired-end reads:
      nextflow run . --reads '*_R{1,2}\\.fastq\\.gz'

  Description:
    A pipeline for performing de novo whole genome assembly on short read
    sequences of eukaryotic origin.

  Pipeline summary:
    1. Raw reads quality control (FastQC)
    2. Adapter trimming (Trim Galore!)
    3. Trimmed reads quality control (FastQC)
    4. De novo whole genome assembly (SPAdes)
    5. Assembly quality assessment (QUAST)

  Mandatory arguments:
    --reads             path to one or more sets of paired-ended reads (valid
                        file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq')

  Input/output options:
    --single_end        when specified, input reads are single-end reads
                        (default: $params.single_end)
    --output            path to a directory which the results are written to
                        (default: $params.output)

  Quality control:
    --qc_adapters       path to adapters files, if any (default: $params.qc_adapters)

  Trimming (Trim Galore!):
    --trim_min_length   discards reads shorter than this (default: $params.trim_min_length)
    --trim_quality      Phred score threshold for quality trimming (default: $params.trim_quality)
    --trim_adapter      adapter sequence to be trimmed (default: auto-detect)
    --trim_phred64      use Phred+64 (i.e., Illumina 1.5) encoding for quality scores
                        (default: Phred+33; Sanger/Illumina 1.8)
    --trim_forward_leading
                        cut off bases at the start of the forward reads or all
                        reads if single-end reads (default: $params.trim_forward_leading)
    --trim_forward_trailing
                        cut off bases at the end of the forward reads or all
                        reads if single-end reads (default: $params.trim_forward_trailing)
    --trim_reverse_leading
                        cut off bases at the start of the forward reads or all
                        reads if single-end reads (default: $params.trim_reverse_leading)
    --trim_reverse_trailing
                        cut off bases at the end of the forward reads or all
                        reads if single-end reads (default: $params.trim_reverse_trailing)
    --trim_forward_cutoff POSITION
                        remove all bases past this position
                        cut off bases at the end of the forward reads or all
                        reads if single-end reads (default: $params.trim_forward_trailing)
    --trim_leading_cutoff
                        INSTEAD of trimming, remove all bases past this position
                        (default: '$params.trim_leading_cutoff')
    --trim_trailing_cutoff
                        INSTEAD of trimming, remove all bases such that there
                        are this many bases from the 3' end (default: '$params.trim_trailing_cutoff')

  Miscellaneous:
    --help              display this help message and exit
    --version           display this pipeline's version number and exit

  """.stripIndent()
}

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

// Display the version number upon request
if ( params.version ) exit 0, versionNumber()

// Input validation
if ( params.reads == null ) {
  exit 1, "Missing mandatory argument '--reads'\n" +
          "Launch this workflow with '--help' for more info"
}

// // Check whether the provided path exists
// def inputPath = new File(params.reads)
// assert inputPath.exists() : "Provided reads, $params.reads, not found"

/*
 * Read quality control using FastQC
 */
process qualityControl {
  publishDir "${params.output}/qualityControl_pre-trimming", mode: 'copy'

  input:
  tuple val(name), path(control)

  output:
  path "*fastqc.{html,zip}"
  path "fastqc_command.txt"

  script:
  """
  fastqc_command="fastqc --threads ${task.cpus} --quiet $control"
  \$fastqc_command
  echo "\$fastqc_command" > 'fastqc_command.txt'
  rename 's/_fastqc\\.zip\$/_pre-trimming_fastqc.zip/' *_fastqc.zip
  rename 's/_fastqc\\.html\$/_pre-trimming_fastqc.html/' *_fastqc.html
  """
}

/*
 * Adapter trimming and using Trim Galore! Includes quality control of trimmed
 * adapters as well.
 */
process trimming {
  publishDir "${params.output}/trimmed_reads", mode: 'copy'

  input:
  tuple val(name), path(read)

  output:
  tuple val(name), path("*{1,2}.fq.gz"), emit: trimmedReads
  path "*.txt"
  path "*.{zip,html}"
  path "*.fq.gz"

  script:

  flagsTrimming = "--fastqc --gzip --quality $params.trim_quality \
--length $params.trim_min_length --cores $task.cpus"
  if ( params.trim_phred64 )
    flagsTrimming += " --phred64"
  if ( params.trim_forward_leading )
    flagsTrimming += " --clip_R1 $params.trim_forward_leading"
  if ( params.trim_forward_trailing )
    flagsTrimming += " --three_prime_clip_R1 $params.trim_forward_trailing"
  if ( params.trim_reverse_leading )
    flagsTrimming += " --clip_R2 $params.trim_reverse_leading"
  if ( params.trim_reverse_trailing )
    flagsTrimming += " --three_prime_clip_R2 $params.trim_reverse_trailing"
  if ( !params.single_end )
    flagsTrimming += " --paired --retain_unpaired"
  commandTrimming = "trim_galore $flagsTrimming $read"

  """
  $commandTrimming
  echo "$commandTrimming" > 'trim_galore_command.txt'
  """
}

/*
 * De novo assembly using the SPAdes assembler.
 */
process assembly {
  echo true
  publishDir "${params.output}/assembly", mode: 'copy'

  input:
  tuple val(name), path(trimmed_reads)
  val kmers

  output:
  path "${name}_assembly" optional true

  script:

  // Make list of kmers SPAdes-compatible ([a, b, c] -> "a,b,c")
  kmers_formatted = kmers.toString().replaceAll("[ \\[\\]]", "")

  """
  cat "${trimmed_reads[0]}" "${trimmed_reads[2]}" > unpaired_reads.fq.gz
  spades.py \
    -1 ${trimmed_reads[1]} \
    -2 ${trimmed_reads[3]} \
    -s unpaired_reads.fq.gz \
    -k "$kmers_formatted" \
    --careful \
    --cov-cutoff auto \
    -o ${name}_spades_out
  """
}

workflow {
  qualityControl(ch_rawReads)
  trimming(ch_rawReads)
  assembly(trimming.out.trimmedReads, params.kmers) | view
}
