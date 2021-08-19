#!/usr/bin/env nextflow

nextflow.enable.dsl=2

rawreads = Channel
  .fromFilePairs( params.input, size: params.single_end ? 1 : 2, type: 'file' )
  .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
  .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.input}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'\n" +
             "For single-end reads, specify '--single-end'" }

def helpMessage() {
log.info"""
  Start: 
  nextflow run . --input '/full/path/*{R1,2}.fastq'

  Options:
  --trimming_cutoff [Number]     Define Quality trimming threshold. Default: 20

  """.stripIndent()
}

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

// Display the version number upon request
if ( params.version ) exit 0, versionNumber()

// Input validation
if ( params.input == null ) {
  exit 1, "Missing mandatory argument '--input'\n" +
          "Launch this workflow with '--help' for more info"
}

process quality_control {
  publishDir "${params.output}/quality_control_pre-trimming", mode: 'copy'

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

process trimming {
publishDir "${params.output}/trimmed_reads", mode: 'copy'

  input:
  tuple val(name), path(read)

  output:
  tuple val(name), path("*val_{1,2}.fq.gz"), emit: trimmedReads
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
  if ( ! params.single_end )
    flagsTrimming += " --paired --retain_unpaired"
  commandTrimming = "trim_galore $flagsTrimming $read"

  """
  $commandTrimming
  echo "$commandTrimming" > 'trim_galore_command.txt'
  """
} 


process assembly {
publishDir "${params.output}/assembly", mode: 'copy'

  input:
  tuple val(name), path(trimmed_read)

  output:
  path "${name}_assembly"

  script:
  """
  spades.py -1 ${trimmed_read[0]} -2 ${trimmed_read[1]} -k 51 --careful --cov-cutoff auto -o ${name}_assembly
  """
}

workflow {
  quality_control(rawreads)
  trimming(rawreads)
  assembly(trimming.out.trimmedReads) | view
}