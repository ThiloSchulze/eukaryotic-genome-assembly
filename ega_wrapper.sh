#!/usr/bin/env bash
#
# Generate a job file for the eukaryotic genome assembly (EGA) pipeline

set -o errtrace
set -o errexit
set -o nounset
set -o pipefail

# Show line numbers and function names when debugging
export PS4='+ ${LINENO}:${FUNCNAME[0]:-}() '

readonly ARG_COUNT="$#"
readonly VERSION="0.1.0"
readonly EGA_BIN="/usr/users/${USER}/pipelines/eukaryotic-genome-assembly"
readonly TIMESTAMP=$( date +%Y-%m-%d_%H-%M-%S )
basedir=''
name=''
outdir='ega_assembly_out'
reads='raw_reads'
mail='alexander.waehling@gmx.de'

usage() { printf "%s" "\
usage:
  ega_wrapper.bash [--help] [--version] DIRECTORY

description:
  Generate a batch job for the eukaryotic genome assembly (EGA) pipeline.

options:
  -h, --help       display this help message and exit
  -v, --version    display the version number and exit
  -m, --mail       report change of job status (started, finished, etc.) to this
                   email address (default: 'alexander.waehling@gmx.de')
  -n, --name       give your batch job a custom name (default: first 8 letters
                   of \`DIRECTORY\`')
  -o, --outdir     output directory where (default: '$outdir')
  -r, --reads      name of the subdirectory where raw reads are located
                   (default: '$reads')

expected input directory structure:
  <basedir>
  └── <reads>
      ├── <prefix>1.fq
      └── <prefix>2.fq
"
  exit 1
}

# Display the version number and exit.
version_info() {
  echo "$VERSION"
  exit 1
}

# Display the provided error message before exiting with the provided status
# (default: 1).
error() {
  local message="$1"
  local status="${2-1}" # default exit status: 1
  echo "ega_wrapper: error: $message"
  exit "$status"
}

# Returns the first 8 letters of the basename of the provided directory
get_name() {
  local path="$1"
  basename "$path" |\
    cut -c1-8
}

batch_job() {
  local name="$1"
  local basedir="${2}"
  local raw_reads_pattern="${3}"
  local dir_out="${4}"
  tee "${dir_out}/batch_job.sh" << EOF
#!/bin/bash
#SBATCH --job-name=$name
#SBATCH --partition=medium
##SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --constraint=scratch2
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --output=${dir_out}/${name}_stdout.txt
#SBATCH --error=${dir_out}/${name}_stderr.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$mail

module purge
module load nextflow
module load singularity

readonly SINGULARITY_CACHEDIR="${dir_out}/singularity_cachedir"
export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"
mkdir -p "$SINGULARITY_CACHEDIR"

nextflow run "${EGA_BIN}"\\
  -profile cluster,singularity\\
  -resume\\
  --reads "${raw_reads_pattern}"\\
  --max_memory '256.GB'
EOF
}

# Takes the path to a FASTQ file with forward reads and the path to another
# FASTQ file with reverse reads as an input. Returns a pattern that captures
# both of these reads. Also escapes dot to adhear to the expected Nextflow
# input format.
#
# example usage:
#   $ get_pattern NG-26745_1.fastq.gz NG-26745_2.fastq.gz
#   NG-26745_{1,2}\.fastq\.gz
get_pattern() {
  local file_a="$1"
  local file_b="$2"
  [[ "${#file_a}" -eq "${#file_b}" ]] || error "unequal length of filenames in \
raw reads: $file_a $file_b"
  mismatch_found=false
  before=''
  after=''
  pattern=''
  for (( i=0; i<${#file_a}; i++ ))
  do
    if [[ "${file_a:$i:1}" == "${file_b:$i:1}" ]]
    then
      if $mismatch_found
      then
        after+="${file_a:$i:1}"
      else
        before+="${file_a:$i:1}"
      fi
    else
      if $mismatch_found
      then
        error "found more than one difference in FASTQ files at position: $i"
      fi
      pattern="{${file_a:$i:1},${file_b:$i:1}}"
      mismatch_found=true
    fi
  done
  echo "${before}${pattern}${after}" | sed 's/\./\\\./g'
}

# Returns all FASTQ files found in the provided `path`
fastq_files() {
  local path="$1"
  find "$path" -type f \( -iname \*.fastq.gz -o -iname \*.fq.gz \) |\
    sort
}

# Parse user-provided arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -h | --help)
      usage
      ;;
    -v | --version)
      version_info
      ;;
    -m | --mail)
      mail="$2"
      shift # past argument
      shift # past value
      ;;
    -n | --name)
      name="$2"
      shift # past argument
      shift # past value
      ;;
    -r | --reads)
      reads="$2"
      shift # past argument
      shift # past value
      ;;
    --) # end of all options
      break
      ;;
    -*) # unknown option
      error "unknown option: $key"
      ;;
    *) # end of options
      [[ -n "$basedir" ]] && error "unrecognized positional argument: $key"
      basedir="$key"
      shift # past argument
      ;;
  esac
done

main() {
  # Display an error message if no arguments were provided
  [[ ARG_COUNT -lt 1 ]] && usage
  [[ -z "$basedir" ]] && error "missing mandatory argument DIRECTORY"
  [[ ! -d "$basedir" ]] && error "provided directory not found: $basedir"
  basedir_abs=$( realpath "$basedir" )
  raw_reads_dir="${basedir_abs}/${reads}"
  [[ ! -d "$raw_reads_dir" ]] && error "missing subdirectory 'raw_reads' in \
provided directory: $basedir"
  # store the output of `fastq_files` to an array, `files`
  mapfile -t files < <( fastq_files "$raw_reads_dir" )
  [[ ${#files[@]} -eq 2 ]] || error "more or less than 2 FASTQ files found in: \
$raw_reads_dir"
  raw_reads_pattern=$( get_pattern "${files[0]}" "${files[1]}" )
  # echo "$raw_reads_pattern"
  outdir_full="${basedir_abs}/${outdir}"
  mkdir -p "$outdir_full"
  [[ -z "$name" ]] && name=$( get_name "$basedir_abs" )

  echo -e "$( basename "${0}" )"
  echo -e "developer: Felix Thalen <felix.thalen@uni-goettingen.de>\n"
  echo -e "\e[4mBatch job:\e[0m"
  batch_job "$name" "$basedir_abs" "$raw_reads_pattern" "$outdir_full"

  echo -e "\nwrote batch job file to ${outdir_full}/batch_job.sh"
  echo -e "type \`sbatch ${outdir_full}/batch_job.sh\` to submit the batch job"
  echo -e "type \`scancel <job_id>\` to cancel an incorrect job submission"
}

main
