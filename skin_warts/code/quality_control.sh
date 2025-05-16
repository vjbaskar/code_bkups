#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# title: Build QC jobs.
# purpose: Writing a simple bash script to run papermill with the QC template.
#
# created: 2025-02-05 Wed 13:27:02 GMT
# updated: 2025-02-05
# version: 0.0.9
# status: Prototype
#
# maintainer: Ciro Ramírez-Suástegui
# author:
#   - name: Ciro Ramírez-Suástegui
#     affiliation: The Wellcome Sanger Institute
#     email: cs59@sanger.ac.uk, cramsuig@gmail.com
# ------------------------------------------------------------------------------
#% SYNOPSIS
#+    ${SCRIPT_NAME} [-hv] [-t[--template]] [-r[--run_path]] [-o[--output]] ...
#%
#+ FNM=code/quality_control.sh
#+ bash ${FNM}.ext 2>&1 | tee -a .logs/${FNM/\//.}_$(date +%Y%m%d-%H%M%S)
# END_OF_HEADER

## Environment setup ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Quick installations --------------------------------------
# Basic packages -------------------------------------------
# Logging configuration ------------------------------------
set -e
trap cleanup SIGINT SIGTERM ERR EXIT
# In-house/developing --------------------------------------
#source code/misc.sh
source code/utils.sh

SCRIPT_HEADSIZE=$(head -200 ${0} |grep -n "^# END_OF_HEADER" | cut -f1 -d:)
SCRIPT_NAME=$(basename ${0})

function usage() {
  printf "Usage: "
  head -${SCRIPT_HEADSIZE:-99} $(realpath ${0}) | grep -e "^#+" \
  | sed -e "s|^#+[ ]*||g" -e "s|\${SCRIPT_NAME}|${SCRIPT_NAME}|g"
}
function usagefull() {
  head -${SCRIPT_HEADSIZE:-99} $(realpath ${0}) | grep -e "^#[%+-]" \
  | sed -e "s|^#[%+-]||g" -e "s|\${SCRIPT_NAME}|${SCRIPT_NAME}|g"
}
function scriptinfo() {
  head -${SCRIPT_HEADSIZE:-99} $(realpath ${0}) | grep -e "^#-" \
  | sed -e "s|^#-||g" -e "s|\${SCRIPT_NAME}|${SCRIPT_NAME}|g"
}

function cleanup() {
  trap - SIGINT SIGTERM ERR EXIT
}

function logger() {
    LOGLINE="[$(date '+%Y-%m-%d %T')] \033[1;32mINFO\033[0m "
    LOGLINE="${LOGLINE}[${0##*/}:${FUNCNAME}:${BASH_LINENO}]"
    echo -e "${LOGLINE} :: $*"  >&2
}

# Tool (packaged) modules ----------------------------------

## Global configuration ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while :; do
  case "${1-}" in
    -h | --help) usage; exit ;; # Display a usage synopsis.
    -v | --verbose) set -x ;;
    --no-color) NO_COLOR=1 ;;
    -t | --template) TEMPLATE="${2-}"; shift;;
    -r | --run_path) RUN_PATH="${2-}"; shift;;
    -o | --output) OUTPUTS="${2-}"; shift;;
    -?*) echo "Unknown option: $1" ;;
    *) break ;;
  esac
  shift
done


# Scripts folders
OUTCODE=$(echo ${TEMPLATE%.*} | sed 's/analysis/code/')
ROUTINE=$(basename ${TEMPLATE%.*})

logger "Present parameters:"
echo "Working at: '$(pwd)'"
echo "Template: ${TEMPLATE}"
echo "Path: ${RUN_PATH}"
echo "Step: ${ROUTINE}"
echo "Scripts: ${OUTCODE}"
echo "Output: ${OUTPUTS}"

## Loading data ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNAMES=($(ls -d ${RUN_PATH}/*))
STRING=($(find_substring_names ${SNAMES[@]}))

## Pre-processing ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir -p ${OUTPUTS}
mkdir -p ${OUTCODE}

## Main ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for SNAME in ${SNAMES[@]}; do
  # Extracting sample name
  NNAME=$(echo ${SNAME} | sed 's|'"${STRING[0]}"'||g; s|'"${STRING[1]}"'||g')
  EXECF=${OUTCODE}/${NNAME}.sh
  LOGGF=${WORKDIR}/../.logs/${ROUTINE}_${NNAME}
  logger " \033[0;32m${NNAME}\033[0m"
  if [[ -f ${LOGGF}*.err ]]; then
    rm ${LOGGF}* # remove previous logs
    rm ${OUTPUTS}/${NNAME}.${TEMPLATE##*.} # remove previous QC
  fi
  # Creating bash script to run the QC # -------------------
  echo "#!/usr/bin/env bash" > ${EXECF}
  echo "" >> ${EXECF}
#  echo "if [[ ~/.conda/init.sh ]]; then" >> ${EXECF}
#  echo "  . ~/.conda/init.sh" >> ${EXECF}
#  echo -e "fi\n" >> ${EXECF}
#  echo -e "mamba activate squidpy_v1.4.1\n" >> ${EXECF}
  echo "NBTMP=${TEMPLATE}" >> ${EXECF}
  echo "NBOUT=${OUTPUTS}/${NNAME}.${TEMPLATE##*.}" >> ${EXECF}
  echo -e "SPATH=${SNAME}\n" >> ${EXECF}
  printf "papermill \${NBTMP} \${NBOUT} " >> ${EXECF}
  echo -n "-p sample_path \"\${SPATH}\" " >> ${EXECF}
  echo -n "-p indata_name \"${NNAME}\"" >> ${EXECF}
  chmod +x ${EXECF}
  #bsub -G team298 \
  #  -o ${LOGGF}_$(date '+%Y-%m-%d').out \
  #  -e ${LOGGF}_$(date '+%Y-%m-%d').err \
  #  -Is -n1 -q normal \
  #  -R "select[mem>6000] rusage[mem=6000]" -M6000 \
   # -W7:30 -J "${ROUTINE}_${NNAME}" ${EXECF} &
  echo "${EXECF}"
  ${EXECF}
done

## Conclusions ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Save ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
