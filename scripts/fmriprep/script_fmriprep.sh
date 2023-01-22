#!/bin/bash
#Template provided by Daniel Levitas of Indiana University
#Edits by Andrew Jahn, University of Michigan, 07.22.2020
#Edits TErnst 2020-12-04

#User inputs:
bids_root_dir=$1
#/media/diskEvaluation/Evaluation/sfb1280a05study2b/rawdata
bids_derv_dir=$2
#/media/diskEvaluation/Evaluation/sfb1280a05study2b/derivatives/fmriprep
fmriprep_work_dir=$3
subj=$4
nthreads=12
mem=60 #gb
container=docker #docker or singularity

#Begin:

#Convert virtual memory from gb to mb
mem=`echo "${mem//[!0-9]/}"` #remove gb at end
mem_mb=`echo $(((mem*1000)-5000))` #reduce some memory for buffer space during pre-processing

#export TEMPLATEFLOW_HOME=$HOME/.cache/templateflow
export FS_LICENSE=/opt/freesurfer/license.txt
fs_license=/opt/freesurfer/license.txt

#Run fmriprep
if [ $container == singularity ]; then
  unset PYTHONPATH; singularity run -B $HOME/.cache/templateflow:/opt/templateflow $HOME/fmriprep.simg \
    $bids_root_dir $bids_root_dir/derivatives \
    participant \
    --participant-label $subj \
    --skip-bids-validation \
    --md-only-boilerplate \
    --fs-license-file $HOME/Desktop/Flanker/derivatives/license.txt \
    --fs-no-reconall \
    --output-spaces MNI152NLin2009cAsym:res-2 \
    --nthreads $nthreads \
    --stop-on-first-crash \
    --mem_mb $mem_mb \
    -w $HOME
else
  fmriprep-docker $bids_root_dir $bids_derv_dir \
    participant \
    --participant-label $subj \
    --skip-bids-validation \
    --md-only-boilerplate \
    --bold2t1w-dof 12 \
    --fs-license-file $fs_license \
    --fs-no-reconall \
    --output-spaces MNI152NLin2009cAsym:res-2 \
    --nthreads $nthreads \
    --stop-on-first-crash \
    --mem_mb $mem_mb \
    -w $fmriprep_work_dir
  #  -w $HOME  
  # -w $fmriprep_work_dir
  ## known errors with clean workdir
  # --clean-workdir \
  #docker system prune
fi
