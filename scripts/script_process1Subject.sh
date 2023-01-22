#!/bin/bash


export ANTSPATH=/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/ANTs/install/bin/
export SKSCRIPTS=/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/sk_ants_scripts/
export PATH=${SKSCRIPTS}:${ANTSPATH}:$PATH

pID=$1
fpath="/media/dataVol0/Evaluation/temporary_files/sfb1280a05study7/"$pID"/"
tr=1.99


echo $pID
fn_t1map=$fpath$pID"_ses-01_calcT1map_clean.nii.gz"

#sub-SE05NS19_ses-01_task-fear_run-1_bold_brain_DistCorr_template0
for fn_func_tmpl in $(find $fpath -regextype posix-extended -regex '^.*task-fear.*bold_brain_DistCorr_template0\.nii\.gz')
do
  fn_func=${fn_func_tmpl/_DistCorr_template0}
  ses=`expr "$fn_func" : '.*\(ses-[0-9]\+\)'`
  run=`expr "$fn_func" : '.*\(run-[0-9]\)'`
  prefix=$fpath$pID"_reg-t1map2func_"$ses"_"$run"_"
  fn_coreg=$prefix"0GenericAffine.mat"
  echo
  echo template: $fn_func_tmpl
  echo func:     $fn_func
  echo T1map:    $fn_t1map
  echo prefix:   $prefix
  echo coreg:    $fn_coreg
  echo  
  
  fn_final="${fn_func/.nii.gz/_MoCorr_DistCorr.nii.gz}"
  if [ -f $fn_final ]; then
    echo already finished, exists: $fn_final
  else
    antsRegistrationSyN.sh \
    -e 42            `# Fix random seed to an int value` \
    -d 3             `# ImageDimension: 2 or 3, here 3 obviously` \
    -f $fn_func_tmpl `# Fixed image(s) or source image(s) or reference image(s)` \
    -m $fn_t1map     `# Moving image(s) or target image(s)`\
    -o $prefix       `# OutputPrefix: A prefix that is prepended to all output files.` \
    -t r             `# transform type (default = 's'), r=rigid(1 stage)`\
    -p f             `# precision type (default = 'd'), f: float, d: double` \
    -j 1             `# use histogram matching, 0=false or 1=true (default = 0)` \
    -n 12            `# Number of threads`
  
    sk_ants_Realign_Reslice_editTE.sh -t $tr -n 12 -a $fn_func -y $fn_coreg
  fi

done
