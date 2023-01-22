#!/bin/bash
bids_root_dir=/media/diskEvaluation/Evaluation/sfb1280a05study2b/rawdata
bids_derv_dir=/media/diskEvaluation/Evaluation/sfb1280a05study2b/derivative
fmriprep_work_dir=/media/dataVol0/Evaluation/temporary_files/sfb1280a05study2b

if [ ! -d "$fmriprep_work_dir" ]; then
        echo creating dir $fmriprep_work_dir
	mkdir -p "$fmriprep_work_dir"
fi

echo Using work direktory "$fmriprep_work_dir"

for pID in $bids_root_dir/sub-*/
do
        pID=${pID%*/}
        pID=${pID##*sub-}
        if [ -f "$bids_derv_dir"+"/fmriprep/sub-"+$pID+".html" ]; then	
		echo Skipping $pID, participant already run
	else
		echo $pID
		./script_fmriprep.sh $bids_root_dir $bids_derv_dir $fmriprep_work_dir $pID
	fi
done