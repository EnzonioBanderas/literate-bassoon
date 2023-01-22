#!/bin/bash
#bids_root_dir=/media/diskEvaluation/Evaluation/sfb1280a05study7/rawdata
#bids_derv_dir=/media/diskEvaluation/Evaluation/sfb1280a05study7/derivatives
bids_root_dir=//media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-08/bids_mods
bids_derv_dir=/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-08/derivatives
fmriprep_work_dir=/media/dataVol0/Evaluation/temporary_files/sfb1280a05study7_pilots


if [ ! -d "$fmriprep_work_dir" ]; then
    echo creating dir $fmriprep_work_dir
    mkdir -p "$fmriprep_work_dir"
fi
if [ ! -d "$bids_derv_dir" ]; then
    echo creating dir $bids_derv_dir
    mkdir -p "$bids_derv_dir"
fi
if [ ! -d "$bids_derv_dir/fmriprep" ]; then
    echo creating dir $bids_derv_dir/fmriprep
    mkdir -p "$bids_derv_dir/fmriprep"
fi


echo
echo =================================================
echo Using work directory "$fmriprep_work_dir"
echo =================================================
echo 

mydate=`date +%s`
count=0
nsub=0
for pID in $bids_root_dir/sub-*/
do
    nsub=$((nsub+1))
done

for pID in $bids_root_dir/sub-*/
do
    count=$((count+1))
    pID=${pID%*/}
    pID=${pID##*sub-}
    fn_log="${bids_derv_dir}/fmriprep/sub-${pID}.html"
    fn_wip="${bids_derv_dir}/fmriprep/sub-${pID}_wip.txt"
    #echo $fn_log
    if [ -f "$fn_log" ]; then	
        echo Skipping $pID "("$count of $nsub")", participant already run
    elif [ -f "$fn_wip" ]; then	
        echo Skipping $pID "("$count of $nsub")", participant marked as work-in-progress
    else
        echo "work in progress" >> $fn_wip
        echo `date -Iseconds` : $pID "("$count of $nsub")" >> $fn_wip
        echo `date -Iseconds` : $pID "("$count of $nsub")"
        mydate2=`date +%s`
        echo $((mydate2-mydate)) seconds after script start
        ./script_fmriprep.sh $bids_root_dir $bids_derv_dir $fmriprep_work_dir $pID
        rm $fn_wip
    fi
done
