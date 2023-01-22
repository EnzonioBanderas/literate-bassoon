#!/bin/bash

fp_source=/media/dataVol0/Evaluation/temporary_files/sfb1280a05study2b

for pdir in $(ls -d $fp_source/sub-*) 
do
  #echo $pid
  pid=$(basename -- "$pdir")
  ./script_process1Subject.sh $pid
done
