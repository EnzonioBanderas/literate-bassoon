fn_rd1 = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-08/raw_mp2rage/MID00123/meas_MID00123_FID06862_anat_mp2rage_ses_1.dat';
fn_rd2 = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-08/raw_mp2rage/MID00169/meas_MID00169_FID06903_anat_mp2rage_ses_1.dat';
fn_rd3 = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-08/raw_mp2rage/MID00202/meas_MID00202_FID06936_anat_mp2rage_ses_1.dat';
fn_rd4 = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-08/raw_mp2rage/MID00234/meas_MID00234_FID06968_anat_mp2rage_ses_1.dat';

addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/spm12
addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/mapVBVD_20150918withFatNavs
addpath(genpath('/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/retroMoCoBox/'));

%reconstructSiemensMP2RAGEwithFatNavs(fp_rd1)
reconstructSiemensMP2RAGEwithFatNavs(fp_rd2)
reconstructSiemensMP2RAGEwithFatNavs(fp_rd3)
reconstructSiemensMP2RAGEwithFatNavs(fp_rd4)