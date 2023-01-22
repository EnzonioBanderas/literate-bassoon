

addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/spm12
addpath /media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/mapVBVD_20150918withFatNavs
addpath(genpath('/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/retroMoCoBox/'));


%reconstructSiemensMP2RAGEwithFatNavs('/path/filename.dat')

fl = essbids_listFiles(fullfile(fp_s,'**','sub-*.dat'))

% - loop over fl
% - check wether reconstructed data is already there, then skip
% - check wether reconstruction is in process
%         fid = fopen('wip.txt','w');
%         fprintf(fid,'currenttime starting process');
%         fclose(fid);
%          ...
%         delete(fn_wip)
%  - execute the reconstruction script if conditions are met


