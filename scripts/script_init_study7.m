fp_pwd = pwd;
fp_scr = fileparts(mfilename('fullpath'));
if isempty(fp_scr)
  % this happens when you do not run this as a script but step by step from
  % console. in this case use a hardwired path, You may need to modify it
  if ispc
    fp_scr = 'R:\Evaluation\sfb1280a05study7\scripts';
  else
    fp_scr = '/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts';
  end
end
fp0    = fileparts(fp_scr);
fp_d   = fullfile(fp0,'rawdata');
fp_s   = fullfile(fp0,'sourcedata');
fp_us  = fullfile(fp0,'sourcefiles_unsorted');
fp_de  = fullfile(fp0,'derivatives');
fp_dis = fullfile(fp0,'discarded');

dl2chk = {fp0,fp_d,fp_s,fp_us,fp_de,fp_dis};
for i=1:numel(dl2chk)
  if not(isfolder(dl2chk{i}))
    mkdir(dl2chk{i});
  end
end

addpath(fullfile(fp_scr,'packages','pspm_v5.0.0'));
addpath(fullfile(fp_scr,'packages','toolbox_ter'));
addpath(fullfile(fp_scr,'packages','export_fig'));
%addpath(fullfile(fp_scr,'packages','SPM','startspm'));
addpath(fullfile(fp_scr,'packages','startspm'));
addpath(fp_scr)
addpath(fullfile(fp_scr,'data_preparation'));
addpath(fullfile(fp_scr,'data_preparation_study7'));
addpath(fullfile(fp_scr,'packages','essbids'));
addpath(fullfile(fp_scr,'packages','GetFullPath_20210620'));
% addpath(genpath(fullfile(fp_scr,'packages','presurfer')));
addpath(genpath(fullfile(fp_scr,'packages','Version3_UniDue')));  % QSM scripts
addpath(genpath(fullfile(fp_scr,'packages','MP2RAGE-related-scripts')));
% addpath(genpath(fullfile(fp_scr,'packages','load_acq_20110222')));
% addpath(genpath(fullfile(fp_scr,'arrington-eyetracker-scripts_v0.6')));
% addpath(genpath(fullfile(fp_scr,'packages','spm12')));
addpath(genpath(fullfile(fp_scr,'packages','spm12', 'toolbox', 'suit')));
% 
% cd(fullfile(fp_scr,'packages', 'tapas-6.0.0'))
% tapas_init
% cd(fp_pwd)

path_2_dcm2niix = fullfile(fp_scr,'packages','dcm2niix','build','bin');
myPath = getenv('PATH');
if not(contains(myPath,path_2_dcm2niix))
  setenv('PATH',[ path_2_dcm2niix pathsep myPath]);
end
clearvars path_2_dcm2niix myPath i dl2chk
