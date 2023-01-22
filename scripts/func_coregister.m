function [fp_func_norm] = func_coregister_normalize(fp_moving, fp_reference, fp_other, useMod)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%   Remaining issues:
%       1. display resampled for CAT12 tool
%       2. fp_other should also work for QSM (fp_func_mean)

if ~exist('useMod', 'var')
    useMod = true;
end

if useMod
    fp_reference = func_changeFileName_PrePostFileext(fp_reference, 'm');
end
fp_func_norm = func_changeFileName_PrePostFileext(fp_other, 'w');

fl_reference = func_fp2fl_Texpand(fp_reference);
fl_other = func_fp2fl_Texpand(fp_other);
fl_func_mean = func_fp2fl_Texpand(fp_moving);

%% Coregister
matlabbatch{1}.spm.spatial.coreg.estimate.ref = fl_reference;
matlabbatch{1}.spm.spatial.coreg.estimate.source = fl_func_mean;
if exist('fl_other', 'var')
    matlabbatch{1}.spm.spatial.coreg.estimate.other = fl_other;
end
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

%% Combine matlabbatch's and save
[fp_anat, ~, ~] = fileparts(fp_reference);

fp_pl = fullfile(fp_anat, '..', 'processing_logs');
mkdir(fp_pl)
cd(fp_pl)
save(fullfile(fp_pl, ...
    ['matlabbatch_', ...
    'coNorm_', ...
    't-', datestr(now, 'yymmddTHHMMSS'), ...
    '.mat']), 'matlabbatch')

%% Run
spm_jobman('run', matlabbatch)

end

