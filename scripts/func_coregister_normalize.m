function [fp_func_norm, fp_other] = func_coregister_normalize(fp_moving, fp_reference, fp_other, useMod)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%   Remaining issues:
%       1. display resampled for CAT12 tool
%       2. fp_other should also work for QSM (fp_func_mean)
%       3. Not verified to work without fp_other
%

if ~exist('useMod', 'var')
    useMod = true;
end

fp_transform = func_changeFileName_PrePostFileext(fp_reference, 'y_');
if useMod
    fp_reference = func_changeFileName_PrePostFileext(fp_reference, 'm');
end

fp_moving_old = fp_moving;
fp_moving = func_changeFileName_PrePostFileext(fp_moving, 'c');
fp_other_old = fp_other;
fp_other = func_changeFileName_PrePostFileext(fp_other, 'c');
fp_func_norm = func_changeFileName_PrePostFileext(fp_other, 'w');

if ~exist(fp_moving, 'file')
    copyfile(fp_moving_old, fp_moving)
    if exist('fp_other', 'var')
        copyfile(fp_other_old, fp_other)
    end

    fl_other = func_fp2fl_Texpand(fp_other);
    fl_reference = func_fp2fl_Texpand(fp_reference);
    fl_func_mean = func_fp2fl_Texpand(fp_moving);
    fl_norm = {fp_transform};

    fp_func_mean_info = niftiinfo(fp_moving);

    %% Coregister moving to reference, other is coregistered in the same way
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = fl_reference;
    matlabbatch{1}.spm.spatial.coreg.estimate.source = fl_func_mean;
    if exist('fl_other', 'var')
        matlabbatch{1}.spm.spatial.coreg.estimate.other = fl_other;
    end
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    %% Normalize to MNI space with CAT12 transform
    if exist('fl_other', 'var') % should also be able to work without fl_other
        matlabbatch{2}.spm.spatial.normalise.write.subj.def = fl_norm;
        matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cfg_dep('Coregister: Estimate: Coregistered Images', ...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('.','cfiles'));
        matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
                                                                  NaN NaN NaN];
    %     matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    %                                                               78 76 85];
        matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = fp_func_mean_info.PixelDimensions(1:3);
        matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';
    else % there should be a without fl_other option here!
    end

    %% Combine matlabbatch's and save
    fp_anat = strsplit(fp_reference, filesep);
    fp_anat = strjoin(fp_anat(1:find(strcmp(fp_anat, 'anat'))), filesep);
%     [fp_anat, ~, ~] = fileparts(fp_reference);
    
    fp_pl = fullfile(fp_anat, '..', 'processing_logs');
    if ~exist(fp_pl, 'dir')
        mkdir(fp_pl)
    end
    save(fullfile(fp_pl, ...
        ['matlabbatch_', ...
        'sliceTimeRealign_', ...
        't-', datestr(now, 'yymmddTHHMMSS'), ...
        '.mat']), 'matlabbatch')

    %% Run
    spm_jobman('run', matlabbatch)
else
    warning([fp_moving, ' already exists, skipping coregistration and normalization'])
end

end

