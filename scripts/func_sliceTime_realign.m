function  [fp_func, fp_func_mean, fp_reg] = func_sliceTime_realign(fp_func)
%FUNC_SLICETIME_REALIGN Do slicetiming and realignment for bold measurement
%   Process bold measurement with optional slice timing correction and
%   non-optional motion correction in the form of a realignment step. The
%   realignment is to the first measurement of the bold measurement, a mean
%   EPI image is also created.

% split based on '.', separate fp from fileExtensions
fp_split = strsplit(fp_func, '.');

% read in json
fp_json = [fp_split{1}, '.json'];
jinfo = jsondecode(fileread(fp_json));

% final output has r or ra prefix depending on whether slicetiming
% correction is done, if this file exists you can skip the rest
if isfield(jinfo, 'SliceTiming')
    fp_func_output = func_changeFileName_PrePostFileext(fp_func, 'ra', '', 'nii');
else
    fp_func_output = func_changeFileName_PrePostFileext(fp_func, 'r', '', 'nii');
end

if ~exist(fp_func_output, 'file')

if strcmp(fp_split{end}, 'gz')
    fp_func = gunzip(fp_func);
end
fp_func = fp_func{1};

fl_func = func_fp2fl_Texpand(fp_func);

%% sliceTime
if isfield(jinfo, 'SliceTiming')
    matlabbatch_sliceTime{1}.spm.temporal.st.scans = fl_func;
    %%
    matlabbatch_sliceTime{1}.spm.temporal.st.nslices = length(jinfo.SliceTiming);
    matlabbatch_sliceTime{1}.spm.temporal.st.tr = jinfo.RepetitionTime;
    matlabbatch_sliceTime{1}.spm.temporal.st.ta = 0;
    matlabbatch_sliceTime{1}.spm.temporal.st.so = jinfo.SliceTiming;
    matlabbatch_sliceTime{1}.spm.temporal.st.refslice = 0;
    matlabbatch_sliceTime{1}.spm.temporal.st.prefix = 'a';
    fp_func = func_changeFileName_PrePostFileext(fp_func, 'a');
    fl_func = func_changeFileName_PrePostFileext(fl_func, 'a', '', 'nii');
else
    matlabbatch_sliceTime = [];
end

%% realignment
matlabbatch{1}.spm.spatial.realign.estwrite.data = {fl_func};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
fp_func_mean = func_changeFileName_PrePostFileext(fp_func, 'mean');
fp_reg = func_changeFileName_PrePostFileext(fp_func, 'rp_', '', 'txt');
fp_func = func_changeFileName_PrePostFileext(fp_func, 'r'); % fp_func = fp_func_output

%% Combine matlabbatch's and save
matlabbatch = [matlabbatch_sliceTime, matlabbatch];
% save(fullfile(fp_func_folder, ...
%     ['matlabbatch_', ...
%     'moCorr_', ...
%     datestr(now, 'yymmddTHHMMSS'), ...
%     '.mat']), 'matlabbatch')
%% Combine matlabbatch's and save
fp_func_folder = strsplit(fp_func);
fp_func_folder = strjoin(fp_func_folder(1:find(strcmp(fp_func_folder, 'func'))), filesep);

fp_pl = fullfile(fp_func_folder, '..', 'processing_logs');
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

%% Warning in case fp_func_output already exists
else
    if isfield(jinfo, 'SliceTiming')
        fp_func = func_changeFileName_PrePostFileext(fp_func, 'a');
    end
    fp_func_mean = func_changeFileName_PrePostFileext(fp_func, 'mean', '', 'nii');
    fp_reg = func_changeFileName_PrePostFileext(fp_func, 'rp_', '', 'txt');
    fp_func = fp_func_output;
    warning([fp_func_output, ' already exists, skipping slicetime_realign'])
end

end

