function fl_suit = func_SUIT(fl_anat)
% FUNC_SUIT: Does SUIT processing for input anatomical files
%   Each anatomical file should have a separate output SUIT directory
%   created for them (including a copy of the original file). This SUIT
%   directory should then be filled with the segmentation files, norm

%% Initialization
fp_packages = func_getPackagesFolder();
addpath(genpath(fullfile(fp_packages,'spm12', 'toolbox', 'suit')));
startspm 12
if ~iscell(fl_anat)
    fl_anat = {fl_anat};
end
fl_anat_parsed = essbids_parseLabel(fl_anat);
for iANAT = 1:length(fl_anat)
    if ~isfield(fl_anat_parsed{iANAT}, 'acq')
        fl_anat_parsed{iANAT}.acq = 'base';
    end
end

%% SUIT segmentation
% prepare SUIT processing (copyfile, reorient, gunzip)
nSUIT = length(fl_anat);
fl_suit = cell(nSUIT, 1);
fl_suit_folder = cell(nSUIT, 1);
fl_suit_run = false(nSUIT, 1);
job_normalize_cell = cell(nSUIT, 1);
job_reslice_cell = cell(nSUIT, 1);
job_resliceWM_cell = cell(nSUIT, 1);
job_resliceInv_cell = cell(nSUIT, 1);
for iSUIT = 1:nSUIT
    fp_anat = fl_anat{iSUIT};
    
    % define folder
    fl_suit_folder{iSUIT} = fullfile(fl_anat_parsed{iSUIT}.fpath, ...
        ['SUIT_ses-', fl_anat_parsed{iSUIT}.ses, '_acq-', fl_anat_parsed{iSUIT}.acq]);
    if ~exist(fl_suit_folder{iSUIT}, 'dir') 
        fl_suit_run(iSUIT) = true;
        [~, ~] = mkdir(fl_suit_folder{iSUIT});
    end
    
    % copy original anat file to SUIT processing folder
    fp_suit = fullfile(fl_suit_folder{iSUIT}, [fl_anat_parsed{iSUIT}.fname, fl_anat_parsed{iSUIT}.extension]);
    copyfile(fp_anat, fp_suit)
    
    % % fslreorient needed?
    % system('fslreorient2std, fl_suit1');
    % cellfun(@(x, y) system(['fslreorient2std ', x, ' ', y]), fl_suit1, fl_suitReo, 'uni', 0);
    % fl_suitReo = func_changeFileName_PrePostFileext(fl_suit1, '', '_reoriented', '.nii');
    
    % gunzip if necessary
    if strcmp(fl_anat_parsed{iSUIT}.extension, '.nii.gz')
        fl_suit(iSUIT) = gunzip(fp_suit);
    else
        fl_suit{iSUIT} = fp_suit;
    end
end

% Prepare SUIT jobs
% Output of isolation used by normalization
fl_suit_gm = func_changeFileName_PrePostFileext(fl_suit, '', '_seg1', '.nii');
fl_suit_wm = func_changeFileName_PrePostFileext(fl_suit, '', '_seg2', '.nii');
fl_suit_pc = func_changeFileName_PrePostFileext(fl_suit, 'c_', '_pcereb', '.nii');
% fl_suitReo_cr = func_changeFileName_PrePostFileext(fl_suit1, 'c_', '', '.nii');

% Output of normalization used by reslicing
fl_suit_af = func_changeFileName_PrePostFileext(fl_suit, 'Affine_', '_seg1', '.mat');
fl_suit_ff = func_changeFileName_PrePostFileext(fl_suit, 'u_a_', '_seg1', '.nii');

% fl_suitReo = gunzip(fl_suitReo)';
% fl_suit1 = cellfun(@(x, y) [x ',1'], fl_suit1, 'uni', 0);
% fl_suitReo = cellfun(@(x, y) [x ',1'], fl_suitReo, 'uni', 0);
% matlabbatch{1}.spm.tools.suit.isolate_seg.source = cellfun(@(x) {{x}}, fl_suit1');
% fl_suit1 = fl_suitReo(end);

for iSUIT = 1:nSUIT
    
    job_normalize_cell{iSUIT}.subjND(1).gray = fl_suit_gm(iSUIT); % Gray matter image.
    job_normalize_cell{iSUIT}.subjND(1).white = fl_suit_wm(iSUIT); % White matter image.
    job_normalize_cell{iSUIT}.subjND(1).isolation = fl_suit_pc(iSUIT); % Mask
    
    job_reslice_cell{iSUIT}.subj.affineTr=fl_suit_af(iSUIT); % Affine transformation
    job_reslice_cell{iSUIT}.subj.flowfield=fl_suit_ff(iSUIT); % This is the flowfield
    job_reslice_cell{iSUIT}.subj.resample{1}=fl_suit_gm{iSUIT}; % Gray matter image.
    job_reslice_cell{iSUIT}.subj(1).mask=fl_suit_pc(iSUIT); % Mask
    job_reslice_cell{iSUIT}.jactransf=1; % Important for VBM
    
    job_resliceWM_cell{iSUIT}.subj.affineTr=fl_suit_af(iSUIT); % Affine transformation
    job_resliceWM_cell{iSUIT}.subj.flowfield=fl_suit_ff(iSUIT); % This is the flowfield
    job_resliceWM_cell{iSUIT}.subj.resample{1}=fl_suit_wm{iSUIT}; % Gray matter image.
    job_resliceWM_cell{iSUIT}.subj(1).mask=fl_suit_pc(iSUIT); % Mask
    job_resliceWM_cell{iSUIT}.jactransf=1; % Important for VBM
    
    job_resliceInv_cell{iSUIT}.Affine = fl_suit_af(iSUIT); % Affine transformation
    job_resliceInv_cell{iSUIT}.flowfield = fl_suit_ff(iSUIT); % This is the flowfield
    job_resliceInv_cell{iSUIT}.resample = {fullfile(fp_packages, 'spm12', 'toolbox', 'suit', 'atlasesSUIT', 'Lobules-SUIT.nii')}; % To resample to native space
    job_resliceInv_cell{iSUIT}.ref = fl_suit(iSUIT); % Important for VBM
end

% Run SUIT
runFail = false(nSUIT, 1);
exception_cell = cell(nSUIT, 1);
for iSUIT = 1:nSUIT
%     disp(iSUIT)
    if fl_suit_run(iSUIT)
        try
            % isolate cerebellar GM and WM
            fprintf('suit_isolate_seg \n')
            suit_isolate_seg(fl_suit(iSUIT)) 

            % normalize to SUIT space
            fprintf('suit_normalize_dartel \n')
            suit_normalize_dartel(job_normalize_cell{iSUIT}) 

            % use transformation to SUIT space for GM VBM
            fprintf('suit_reslice_dartel GM VBM \n')
            suit_reslice_dartel(job_reslice_cell{iSUIT}); 
    %         flatmap_GM = suit_map2surf(); % SUIT flatmap

            % use transformation to SUIT space for WM VBM
            fprintf('suit_reslice_dartel WM VBM \n')
            suit_reslice_dartel(job_resliceWM_cell{iSUIT});
    %         flatmap_WM = suit_map2surf(); % SUIT flatmap 

            % SUIT annotation to native space (inverse SUIT transform)
            fprintf('suit_reslice_dartel_inv \n')
            suit_reslice_dartel_inv(job_resliceInv_cell{iSUIT})

        catch exception
            exception_cell{iSUIT} = exception;
            runFail(iSUIT) = true;
        end
    end
end

%% Smoothing
fprintf('Smoothing \n')
clear matlabbatch

% Settings
smoothing_arr = [2,4,6,8];

matlabbatch_run_bool = false(2*length(smoothing_arr), 1);

% gray matter probability smoothing
fl_wdm_seg1 = func_changeFileName_PrePostFileext(fl_anat, 'wd', '_seg1', '.nii,1');
for i = 1:length(smoothing_arr)
    smoothing_prefix = ['s', num2str(smoothing_arr(i))];
    
    fl_wdm_seg1_smoothed = func_changeFileName_PrePostFileext(fl_anat, ...
        [smoothing_prefix, 'wd'], '_seg1', '.nii,1');
    fl_wdm_seg1_smoothed_run = cellfun(@(x) ~exist(x, 'file'), ...
        fl_wdm_seg1_smoothed, 'uni', 0);
    fl_wdm_seg1_smoothed_run = [fl_wdm_seg1_smoothed_run{:}]';
    
    matlabbatch{i}.spm.spatial.smooth.data = fl_wdm_seg1(fl_wdm_seg1_smoothed_run);
    matlabbatch{i}.spm.spatial.smooth.fwhm = [smoothing_arr(i) smoothing_arr(i) smoothing_arr(i)];
    matlabbatch{i}.spm.spatial.smooth.dtype = 0;
    matlabbatch{i}.spm.spatial.smooth.im = 0;
    matlabbatch{i}.spm.spatial.smooth.prefix = smoothing_prefix;
    
    matlabbatch_run_bool(i) = ~isempty(matlabbatch{i}.spm.spatial.smooth.data);
end

% white matter probability smoothing
fl_wdm_seg2 = func_changeFileName_PrePostFileext(fl_anat, 'wd', '_seg2', '.nii,1');
j = 0;
for i = length(smoothing_arr)+1:2*length(smoothing_arr)
    j = j + 1;
    smoothing_prefix = ['s', num2str(smoothing_arr(j))];
    
    fl_wdm_seg2_smoothed = func_changeFileName_PrePostFileext(fl_anat, ...
        [smoothing_prefix, 'wd'], '_seg2', '.nii,1');
    fl_wdm_seg2_smoothed_run = cellfun(@(x) ~exist(x, 'file'), ...
        fl_wdm_seg2_smoothed, 'uni', 0);
    fl_wdm_seg2_smoothed_run = [fl_wdm_seg2_smoothed_run{:}]';

    matlabbatch{i}.spm.spatial.smooth.data = fl_wdm_seg2(fl_wdm_seg2_smoothed_run);
    matlabbatch{i}.spm.spatial.smooth.fwhm = [smoothing_arr(j) smoothing_arr(j) smoothing_arr(j)];
    matlabbatch{i}.spm.spatial.smooth.dtype = 0;
    matlabbatch{i}.spm.spatial.smooth.im = 0;
    matlabbatch{i}.spm.spatial.smooth.prefix = smoothing_prefix;
    
    matlabbatch_run_bool(i) = ~isempty(matlabbatch{i}.spm.spatial.smooth.data);
end

matlabbatch = matlabbatch(matlabbatch_run_bool);

%% Combine matlabbatch's and save
for iSUIT = 1:length(fl_anat)
    fp_anat = fl_suit{iSUIT};
    fp_anat_folder = strsplit(fp_anat, filesep);
    fp_anat_folder = strjoin(fp_anat_folder(1:find(strcmp(fp_anat_folder, 'anat'))), filesep);

    fp_pl = fullfile(fp_anat_folder, '..', 'processing_logs');
    if ~exist(fp_pl, 'dir')
        [~, ~] = mkdir(fp_pl);
    end
    save(fullfile(fp_pl, ...
        ['matlabbatch_', ...
        'SUITsmooth_', ...
        't-', datestr(now, 'yymmddTHHMMSS'), ...
        '.mat']), 'matlabbatch')
end

end
