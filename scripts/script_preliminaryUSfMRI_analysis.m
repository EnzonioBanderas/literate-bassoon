% This script contains a complete pipeline for 1stlevel fMRI analysis over
% given subjects and sessions
%
% TD:
%   2. Brain mask (only for later visualization)
%   3. 2nd level analysis expansion, CS+>CS- contrast
%   5. ANTs normalization, replace topup
%   6. Some figures and images are generated in what I assume is the
%      working directory. How can this be prevented?
%           1) Change working directory during each function wrapper
%           2) Directly cancel figure generation
% 
% SUIT processing is not really necessary, CAT12 normalization is
% sufficient, SUIT processing likely requires manually drawn cerebellar
% masks in most cases.
%
%% Initialization
script_init_study7
fp_packages = func_getPackagesFolder();

% settings
doStats = true;  
reanalysis = false;
uncertaintyUS_time = 0.5; % time uncertainty around US in s (used for searching for US around a time point)
% SUBoI = 'sub-Z7T7185';
% SESOI = 'ses-1';
% fp_s = fullfile(fp0, 'dumpHereForSorting', 'sourcedata2');
% fp_d = fullfile(fp0, 'dumpHereForSorting', 'rawdata');
% fp_de = fullfile(fp0, 'derivatives', 'prelimUSfMRI_sub-Z7T7285');
fp_deri = fullfile(fp_de, ['prelimUSfMRI_sub-Z7T7268_datetime-', datestr(now, 'yyyymmddTHHMMSS')]);
fp_deri = ['/media/diskEvaluation/Evaluation/sfb1280a05study7/derivatives/', ...
    'prelimRerunTopupfMRI'];
fp_2ndLevel = fullfile(fp_de, 'prelimRerunTopupfMRI_2ndLevel');
[~, ~] = mkdir(fp_2ndLevel);
firstLevelName = '1stLevel_4p5_addConDebug2';
smoothing_FWHM = 4.5;

% % packages
% addpath(genpath(fullfile(fp_scr,'packages','spm12')));

% Clear fp_deri if it exists, go to fp_deri 
% if exist(fp_de, 'dir')
%     rmdir(fp_de, 's')
% end
[~, ~] = mkdir(fp_deri);
cd(fp_deri)

% identify fear data fp_rawd and copy whole sessions to fp_deri
fl_fear = dir(fullfile(fp_d, '*', '*', '*', '*task-fear*bold*'));
fl_fear_subses = cell(length(fl_fear), 2);
for iFL=1:length(fl_fear)
    fl_fear_folderSplit = strsplit(fl_fear(iFL).folder, filesep);
    fl_fear_subses(iFL, :) = fl_fear_folderSplit([end-2, end-1]);
end

% optional filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub_IDs = {'sub-Z7T7333', 'sub-Z7T7337', 'sub-Z7T7338'}; % debug2 % MAYBE DONE
% sub_IDs = {'sub-Z7T7330', 'sub-Z7T7331', 'sub-Z7T7332'}; % debug2 % MAYBE DONE
% sub_IDs = {'sub-Z7T7316', 'sub-Z7T7317', 'sub-Z7T7318'}; % debug2 % MAYBE DONE
% sub_IDs = {'sub-Z7T7311', 'sub-Z7T7312', 'sub-Z7T7313'}; % debug2 % MAYBE DONE
% sub_IDs = {'sub-Z7T7300', 'sub-Z7T7301', 'sub-Z7T7302'}; % debug2 % MAYBE DONE
sub_IDs = {'sub-Z7T7311', 'sub-Z7T7312', 'sub-Z7T7313', ...
    'sub-Z7T7316', 'sub-Z7T7317', 'sub-Z7T7318', ...
    'sub-Z7T7330', 'sub-Z7T7331', 'sub-Z7T7332', ...
    'sub-Z7T7333', 'sub-Z7T7337', 'sub-Z7T7338'}; % debug2 % MAYBE DONE
% sub_IDs = {'sub-Z7T7300', 'sub-Z7T7301', 'sub-Z7T7302'}; % debug2 % MAYBE DONE
% sub_IDs = {'sub-Z7T7295', 'sub-Z7T7296'}; % debug2 % DONE 297 error volatile missing images {'sub-Z7T7295', 'sub-Z7T7296', 'sub-Z7T7297'};
% sub_IDs = {'sub-Z7T7268', 'sub-Z7T7269', 'sub-Z7T7284', 'sub-Z7T7285'}; % DONE
% sub_IDs = {'sub-Z7T7267', 'sub-Z7T7262', 'sub-Z7T7253', 'sub-Z7T7252'}; % DONE
% sub_IDs = {'sub-Z7T7233', 'sub-Z7T7232', 'sub-Z7T7229', 'sub-Z7T7228'}; % DONE
% sub_IDs = {'sub-Z7T7297'}; % DONE likely
% sub_IDs = {'sub-Z7T7216', 'sub-Z7T7217'}; % DONE likely
% sub_IDs = {'sub-Z7T7219', 'sub-Z7T7220'}; % DONE likely
% sub_IDs = {'sub-Z7T7227', 'sub-Z7T7221'}; % DONE likely
% 215 214 213 205 201 186 185 184 
% sub_IDs = {'sub-Z7T7215', 'sub-Z7T7214'}; % DONE, although not fully
% sub_IDs = {'sub-Z7T7213', 'sub-Z7T7205'}; % DONE, although not fully (213 lacks ses-3)
% sub_IDs = {'sub-Z7T7186'}; % 'sub-Z7T7201' has a unique CAT12 error, reprocess MP2RAGE with FATNAV corrections instead
% sub_IDs = {'sub-Z7T7201'};
% sub_IDs = {'sub-Z7T7185', 'sub-Z7T7184'};
% sub_IDs = {'sub-Z7T7180', 'sub-Z7T7179'};
% sub_IDs = {'sub-Z7T7175'};
% 180 179 175 
% 173 172 162 161 % These have a lot of problems!
% sub_IDs = {'sub-Z7T7180'};
%
% 
% sub_IDs = {'sub-Z7T7185', 'sub-Z7T7184', 'sub-Z7T7180', 'sub-Z7T7179', 'sub-Z7T7175'};
fl_fear_subses = fl_fear_subses(ismember(fl_fear_subses(:,1), sub_IDs), :);
% fl_fear_subses = fl_fear_subses(ismember(fl_fear_subses(:,2), 'ses-1'), :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% copy files and only use copied rawdata
[~, fl_fear_subses_ind] = unique(cellfun(@(x, y) [x, y], fl_fear_subses(:, 1), fl_fear_subses(:, 2), 'uni', 0));
fl_fear_subses = fl_fear_subses(fl_fear_subses_ind, :);
dataAlreadyExists = false(size(fl_fear_subses, 1), 1);
fp_rawd_subses = cell(size(fl_fear_subses, 1), 1);
fp_deri_subses = cell(size(fl_fear_subses, 1), 1);
for iSS = 1:size(fl_fear_subses, 1)
    fp_rawd_subses{iSS} = fullfile(fp_d, fl_fear_subses{iSS, 1}, fl_fear_subses{iSS, 2});
    fp_deri_subses{iSS} = fullfile(fp_deri, fl_fear_subses{iSS, 1});
    dataAlreadyExists(iSS) = exist(fp_deri_subses{iSS}, 'dir');
end
if all(dataAlreadyExists) % no data already exists
    warning('skipping copyfile, data of at least some of all folders already exists')
else
    for iSS = 1:size(fl_fear_subses, 1)
%         if ~exist(fp_deri_subses{iSS}, 'dir')
%             [~, ~] = mkdir(fp_deri_subses{iSS});
%             copyfile(fp_rawd_subses{iSS}, fp_deri_subses{iSS});
%         end
        % come up with algorithm that skips already existing files
        % for each of the subfolders in source, create a subfolder in
        % destination if it does not exist already and after looping 
        % through subfolders, copy file from source to destination if the 
        % file does not exists in there already
        copyfile(fp_rawd_subses{iSS}, fp_deri_subses{iSS});
    end
end

% for each subject-session which contains MRI func *task-fear* data and events table
fp_deri_subses = unique(fp_deri_subses);
table_subses = cell(length(fp_deri_subses), 1);
%% Loop over subjects (which now contain all ses data in one folder)
for iSS = 1:length(fp_deri_subses)
    
    %% Define
    fp_anat = fullfile(fp_deri_subses{iSS}, 'anat');
    fp_1stLevel = fullfile(fp_deri_subses{iSS}, firstLevelName);
    [~, ~] = mkdir(fp_1stLevel);
    fp_pl = fullfile(fp_deri_subses{iSS}, 'processing_logs');
    [~, ~] = mkdir(fp_pl);
    fp_fig = fullfile(fp_deri_subses{iSS}, 'figures');
    [~, ~] = mkdir(fp_fig);

    %% CAT12 anat
    % clean MP2RAGE
    fp_anat_clean = func_MP2RAGE_clean(fp_anat);

    % run CAT12
%     cat12('expert')
    fp_anat_CAT12 = cell(length(fp_anat_clean), 1);
    fp_anat_CAT12_acq = cell(length(fp_anat_clean), 1);
    catlog_IQR = zeros(length(fp_anat_clean), 1);
    for iANAT = 1:length(fp_anat_clean)
        % parse cleaned MP2RAGE and create its CAT12 output folder
        fp_anat_clean_parsed = essbids_parseLabel(fp_anat_clean{iANAT}, 'silent', true);
        if ~isfield(fp_anat_clean_parsed, 'acq')
            fp_anat_clean_parsed.acq = 'base';
            fp_anat_CAT12_acq{iANAT} = 'base';
        else
            fp_anat_CAT12_acq{iANAT} = fp_anat_clean_parsed.acq;
        end
        fp_anat_CAT12{iANAT} = fullfile(fp_anat_clean_parsed.fpath, ['CAT12_', 'ses-', fp_anat_clean_parsed.ses, '_acq-' fp_anat_clean_parsed.acq]);

        % updated cleaned MP2RAGE location within CAT12 output folder
        fp_anat_clean{iANAT} = fullfile(fp_anat_CAT12{iANAT}, [fp_anat_clean_parsed.fname, '.nii']);

%         % if the directory already exists, skip CAT12 processing
%         if ~exist(fp_anat_CAT12{iANAT}, 'dir')
        % if the modulated anat file already exists within the CAT12
        % directory, skip CAT12 processing
        if ~exist(func_changeFileName_PrePostFileext(fp_anat_clean{iANAT}, 'm'), 'file')
            [~, ~] = mkdir(fp_anat_CAT12{iANAT});

            % copy cleaned MP2RAGE nii and json to CAT12, update fp_anat_clean variable
            copyfile(fullfile(fp_anat_clean_parsed.fpath, [fp_anat_clean_parsed.fname, '.nii']), ...
                     fullfile(fp_anat_CAT12{iANAT}))
            copyfile(fullfile(fp_anat_clean_parsed.fpath, [fp_anat_clean_parsed.fname, '.json']), ...
                     fullfile(fp_anat_CAT12{iANAT}))

            % CAT12 processing 
            func_CAT12(fp_anat_clean(iANAT))
        end

        % Get IQR
        fp_anat_catlog = func_changeFileName_PrePostFileext(fp_anat_clean{iANAT}, 'catlog_', '', 'txt');
        catlog = fileread(fp_anat_catlog);
        catlog = regexp(catlog, 'IQR(.*?)%', 'match');
        catlog_IQR(iANAT) = str2double(catlog{1}(end-5:end-1));
    end

    % choose CAT12 which is not base, if present, else choose the one
    % which has highest IQR
    nonBase_ind = find(~ismember(fp_anat_CAT12_acq, 'base'));
    if isempty(nonBase_ind) % if there are only base MP2RAGE reconstructions
        [~, anatInd] = max(catlog_IQR); % use the one with the highest IQR
        fp_anat_clean_chosen = fp_anat_clean{anatInd};
        fp_anat_mod = func_changeFileName_PrePostFileext(fp_anat_clean_chosen, ...
            'm', '', 'nii');
    else % if there are non-base MP2RAGE reconstructions
        [~, catlog_IQR_maxInd] = max(catlog_IQR(nonBase_ind)); % use the one with the highest IQR
        anatInd = nonBase_ind(catlog_IQR_maxInd);
        fp_anat_clean_chosen = fp_anat_clean{anatInd};
        fp_anat_mod = func_changeFileName_PrePostFileext(fp_anat_clean_chosen, ...
            'm', '', 'nii');
    end
    
    %% Do SUIT processing
    % use fp_anat_mod, if something goes wrong with this might want to
    % switch to cleaned (or even before cleaning)
    fl_anat_SUIT = func_SUIT({fp_anat_clean_chosen});
    fp_anat_SUIT = fl_anat_SUIT{1};
    
    %% func loop
    % sliceTiming_realign matlabbatch
    fl_func = dir(fullfile(fp_deri_subses{iSS}, 'func', 'sub*bold.nii.gz'));
    fl_ev = func_dirl2fl(dir(fullfile(fp_deri_subses{iSS}, 'func', '*_events.tsv')));
    
    fl_func_preprocessed = cell(length(fl_func), 1);
    fl_func_mean = cell(length(fl_func), 1);
    fp_reg = cell(length(fl_func), 1);
    fl_func_dc = cell(length(fl_func_preprocessed), 1);
    fl_func_mean_dc = cell(length(fl_func_preprocessed), 1);
    fl_func_coreg = cell(length(fl_func_preprocessed), 1);
    fl_func_norm = cell(length(fl_func_preprocessed), 1);
    fp_func = cell(length(fl_func_preprocessed), 1);
    fp_func_parsed = cell(length(fl_func_preprocessed), 1);
    fp_func_phase = cell(length(fl_func_preprocessed), 1);
    fn_func = cell(length(fl_func_preprocessed), 1);
    LCEg_spmT_mean = zeros(length(fl_func_preprocessed), 1);
    RCEg_spmT_mean = zeros(length(fl_func_preprocessed), 1);
    LCEg_tsnr_mean = zeros(length(fl_func_preprocessed), 1);
    RCEg_tsnr_mean = zeros(length(fl_func_preprocessed), 1);
%     for iFUNC = 3
    for iFUNC = 1:length(fl_func)
        %% fl_func parse
        fp_func{iFUNC} = fullfile(fl_func(iFUNC).folder, fl_func(iFUNC).name);
        fp_func_parsed{iFUNC} = essbids_parseLabel(fp_func{iFUNC});
        [~, fn_func{iFUNC}] = fileparts(fp_func{iFUNC});
%         fp_func_niftiinfo = func_niftiinfo(fp_func{iFUNC});
        fp_func_json = func_changeFileName_PrePostFileext(fp_func{iFUNC}, '', '', 'json');
        fp_func_jinfo = jsondecode(fileread(fp_func_json));
        switch [fp_func_parsed{iFUNC}.task, fp_func_parsed{iFUNC}.ses, fp_func_parsed{iFUNC}.run]
            case 'fear11'
                fp_func_phase{iFUNC} = 'Habituation'; 
            case 'fear12'
                fp_func_phase{iFUNC} = 'Acquisition';
            case 'fear21'
                fp_func_phase{iFUNC} = 'Extinction';
            case 'fear31'
                fp_func_phase{iFUNC} = 'Recall';
            case 'fear32'
                fp_func_phase{iFUNC} = 'Volatile';
            case 'rest11'
                fp_func_phase{iFUNC} = 'preAcquisition'; 
            case 'rest12'
                fp_func_phase{iFUNC} = 'postAcquisition';
            case 'rest21'
                fp_func_phase{iFUNC} = 'preExtinction';
            case 'rest22'
                fp_func_phase{iFUNC} = 'postExtinction';
            case 'rest31'
                fp_func_phase{iFUNC} = 'postVolatile';
            otherwise
                fp_func_phase{iFUNC} = 'Unknown';
        end

        %% slice-time correction and movement: preprocess func
        [fl_func_preprocessed{iFUNC}, fl_func_mean{iFUNC}, fp_reg{iFUNC}] = ...
            func_sliceTime_realign(fp_func{iFUNC});

        %% distortion correction: topup func
        % get preprocessed func outputs
        fp_func_preprocessed = fl_func_preprocessed{iFUNC};
        fp_func_mean = fl_func_mean{iFUNC};
        fp_func_preprocessed_parsed = essbids_parseLabel(fp_func_preprocessed, 'silent', true);

        % get topup output files
        fp_func_dc = func_changeFileName_PrePostFileext(fp_func_preprocessed, '', '_dc');
        fp_func_mean_dc = func_changeFileName_PrePostFileext(fp_func_mean, '', '_dc');

        if ~exist(fp_func_dc, 'file')

            % create topup folder
            fp_topup = fullfile(fp_deri_subses{iSS}, 'func', 'topup');
            if isfield(fp_func_preprocessed_parsed, 'acq')
                fp_topup = [fp_topup, '_acq-', fp_func_preprocessed_parsed.acq];
            end
            if isfield(fp_func_preprocessed_parsed, 'ses')
                fp_topup = [fp_topup, '_ses-', fp_func_preprocessed_parsed.ses];
            end
            if isfield(fp_func_preprocessed_parsed, 'run')
                fp_topup = [fp_topup, '_run-', fp_func_preprocessed_parsed.run];
            end
            mkdir(fp_topup)

            % match func to fmap (update func_matchFiles with fmap option?)
            fl_fmap = dir(fullfile(fp_deri_subses{iSS}, 'fmap', 'sub*epi.nii.gz'));
            for iFMAP = 1:length(fl_fmap)
                fp_fmap = fullfile(fl_fmap(iFMAP).folder, fl_fmap(iFMAP).name);
                fp_fmap_parsed = essbids_parseLabel(fp_fmap);
                try
                    if strcmp(fp_fmap_parsed.run, fp_func_preprocessed_parsed.run) && ...
                       strcmp(fp_fmap_parsed.acq, [fp_func_preprocessed_parsed.acq, ...
                                                   fp_func_preprocessed_parsed.task])
                        disp('break');
                        break
                    end
                catch
                end
            end

            % copyfiles to topup folder
            copyfile(fullfile(fp_packages, 'topup', 'acq_param.txt'), fp_topup)
            copyfile(fp_func_preprocessed, fp_topup)
            copyfile(fp_func_mean, fp_topup)
            copyfile(fp_fmap, fp_topup)

            %
            [~, fp_name, fp_ext] = fileparts(fp_func_preprocessed);
            fp_func_preprocessed_topup = fullfile(fp_topup, [fp_name, fp_ext]);
            [~, fp_name, fp_ext] = fileparts(fp_func_mean);
            fp_func_mean_topup = fullfile(fp_topup, [fp_name, fp_ext]);
            [~, fp_name, fp_ext] = fileparts(fp_fmap);
            fp_fmap_topup = fullfile(fp_topup, [fp_name, fp_ext]);

            % make sure there are even number of voxels in all spatial
            % dimensions for both bold, bold mean and epi
            fp_func_preprocessed_info = func_niftiinfo(fp_func_preprocessed_topup);
            fp_func_preprocessed_info_IS = fp_func_preprocessed_info.ImageSize;
            for iDIM=1:3
                % if DIM is odd, subtract 1 and make it even
                if rem(fp_func_preprocessed_info.ImageSize(iDIM), 2) == 1
                    fp_func_preprocessed_info_IS(iDIM) = fp_func_preprocessed_info_IS(iDIM) - 1;
                end
            end
            cmd = ['fslroi ', fp_func_preprocessed_topup, ' ', fp_func_preprocessed_topup, ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(1)), ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(2)), ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(3))];
            system(cmd);
            cmd = ['fslroi ', fp_func_mean_topup, ' ', fp_func_mean_topup, ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(1)), ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(2)), ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(3))];
            system(cmd);
            cmd = ['fslroi ', fp_fmap_topup, ' ', fp_fmap_topup, ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(1)), ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(2)), ...
                ' 0 ', num2str(fp_func_preprocessed_info_IS(3))];
            system(cmd);
            fp_func_preprocessed_topup = gzip(fp_func_preprocessed_topup);
            fp_func_preprocessed_topup = fp_func_preprocessed_topup{1};
            fp_func_mean_topup = gzip(fp_func_mean_topup);
            fp_func_mean_topup = fp_func_mean_topup{1};

            [~, fp_name, fp_ext] = fileparts(fp_func_preprocessed_topup);
            fp_func_preprocessed_topup = [fp_name];
            [~, fp_name, fp_ext] = fileparts(fp_func_mean_topup);
            fp_func_mean_topup = [fp_name];
            [~, fp_name, fp_ext] = fileparts(fp_fmap_topup);
            fp_fmap_topup = [fp_name];
            [~, fp_name, fp_ext] = fileparts(fp_func_preprocessed_topup);
            fp_func_preprocessed_topup = [fp_name];
            [~, fp_name, fp_ext] = fileparts(fp_func_mean_topup);
            fp_func_mean_topup = [fp_name];
            [~, fp_name, fp_ext] = fileparts(fp_fmap_topup);
            fp_fmap_topup = [fp_name];

            % topup script
            topup_cmd = ['bash ', ...
                fullfile(fp_packages, 'topup', 'topup'), ' ', ...
                fp_topup, ' ', ...
                fp_func_preprocessed_topup, ' ', ...
                fp_fmap_topup, ' ', ...
                fp_func_mean_topup];
            system(topup_cmd);

            % move distortion corrected output to upper folder
            fp_func_dc = fullfile(fp_topup, func_changeFileName_PrePostFileext(fp_func_preprocessed_topup, ...
                '', '_dc', 'nii.gz'));
            fl_func_dc(iFUNC) = gunzip(fp_func_dc, fullfile(fp_deri_subses{iSS}, 'func'));
            fp_func_mean_dc = fullfile(fp_topup, func_changeFileName_PrePostFileext(fp_func_mean_topup, ...
                '', '_dc', 'nii.gz'));
            fl_func_mean_dc(iFUNC) = gunzip(fp_func_mean_dc, fullfile(fp_deri_subses{iSS}, 'func'));
        else % get filepaths for next parts, without having to rerun topup!
            fl_func_dc{iFUNC} = fp_func_dc;
            fl_func_mean_dc{iFUNC} = fp_func_mean_dc;
        end

        %% coregister
        %% brain mask
        %% stats
        %% normalize
        % coregister func to anat and normalize
%         coregister_normalize matlabbatch
        [fl_func_norm{iFUNC}, fl_func_coreg{iFUNC}] = func_coregister_normalize(...
            fl_func_mean_dc{iFUNC}, fp_anat_clean{anatInd}, fl_func_dc{iFUNC});
%         
%         % SUIT normalization
%         job_reslice.subj(1).affineTr = ...
%             {func_changeFileName_PrePostFileext(fp_anat_SUIT, ...
%             'Affine_', '_seg1', '.mat')}; % Affine transformation
%         job_reslice.subj(1).flowfield = ... 
%             {func_changeFileName_PrePostFileext(fp_anat_SUIT, ...
%             'u_a_', '_seg1', '.nii')}; % Flow field
%         job_reslice.subj(1).resample = ...
%             fl_func_coreg(iFUNC); % Func image to be resliced
%         job_reslice.subj(1).mask = ...
%             {func_changeFileName_PrePostFileext(fp_anat_SUIT, ...
%             'c_', '_pcereb', '.nii')}; % Mask
%         fprintf('suit_reslice_dartel coregistered to SUIT space\n')
%         suit_reslice_dartel(job_reslice)

        %% func now coregistered and normalized, do brain mask (for visualization)
        
        %% run everything for both SUIT and CAT12 normalizations!
        
        %% 1stLevel - calculate tSNR maps for native, native-distortion-corrected and normalized func images
        if ~exist(func_changeFileName_PrePostFileext(fullfile(fl_func(iFUNC).folder, fl_func(iFUNC).name), 'tsnr_', '', 'nii'), 'file')
            spm_imcalc(func_changeFileName_PrePostFileext(fullfile(fl_func(iFUNC).folder, fl_func(iFUNC).name), '', '', 'nii'), ...
                func_changeFileName_PrePostFileext(fullfile(fl_func(iFUNC).folder, fl_func(iFUNC).name), 'tsnr_', '', 'nii'), ...
                'nanmean(X)./nanstd(X)', struct('dmtx', 1))
        end
        if ~exist(func_changeFileName_PrePostFileext(fl_func_dc{iFUNC}, 'tsnr_', '', 'nii'), 'file')
            spm_imcalc(func_changeFileName_PrePostFileext(fl_func_dc{iFUNC}, '', '', 'nii'), ...
                func_changeFileName_PrePostFileext(fl_func_dc{iFUNC}, 'tsnr_', '', 'nii'), ...
                'nanmean(X)./nanstd(X)', struct('dmtx', 1))
        end
        if ~exist(func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'tsnr_', '', 'nii'), 'file')
        spm_imcalc(func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, '', '', 'nii'), ...
            func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'tsnr_', '', 'nii'), ...
            'nanmean(X)./nanstd(X)', struct('dmtx', 1))
        end
    
        %% 1stLevel - get multcond
        fp_ev = func_matchFiles(fl_ev, fullfile(fl_func(iFUNC).folder, fl_func(iFUNC).name));    
        if ~strcmp(fp_func_parsed{iFUNC}.task, 'rest') && ~isempty(fp_ev)
            [~, fn_multcond] = fileparts(func_changeFileName_PrePostFileext(fp_ev, ...
                '', '_multcond', '.mat'));
            fp_multcond = fullfile(fullfile(fp_1stLevel, [fn_multcond, '.mat']));
            if ~exist(fp_multcond, 'file')
                % Read in event table
                data_ev = essbids_readTsv(fp_ev);
                data_ev = table(data_ev.onset, data_ev.duration, data_ev.trial_type, ...
                    'VariableNames', {'onset', 'duration', 'trial_type'});
                
%                 % divide US: USpostCSplus, USpostCSminus
%                 logical_US = ismember(data_ev.trial_type, {'US'});
%                 ind_US = find(logical_US);
%                 for iUS = ind_US
%                     iCS = iUS - 1;
% %                     while ~ismember(data_ev.trial_type{iCS}, {'CSplus', 'CSminus'})
% %                         iCS = iCS - 1;
% %                     end
%                     data_ev.trial_type{iUS} = ...
%                         [data_ev.trial_type{iUS}, 'post', data_ev.trial_type{iCS}]; 
%                 end
                
                % divide US: USpostCSplus, noUSpostCSplus and noUSpostCSminus
                logical_US = ismember(data_ev.trial_type, {'US'});
                data_ev_US = data_ev(logical_US, :);
                if any(logical_US)
                    data_ev.trial_type(logical_US) = {'USpostCSplus'};
                    data_ev_US_time = data_ev_US.onset;
                else
                    data_ev_US_time = [];
                end
                ind_CS = find(ismember(data_ev.trial_type, {'CSplus', 'CSminus'}));
                for iCS = ind_CS'
                    expectedUS_time = data_ev.onset(iCS) + data_ev.duration(iCS);
                    if isempty( data_ev_US( ... % if no US found after CS ...
                            data_ev_US_time > expectedUS_time - uncertaintyUS_time & ...
                            data_ev_US_time < expectedUS_time + uncertaintyUS_time, :))
                        % ... then add noUSpostCS(+/-)
                        data_ev = [data_ev; cell2table( ...
                            {expectedUS_time, 0, ['noUSpost', data_ev.trial_type{iCS}]}, ...
                            'VariableNames', {'onset', 'duration', 'trial_type'})];
                    end
                end
                % resort data_ev from early to late onsets
                [~, ind_sort] = sort(data_ev.onset);
                data_ev = data_ev(ind_sort, :);
                
                if strcmp(fp_func_phase{iFUNC}, 'Volatile')
                    % first three noUSpostCSplus are unexpected omissions (+two?)
                    ind_unexpNoUS = find(strcmp(data_ev.trial_type, 'noUSpostCSplus'));
                    ind_unexpNoUS = ind_unexpNoUS(1:3);
                    data_ev.trial_type(ind_unexpNoUS) = {'noUSpostCSplus_first3'};
                    % last three US are unexpected stimulations
                    ind_unexpUS = find(strcmp(data_ev.trial_type, 'USpostCSplus'));
                    ind_unexpUS = ind_unexpUS(end-2:end);
                    data_ev.trial_type(ind_unexpUS) = {'USpostCSplus_last3'};
                    
                    % early and late reacquisition
                    % 18 blocks, 6+5+4+3*1 = 18, note that the last CS-
                    % might arguably be changed!
                    ind_CSplus = find(strcmp(data_ev.trial_type, 'CSplus'));
                    ind_CSminus = find(strcmp(data_ev.trial_type, 'CSminus'));
                    data_ev.trial_type(ind_CSminus(1:9)) = {'CSminus_reacq_early'};
                    data_ev.trial_type(ind_CSplus(1:9)) = {'CSplus_reacq_early'};
                    data_ev.trial_type(ind_CSminus(10:18)) = {'CSminus_reacq_late'};
                    data_ev.trial_type(ind_CSplus(10:18)) = {'CSplus_reacq_late'};
                    
                    % early and late reextinction (14 blocks) and early late "reinstatement" (13 blocks with 3 unexpected US)
                    % again for reextinction not that last CS- might
                    % arguably be changed!(after unexpected US)
                    data_ev.trial_type(ind_CSminus(19:25)) = {'CSminus_reext_early'};
                    data_ev.trial_type(ind_CSplus(19:25)) = {'CSplus_reext_early'};
                    data_ev.trial_type(ind_CSminus(26:32)) = {'CSminus_reext_late'};
                    data_ev.trial_type(ind_CSplus(26:32)) = {'CSplus_reext_late'};
                    data_ev.trial_type(ind_CSminus(32:37)) = {'CSminus_reins_early'};
                    data_ev.trial_type(ind_CSplus(32:37)) = {'CSplus_reins_early'};
                    data_ev.trial_type(ind_CSminus(38:44)) = {'CSminus_reins_late'};
                    data_ev.trial_type(ind_CSplus(38:44)) = {'CSplus_reins_late'};
                    
%                     % possibly divide USes and noUSes?
%                     ind_USpostCSplus = find(strcmp(data_ev.trial_type, 'USpostCSplus'));
%                     ind_noUSpostCSplus = find(strcmp(data_ev.trial_type, 'noUSpostCSplus'));
%                     ind_noUSpostCSminus = find(strcmp(data_ev.trial_type, 'noUSpostCSminus'));
                end
                
                % Convert event table to multcond
                [names, ~, uniq_ind] = unique(data_ev.trial_type');
                onsets = cell(1, length(names));
                durations = cell(1, length(names));
                for iN = 1:length(names)
                    onsets{iN} = data_ev.onset(uniq_ind == iN);
                    durations{iN} = zeros(sum(uniq_ind == iN), 1); % all durations set to 0 for event based analysis
                end
                
                % split US name into late US and early US
                if ~ismember(fp_func_phase{iFUNC}, {'Volatile'})
                    condition = unique(data_ev.trial_type);
%                     condition = {'CSplus', 'CSminus', 'USpost'};
                    for iCon = 1:length(condition)
                        logical_Condition = ismember(names, condition{iCon});
                        if any(logical_Condition)
                            % get onsets and durations
                            onsets_con = onsets{logical_Condition};
                            durations_con = durations{logical_Condition};

                            % remove old condition
                            names(logical_Condition) = [];
                            onsets(logical_Condition) = [];
                            durations(logical_Condition) = [];

                            % add new early and late conditions
                            names = [names, {[condition{iCon}, '_early']}, ...
                                            {[condition{iCon}, '_late']}];
                            onsets = [onsets, {onsets_con(onsets_con < median(onsets_con))}];
                            onsets = [onsets, {onsets_con(onsets_con >=median(onsets_con))}];
                            durations = [durations, {durations_con(onsets_con < median(onsets_con))}];
                            durations = [durations, {durations_con(onsets_con >=median(onsets_con))}];
                        end
                    end
                end
                
                % save final multcond
                save(fp_multcond, 'names', 'onsets', 'durations')
            end

            %% 1stLevel - task func
    %         fl_multcond = func_dirl2fl(dir(fullfile(fp_1stLevel, '*multcond.mat')));
    %         fp_multcond = func_matchFiles(fl_multcond, fp_func);

            if isempty(fp_multcond)
                doStats = false;
            else
                doStats = true;
            end
            if doStats
                fp_1stLevel_singleFunc = fullfile(fp_1stLevel, ['stats_task-', fp_func_parsed{iFUNC}.task, '_', ...
                    'ses-', fp_func_parsed{iFUNC}.ses, '_', ...
                    'run-', fp_func_parsed{iFUNC}.run, '_', ...
                    'acq-', fp_func_parsed{iFUNC}.acq]);
                
                % define contrastTable here!
                load(fp_multcond, 'names')
                weights = cell(length(names), 1);
                for i=1:length(names)
                    weights{i} = zeros(1, length(names));
                    weights{i}(i) = 1;
                end
                contrastTable = cell2table([names', weights], 'VariableNames', {'names', 'weights'});
                
                %% Adding contrast to contrastTable 
                %% CSplus and CSminus contrast
                % early vs late
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'CSplus_early_gt_late', 'CSplus_late_gt_early'}, ...
                    {'CSplus_early'}, {'CSplus_late'});

                % plus vs minus for different periods
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'CSplus_early_gt_CSminus_early', 'CSminus_early_gt_CSplus_early'}, ...
                    {'CSplus_early'}, {'CSminus_early'});

                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'CSplus_late_gt_CSminus_late', 'CSminus_late_gt_CSplus_late'}, ...
                    {'CSplus_late'}, {'CSminus_late'});

                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'CSplus_gt_CSminus', 'CSminus_gt_CSplus'}, ...
                    {'CSplus_early', 'CSplus_late'}, ...
                    {'CSminus_early', 'CSminus_late'});

                %% habituation of US
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'USpostCSplus_early_gt_late', 'USpostCSplus_late_gt_early'}, ...
                    {'USpostCSplus_early'}, {'USpostCSplus_late'});

                %% sensation of US
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'USsensation', 'USsensation_reverse'}, ...
                    {'USpostCSplus_early', 'USpostCSplus_late'}, ...
                    {'noUSpostCSminus_early', 'noUSpostCSminus_late'});

                %% Volatile contrasts
                % unexpected omission of US
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'unexpectedOmission', 'unexpectedOmission_reverse'}, ...
                    {'noUSpostCSplus_first3'}, ...
                    {'noUSpostCSminus'});
                % unexpected application of US
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'unexpectedUS', 'unexpectedUS_reverse'}, ...
                    {'USpostCSplus_last3'}, ...
                    {'USpostCSplus'});
                % sensation of US
                contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                    {'USsensation', 'USsensation_reverse'}, ...
                    {'USpostCSplus', 'USpostCSplus_last3'}, ...
                    {'noUSpostCSminus'});
                % CS+ and CS-
                for volPhase = {'reacq', 'reext', 'reins'}
                    % early vs late
                    contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                        {['CSplus_', volPhase{1}, '_early_gt_late'], ...
                         ['CSplus_', volPhase{1}, '_late_gt_early']}, ...
                        {['CSplus_', volPhase{1}, '_early']}, ...
                        {['CSplus_', volPhase{1}, '_late']});
    
                    % plus vs minus for different periods
                    contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                        {['CSplus_', volPhase{1}, '_early_gt_CSminus_early'], ...
                         ['CSminus_', volPhase{1}, '_early_gt_CSplus_early']}, ...
                        {['CSplus_', volPhase{1}, '_early']}, ...
                        {['CSminus_', volPhase{1}, '_early']});
    
                    contrastTable = contrastTable_addPosNegCon(contrastTable, ...
                        {['CSplus_', volPhase{1}, '_late_gt_CSminus_late'], ...
                         ['CSminus_', volPhase{1}, '_late_gt_CSplus_late']}, ...
                        {['CSplus_', volPhase{1}, '_late']}, ...
                        {['CSminus_', volPhase{1}, '_late']});
                end

                % if CSplus_late and CSminus_late exist, create CSplus_late>CSminus_late

                % if CSplus_early+late and CSminus_early+late exist, create CSplus>CSminus

                % CS+>CS-
                % CS+<CS-
                % USpostCSplus>noUSpostCSminus (sensation of US: expected receival of US vs expected omission of US reveals)
                % USpostCSplus<noUSpostCSminus (sensation of the omission of US: expected omission of US vs expected receival of US reveals)
                %
                % divide CS+>CS-in early and late
                % divide USpostCSplus>noUSpostCSminus in early and late
                % 
                % VOLATILE ONLY
                % noUSpostCSplus>noUSpostCSminus (unexpected omission of US)
                % noUSpostCSplus<noUSpostCSminus (? expected omission of US)
                %
                % 2*3* + 2*1*1
                %
                % add first three extinction, first three recall, first
                % three acquisition and first three volatile reacquisition
                % and first three volatile reextinction
                %
                % volatile US events early and late division (excluding unexpected US events)
                %
                % Possibly look at CS events after unexpected event,
                % attentional effects?

                func_singleFuncStats(fl_func_norm{iFUNC}, fp_1stLevel_singleFunc, fp_multcond, fp_reg{iFUNC}, smoothing_FWHM, contrastTable)
                fl_spmT = func_dirl2fl(dir(fullfile(fp_1stLevel_singleFunc, 'spmT*.nii*')));
                fl_spmT_SUITflatmapData = cell(length(fl_spmT), 1);
                for iSPM = 1:length(fl_spmT)
                    % Right Cerebellum gray, cobra atlas, multiplication and
                    % nonzero mean of spmT
                    Vo = func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'RCEg_', '', 'nii');
                    if ~exist(Vo, 'file')
                        spm_imcalc({func_changeFileName_PrePostFileext(fl_spmT{iSPM}, '', '', 'nii'), ...
                            fullfile(fp_packages, 'spm12', 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', 'cobra.nii')}, ...
                            Vo, ...
                            'i1.*(i2==115)')
                        Vo = niftiread(Vo);
                        LCEg_spmT_mean(iFUNC) = mean(nonzeros(Vo));
                    end

                    % Left Cerebellum gray, cobra atlas, multiplication and
                    % nonzero mean of spmT
                    Vo = func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'LCEg_', '', 'nii');
                    if ~exist(Vo, 'file')
                        spm_imcalc({func_changeFileName_PrePostFileext(fl_spmT{iSPM}, '', '', 'nii'), ...
                            fullfile(fp_packages, 'spm12', 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', 'cobra.nii')}, ...
                            Vo, ...
                            'i1.*(i2==15)')
                        Vo = niftiread(Vo);
                        RCEg_spmT_mean(iFUNC) = mean(nonzeros(Vo));
                    end
                    
                    % Show spmT on cerebellar flatmap
                    flatmap_tempFig_path = fullfile(fp_fig, sprintf('SUITflatmap_iSPM-%d.pdf', iSPM));
                    if ~exist(flatmap_tempFig_path, 'file')
                        flatmap_tempFig = figure('Units', 'normalized', 'Position', [0 0 1 1]);
                        fl_spmT_SUITflatmapData{iSPM} = suit_map2surf(fl_spmT{iSPM}, ...
                            'space', 'SUIT', 'stats', @nanmean);
                        suit_plotflatmap(fl_spmT_SUITflatmapData{iSPM})
                        % suit_plotflatmap(Data, 'threshold', 1, 'cscale', [1, 2])
                        saveas(flatmap_tempFig, flatmap_tempFig_path)
                        close(flatmap_tempFig)
                    end
                    
                    % Show spmT on cerebellar flatmap with SUIT
                    % normalization
                    
                    % Transform spmT back to native coregistered space with the inverse CAT12 transformation
                    
                    % 
                end
            end
        
        end
        
                
        % Right Cerebellum gray, cobra atlas, multiplication and
        % nonzero mean of tSNR
        Vo = func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'RCEg_tsnr_', '', 'nii');
        if ~exist(Vo, 'file')
            spm_imcalc({func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'tsnr_', '', 'nii'), ...
                fullfile(fp_packages, 'spm12', 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', 'cobra.nii')}, ...
                Vo, ...
                'i1.*(i2==115)')
            Vo = niftiread(Vo);
            RCEg_tsnr_mean(iFUNC) = mean(nonzeros(Vo));
        end
        
        % Left Cerebellum gray, cobra atlas, multiplication and
        % nonzero mean of tSNR
            Vo = func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'LCEg_tsnr_', '', 'nii');
        if ~exist(Vo, 'file')
            spm_imcalc({func_changeFileName_PrePostFileext(fl_func_norm{iFUNC}, 'tsnr_', '', 'nii'), ...
                fullfile(fp_packages, 'spm12', 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', 'cobra.nii')}, ...
                Vo, ...
                'i1.*(i2==15)')
            Vo = niftiread(Vo);
            LCEg_tsnr_mean(iFUNC) = mean(nonzeros(Vo));
        end
        
        %% Save workspace for each func (right now you throw some stuff out when looping func) 
        save(fullfile(fp_pl, ['workspace_iFUNC-', num2str(iFUNC)]))
    
    end
    
    % after finishing iFUNC loop you can construct subses table
    table_subses{iSS} = table(fn_func, ...
            LCEg_spmT_mean, RCEg_spmT_mean, ...
            LCEg_tsnr_mean, RCEg_tsnr_mean, ...
        'VariableNames', {'fn_func', ...
            'LCEg_spmT_mean', 'RCEg_spmT_mean', ...
            'LCEg_tsnr_mean', 'RCEg_tsnr_mean'});
        
end

% after finishing subses loop you can construct derivative folder table
table_derivatives = vertcat(table_subses{:});

%% 2nd level analysis for each contrast for each phase
% Acquisition_USpostCSplusEarly_mainEffect
fl_stats_AcqUSearly = func_dirl2fl(dir(fullfile(fp_deri, ...
    '*', '1stLevel_EB_4p5', 'task-fear_acq-stxtr1620_run-2', 'con_0005.nii')));
% fl_stats_AcqUSearly = func_dirl2fl(dir(fullfile(fp_deri, '*', '1stLevel_EB_4p5', 'task-fear_ses-1_run-2_acq-stxtr1620', 'con_0005.nii')));
func_2ndLevelStats(fl_stats_AcqUSearly, fullfile(fp_2ndLevel, ...
    'Acquisition_USpostCSplusEarly_mainEffect'));

% Volatile_noUSpostCSplusFirst3_mainEffect
func_2ndLevelStats(fl_stats_AcqUSearly, fullfile(fp_2ndLevel, ...
    'Acquisition_USpostCSplusEarly_mainEffect'));

% % Volatile_USpostCSplusLast3_mainEffect
% func_2ndLevelStats(fl_stats_AcqUSearly, fullfile(fp_2ndLevel, ...
%     'Acquisition_USpostCSplusEarly_mainEffect'));

function contrastTable = contrastTable_addPosNegCon(contrastTable, contrastNames, posNames, negNames)
    if ~iscell(contrastNames)
        contrastNames = {contrastNames};
    end
    if all(ismember([posNames, negNames], contrastTable.names))
        contrastTable_addition = cell2table([contrastNames(1), ...
            {ismember(contrastTable.names, posNames)' - ... 
            ismember(contrastTable.names, negNames)'}], ...
            'VariableNames', {'names', 'weights'});
            contrastTable_addition.weights = contrastTable_addition.weights(1:size(contrastTable.weights, 2));
            contrastTable = [contrastTable; contrastTable_addition];
        if length(contrastNames)==2
            contrastTable_addition = cell2table([contrastNames(2), ...
                {-ismember(contrastTable.names, posNames)' + ... 
                ismember(contrastTable.names, negNames)'}], ...
                'VariableNames', {'names', 'weights'});
            contrastTable_addition.weights = contrastTable_addition.weights(1:size(contrastTable.weights, 2));
            contrastTable = [contrastTable; contrastTable_addition];
        end
    end
end
