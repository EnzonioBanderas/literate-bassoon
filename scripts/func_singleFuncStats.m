function func_singleFuncStats(fp_func_norm, fp_stats, fp_multcond, fp_reg, fwhm, contrastTable)
%FUNC_SINGLEFUNCSTATS Do statistics on single functional file
%   smoothing, model specification, model estimation, contrast definition
%   (just all possible), tSNR computation and figure saving
%
%   TD
%       1) Add contrast parameter input [0 0 0 1 0 0 0 1; 0 1 0 1 0 1 0 1 ...]
%

if ~exist('fwhm', 'var')
    fwhm = 4.5;
end

fl_func_norm = func_fp2fl_Texpand(fp_func_norm);
fp_func_norm_info = niftiinfo(fp_func_norm);
TR = fp_func_norm_info.PixelDimensions(4);

fp_packages = func_getPackagesFolder();

%% smoothing 1
matlabbatch{1}.spm.spatial.smooth.data = fl_func_norm;
matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = ['s', strrep(num2str(fwhm), '.', 'p')];

%% model specification 2
% '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-16/rawdata/sub-z7t7095/func/topup/ses-1'
if ~exist(fp_stats, 'dir')
    mkdir(fp_stats)
%     rmdir(fp_stats, 's')
% end

    % matlabbatch{2}.spm.stats.fmri_spec.dir = {fp_stats};
    % matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
    % matlabbatch{2}.spm.stats.fmri_spec.timing.RT = TR;
    % % matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 1.890;
    % matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t = 16;
    % matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    % matlabbatch{2}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    % % matlabbatch{2}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    % % '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-16/rawdata/sub-z7t7095/func/topup/ses-1/sub-z7t7095_20220316T151607_task-fear_ses-1_multcond.mat'
    % matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {fp_multcond};
    % matlabbatch{2}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    % % '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/2022-03-16/rawdata/sub-z7t7095/func/topup/rp_sub-z7t7095_task-fear1_acq-ep3d1p5mmte20tr1620_dir-PA_bold.txt'
    % matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = {fp_reg};
    % matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;
    % matlabbatch{2}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    % matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    % matlabbatch{2}.spm.stats.fmri_spec.volt = 1;
    % matlabbatch{2}.spm.stats.fmri_spec.global = 'None';
    % matlabbatch{2}.spm.stats.fmri_spec.mthresh = 0.1;
    % matlabbatch{2}.spm.stats.fmri_spec.mask = {''};
    % matlabbatch{2}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % matlabbatch{2}.spm.stats.fmri_spec.dir = {'/media/diskEvaluation/Evaluation/sfb1280a05study7/derivatives/prelimUSfMRI_test_sub-Z7T7186_ses-1/sub-Z7T7186/ses-1/1stLevel/task-fear_acq-stxtr1620_run-1'};
    matlabbatch{2}.spm.stats.fmri_spec.dir = {fp_stats};
    matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
    % matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 1.62000000476837;
    matlabbatch{2}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{2}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{2}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    % matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {'/media/diskEvaluation/Evaluation/sfb1280a05study7/derivatives/prelimUSfMRI_test_sub-Z7T7186_ses-1/sub-Z7T7186/ses-1/1stLevel/sub-Z7T7186_ses-1_task-fear_acq-stxtr1620_run-1_events_multcond.mat'};
    matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {fp_multcond};
    matlabbatch{2}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    % matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = {'/media/diskEvaluation/Evaluation/sfb1280a05study7/derivatives/prelimUSfMRI_test_sub-Z7T7186_ses-1/sub-Z7T7186/ses-1/func/rp_sub-Z7T7186_ses-1_task-fear_acq-stxtr1620_run-1_bold.txt'};
    matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = {fp_reg};
    matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{2}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{2}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{2}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{2}.spm.stats.fmri_spec.mthresh = 0.1;
    matlabbatch{2}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{2}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    
    % matlabbatch{6}.spm.stats.fmri_spec.dir = {'/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/pilots/derivatives/feartapping/sub-Z7T7146/ses-220506/stats/task-feartapping_acq-ptxtr1890_run-1'};
    % matlabbatch{6}.spm.stats.fmri_spec.timing.units = 'secs';
    % matlabbatch{6}.spm.stats.fmri_spec.timing.RT = 1.89;
    % matlabbatch{6}.spm.stats.fmri_spec.timing.fmri_t = 16;
    % matlabbatch{6}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    % matlabbatch{6}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    % matlabbatch{6}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    % matlabbatch{6}.spm.stats.fmri_spec.sess.multi = {''};
    % matlabbatch{6}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    % matlabbatch{6}.spm.stats.fmri_spec.sess.multi_reg = {''};
    % matlabbatch{6}.spm.stats.fmri_spec.sess.hpf = 128;
    % matlabbatch{6}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    % matlabbatch{6}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    % matlabbatch{6}.spm.stats.fmri_spec.volt = 1;
    % matlabbatch{6}.spm.stats.fmri_spec.global = 'None';
    % matlabbatch{6}.spm.stats.fmri_spec.mthresh = 0.8;
    % matlabbatch{6}.spm.stats.fmri_spec.mask = {''};
    % matlabbatch{6}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% model estimation 3
    matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
    
    %% contrast 4
    matlabbatch{4}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    % make contrasts dependent on matlab cond or use input contrast table
    if ~exist('contrastTable', 'var')
        load(fp_multcond, 'names')
        for i=1:length(names)
            tcon_weights = zeros(1, length(names));
            tcon_weights(i) = 1;
            matlabbatch{4}.spm.stats.con.consess{i}.tcon.name = names{i};
            matlabbatch{4}.spm.stats.con.consess{i}.tcon.weights = tcon_weights;
            matlabbatch{4}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
        end
    else % use defined contrastTable
        for i=1:size(contrastTable, 1)
            matlabbatch{4}.spm.stats.con.consess{i}.tcon.name = contrastTable.names{i};
            matlabbatch{4}.spm.stats.con.consess{i}.tcon.weights = contrastTable.weights(i,:);
            matlabbatch{4}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
        end
    end
    matlabbatch{4}.spm.stats.con.delete = 0;
    
    %% ortho check skip
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.spmmat = {fullfile(fp_SPM, fn_SPM)};
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.use_unsmoothed_data = 1;
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.adjust_data = 1;
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.outdir = {''};
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.fname = 'CATcheckdesign_';
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.check_SPM_cov.do_check_cov.save = 0;
    % matlabbatch{4}.spm.tools.cat.tools.check_SPM.check_SPM_ortho = 1;
    
    %%
    matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    matlabbatch{5}.spm.stats.results.conspec(1).titlestr = '';
    matlabbatch{5}.spm.stats.results.conspec(1).contrasts = Inf;
    matlabbatch{5}.spm.stats.results.conspec(1).threshdesc = 'FWE';
    matlabbatch{5}.spm.stats.results.conspec(1).thresh = 0.05;
    matlabbatch{5}.spm.stats.results.conspec(1).extent = 0;
    matlabbatch{5}.spm.stats.results.conspec(1).conjunction = 1;
    matlabbatch{5}.spm.stats.results.conspec(1).mask.image.name = func_fp2fl_Texpand(fullfile(fp_packages, 'spm12', 'tpm', 'mask_ICV.nii'));
    % matlabbatch{5}.spm.stats.results.conspec(1).mask.image.name = {'/media/diskEvaluation/Evaluation/sfb1280a05study2b/scripts/packages/spm12/tpm/mask_ICV.nii,1'};
    matlabbatch{5}.spm.stats.results.conspec(1).mask.image.mtype = 0;
    
    % matlabbatch{5}.spm.stats.results.conspec(1).titlestr = '';
    % matlabbatch{5}.spm.stats.results.conspec(1).contrasts = Inf;
    % matlabbatch{5}.spm.stats.results.conspec(1).threshdesc = 'none';
    % matlabbatch{5}.spm.stats.results.conspec(1).thresh = 0.001;
    % matlabbatch{5}.spm.stats.results.conspec(1).extent = 0;
    % matlabbatch{5}.spm.stats.results.conspec(1).conjunction = 1;
    % matlabbatch{5}.spm.stats.results.conspec(1).mask.image.name = func_fp2fl_Texpand(fullfile(fp_packages, 'spm12', 'tpm', 'mask_ICV.nii'));
    % % matlabbatch{5}.spm.stats.results.conspec(1).mask.image.name = {'/media/diskEvaluation/Evaluation/sfb1280a05study2b/scripts/packages/spm12/tpm/mask_ICV.nii,1'};
    % matlabbatch{5}.spm.stats.results.conspec(1).mask.image.mtype = 0;
    
    matlabbatch{5}.spm.stats.results.units = 1;
    matlabbatch{5}.spm.stats.results.export{1}.pdf = true;
    
    %% Save matlabbatch in processing logs folder
    [fp_func, ~, ~] = fileparts(fp_func_norm);
    fp_pl = fullfile(fp_func, '..', 'processing_logs');
    mkdir(fp_pl)
    cd(fp_pl)
    save(fullfile(fp_pl, ...
        ['matlabbatch_', ...
        'sfStats_', ...
        't-', datestr(now, 'yymmddTHHMMSS'), ...
        '.mat']), 'matlabbatch')
    
    %% Run
    spm_jobman('run', matlabbatch)
else
    warning('1stLevel stats directory already exists')
end

end
