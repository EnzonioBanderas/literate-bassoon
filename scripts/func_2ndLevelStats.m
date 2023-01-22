function func_2ndLevelStats(fl_con, fp_stats, stats_correction_cell, fp_reg)
%func_2ndLevelStats: Do 2ndLevel statistics on list of contrast files
%
% fp_stats = '/media/diskEvaluation/Evaluation/sfb1280a05study7/derivatives/prelimRerunTopupfMRI_2ndLevel';

if ~exist(fp_stats, 'dir')
    mkdir(fp_stats)

    matlabbatch{1}.spm.stats.factorial_design.dir = {fp_stats};
    %%
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fl_con;
    
    %%
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'main';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;
    
    matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{4}.spm.stats.results.conspec.contrasts = Inf;
    if ~exist(stats_correction_cell, 'var')
        matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
        matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
    else
        matlabbatch{4}.spm.stats.results.conspec.threshdesc = stats_correction_cell{1};
        matlabbatch{4}.spm.stats.results.conspec.thresh = stats_correction_cell{2};
    end
    matlabbatch{4}.spm.stats.results.conspec.extent = 0;
    matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{4}.spm.stats.results.conspec.mask.image.name = {'/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/packages/spm12/tpm/mask_ICV.nii,1'};
    matlabbatch{4}.spm.stats.results.conspec.mask.image.mtype = 0;
    matlabbatch{4}.spm.stats.results.units = 1;
    matlabbatch{4}.spm.stats.results.export{1}.pdf = true;
    
    %% Save matlabbatch in processing logs folder
    fp_stats_split = strsplit(fp_stats, filesep);
    fp_stats_2ndLevel = strjoin(fp_stats_split(1:end-1), filesep);
    fp_pl = fullfile(fp_stats_2ndLevel, 'processing_logs');
    mkdir(fp_pl)
    cd(fp_pl)
    save(fullfile(fp_pl, ...
        ['matlabbatch_', ...
        '2ndLevelStats_', ...
        't-', datestr(now, 'yymmddTHHMMSS'), ...
        '.mat']), 'matlabbatch')
    
    %% Run
    spm_jobman('run', matlabbatch)

else
    warning('2ndLevel stats subdirectory already exists')
end

end