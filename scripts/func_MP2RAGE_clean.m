function [fp_clean_cell, MP2RAGEimgRobustPhaseSensitive] = func_MP2RAGE_clean(fp_nii)
%FUNC_MP2RAGE_CLEAN Denoise MP2RAGE nii files in folder
%   You could change this so that it cleans all MP2RAGE nii files and not
%   just the first one listed by dir

fl_INV1 = dir(fullfile(fp_nii, '*inv-1_MP2RAGE.nii.gz'));
MP2RAGEimgRobustPhaseSensitive = cell(1, length(fl_INV1));
fp_clean_cell = cell(1, length(fl_INV1));
for i=1:length(fl_INV1)
    fp_INV1 = fullfile(fl_INV1(i).folder, fl_INV1(i).name);
    fp_INV1_nii = gunzip(fp_INV1);

    fp_INV2 = strrep(fp_INV1, 'inv-1_MP2RAGE', 'inv-2_MP2RAGE');
%     fp_INV2 = dir(fullfile(fp_nii, '*INV2.nii.gz'));
%     fp_INV2 = fullfile(fp_INV2(i).folder, fp_INV2(i).name);
    fp_INV2_nii = gunzip(fp_INV2);

    fp_UNI = strrep(fp_INV1, 'inv-1_MP2RAGE', 'MP2RAGE');
    if ~exist(fp_UNI, 'file')
        fp_UNI = strrep(fp_INV1, 'inv-1_MP2RAGE', 'T1w');
        if ~exist(fp_UNI, 'file')
            fp_UNI = strrep(fp_INV1, 'inv-1_MP2RAGE', 'UNIT1');
        end
    end
%     fp_UNI = dir(fullfile(fp_nii, '*UNI_Images.nii.gz'));
%     fp_UNI = fullfile(fp_UNI(i).folder, fp_UNI(i).name);
    fp_UNI_nii = gunzip(fp_UNI);

    fp_clean_nii = strrep(fp_INV1_nii, 'inv-1_MP2RAGE', 'clean');

    MP2RAGE.filenameINV1=fp_INV1_nii{1};
    MP2RAGE.filenameINV2=fp_INV2_nii{1};
    MP2RAGE.filenameUNI=fp_UNI_nii{1};
    MP2RAGE.filenameOUT=fp_clean_nii{1};
    
    fp_clean_cell{i} = fp_clean_nii{1};
    if ~exist(fp_clean_cell{i}, 'file')
        [MP2RAGEimgRobustPhaseSensitive{i}] = RobustCombination(MP2RAGE, 4);
    else
        warning([fp_clean_cell{i}, ' already exists, skipping MP2RAGE cleaning'])
    end

end

end

