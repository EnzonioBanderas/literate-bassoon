function [min_time_diff_pulse, min_time_ind_pulse] = func_synTables_timeDiff(syn1, syn2, field)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('field', 'var')
        field = 'pulse_table';
    end

    syn1_merged = vertcat(syn1(:).(field));
%     syn_MRI_fp_merged = {syn1(:).fp}';
    syn2_merged = vertcat(syn2(:).(field));
    
    % what pulse in syn_SIEMENS is closest to what pulse in syn_MRI?
    min_time_diff_pulse = zeros(size(syn1_merged, 1), 1);
    min_time_ind_pulse = zeros(size(syn1_merged, 1), 1);
    for iM = 1:size(syn1_merged, 1)
        [min_time_diff_pulse(iM), min_time_ind_pulse(iM)] = ...
            min(abs(syn2_merged.time - syn1_merged.time(iM)));
    end
end

