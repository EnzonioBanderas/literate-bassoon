function mr_syn = func_readMRItriggers(fp_bids)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

fl_nii = func_dirl2fl(dir(fullfile(fp_bids, '*', '*.nii.gz')));
if ~isempty(fl_nii)
fprintf('MRI syn structure definition\n')

% get MRI files lists
fl_json = cell(length(fl_nii), 1);
for iFL = length(fl_nii):-1:1
    fl_json{iFL} = func_changeFileName_PrePostFileext(fl_nii{iFL}, '', '', 'json');
    
    jinfo = jsondecode(fileread(fl_json{iFL}));
    startTime_nii = func_DICOMtime2MPCU(jinfo.AcquisitionTime);
    try
        info_nii = niftiinfo(fl_nii{iFL});
    catch
        fp_nii_unzipped = gunzip(fl_nii{iFL});
        info_nii = niftiinfo(fp_nii_unzipped{1});
        delete(fp_nii_unzipped{1});
    end
    
    if length(info_nii.ImageSize)<4
        nT_nii = 1;
    else
        nT_nii = info_nii.ImageSize(4);
    end
    if isfield(jinfo, 'RepetitionTime')
        pulse_timings = (startTime_nii:jinfo.RepetitionTime*1e3:...
            startTime_nii + jinfo.RepetitionTime*1e3*(nT_nii-1))';
        pulse_duration = NaN(length(pulse_timings), 1);
%         pulse_timings = startTime_nii;
%         pulse_duration = NaN;
    else
        pulse_timings = startTime_nii;
        pulse_duration = NaN;
    end
    
    pulse_index = (1:length(pulse_timings))';
    
    pulse_table = table(pulse_index, pulse_timings, pulse_duration, ...
        'VariableNames', {'index', 'time', 'duration'});
    
    mr_syn(iFL).fp = fl_nii{iFL};
    mr_syn(iFL).pulse_table = pulse_table;
    mr_syn(iFL).data_modality = 'MRI';
end

% Match files to each other
% [fp_bold, fp_eye] = func_matchFiles(fl_bold, fl_eyetracker_table(1));

end

end

