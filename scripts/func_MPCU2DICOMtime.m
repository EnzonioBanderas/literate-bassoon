function DICOMtime = func_MPCU2DICOMtime(MPCU)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if isstr(MPCU)
    MPCU = str2double(MPCU);
end

% MPCU_hours = floor(MPCU/(60*60*1000));
% MPCU_remaining = MPCU - MPCU_hours*60*60*1000;
% MPCU_minutes = floor(MPCU_remaining/(60*1000));
% MPCU_remaining = MPCU_remaining - MPCU_minutes*60*1000;
% MPCU_seconds = floor(MPCU_remaining/(1000));
% MPCU_remaining = MPCU_remaining - MPCU_seconds*1000;
% MPCU_miliseconds = MPCU_remaining;
% 
% DICOMtime = [num2str(MPCU_hours, '%02.f'), ':', ...
%     num2str(MPCU_minutes, '%02.f'), ':', ...
%     num2str(MPCU_seconds, '%02.f'), '.', ...
%     num2str(MPCU_miliseconds), ...
%     repmat('0', [1, 6-length(num2str(MPCU_miliseconds))])];

DICOMtime = datestr(func_MPCU2datetime(MPCU), 'HH:MM:ss.FFF');

end

