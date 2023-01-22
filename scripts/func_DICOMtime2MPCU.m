function MPCU = func_DICOMtime2MPCU(DICOMtime)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% if ~isstr(DICOMtime)
%     DICOMtime = num2str(DICOMtime);
% end
% MPCU = str2double(DICOMtime(1:2))*60*60*1000 + ...
%        str2double(DICOMtime(3:4))*60*1000 + ...
%        str2double(DICOMtime(5:6))*1000 + ...
%        str2double(DICOMtime(7:8));

MPCU = milliseconds(timeofday(datetime(DICOMtime, 'InputFormat', 'HH:mm:ss.SSSSSS')));

end

