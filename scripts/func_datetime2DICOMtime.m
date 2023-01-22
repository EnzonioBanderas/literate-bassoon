function DICOMtime = func_datetime2DICOMtime(input_datetime)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% datestr(datetime('10:20:30.155', 'InputFormat', 'HH:mm:ss.SSS'), 'HH:mm:ss:FFF')
DICOMtime = [datestr(input_datetime, 'HH:MM:ss.FFF'), '000'];
end

