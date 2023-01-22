function output_datetime = func_DICOMtime2datetime(DICOMtime, input_date)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% if ~isstr(DICOMtime)
%     DICOMtime = num2str(DICOMtime);
% end

if ~exist('input_date', 'var')
    datetime(DICOMtime, 'InputFormat', 'HH:mm:ss.SSSSSS')
else
    datetime([input_date, ' ', DICOMtime], 'InputFormat', 'yy-MM-dd HH:mm:ss.SSSSSS')
end

end

