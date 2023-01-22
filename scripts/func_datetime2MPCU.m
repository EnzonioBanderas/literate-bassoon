function MPCU = func_datetime2MPCU(input_datetime)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
% MPCU = milliseconds(input_datetime - datetime('today'));
MPCU = milliseconds(timeofday(input_datetime));

end

