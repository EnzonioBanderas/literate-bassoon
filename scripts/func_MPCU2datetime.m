function output_datetime = func_MPCU2datetime(MPCU, input_date)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
if exist('input_date', 'var')
    [YY, MM, dd] = datevec(input_date);
    output_datetime = datetime(YY, MM, dd)+milliseconds(MPCU);
else
    output_datetime = datetime('today')+milliseconds(MPCU);
end

end

