function [mea_path, mes_path, met_path, time_actual, time_estimate, ...
    output_struct] = ...
arringtonEyetracker2bids(fp_eyetracker_table, fp_output)
% =========================================================================
% Parse Arrington Eyetracker output data into 3 BIDS compatible files
% =========================================================================
%   Parse Arrington Eyetracker input table into three separate files:
%   1) One data table (measurement.tsv)
%   2) One message table (messages.tsv)
%   3) One header metadata file (metadata.json)
%
%   With as input: 
%   1) Path to eyetracker file output 
%   2) Optional path to output folder of the 3 output files
%   
%   Outputs the 3 file path strings in the order of measurement, messages 
%   and metadata.
% ========================================================================= 

% If optional output folder path input parameter is not given default to
% current folder
%
% Remaining issues:
%   1) Current warning condition does not make sense, instead of looking at
%       the difference between time and estimates it should already be
%       sufficient to set a threshold for succesive time differences 
%       (at for example 5 times the sampling time)
%   2) If the warning condition is fulfilled, the time estimate no longer
%       makes sense and instead a time estimate time point array should be
%       produced using start and end times, after which pupillometry data is
%       interpolated for missing values. Preferably the values will not be
%       missing at all, but this would be a way to deal with the problem for
%       the many incomplete eyetracker datasets. Afterwards, timing of events
%       should be checked.
%
if ~exist('fp_output','var')
    fp_output = pwd;
else
    if ~exist(fp_output, 'dir')
        mkdir(fp_output)
    end
end

% Separate header from data
unparsed_table = fileread(fp_eyetracker_table);
unparsed_table = split(unparsed_table, '--------------------------------------------------------');
header_str = unparsed_table(1:4);
data_str = unparsed_table{5};

% Separate data into messages and measurements
data_str = splitlines(data_str);
nLines = length(data_str);
nCol = zeros(1, nLines);
for iLine = 1:nLines
% for iLine = 1:1
    data_str{iLine} = split(data_str{iLine});
    nCol(iLine) = length(data_str{iLine});
end
measurement_str = data_str(nCol==28);
measurement_colNames1 = measurement_str{1}(2:end);
measurement_colNames2 = measurement_str{2}(2:end);
measurement_colNames = cell(1, length(measurement_colNames2));

for iCol = 1:length(measurement_colNames2)
    if measurement_colNames1{iCol}(1)=='A'
        measurement_colNames{iCol} = ['EyeA-', measurement_colNames2{iCol}];
    elseif measurement_colNames1{iCol}(1)=='B'
        measurement_colNames{iCol} = ['EyeB-', measurement_colNames2{iCol}];
    else
        measurement_colNames(iCol) = measurement_colNames2(iCol);
    end
end
measurement_table = [measurement_str{3:end}]';
measurement_table = measurement_table(:, 2:end); 
messages_str = data_str(nCol==3);
messages_table = [messages_str{:}]';
messages_table = [{'TotalTime', 'Message'}; messages_table(:, 2:end)];
framerate_str = data_str(nCol==5);
freq_sampling1 = str2double(framerate_str{1}{5});
freq_sampling2 = str2double(framerate_str{2}{5});
if freq_sampling1 ~= freq_sampling2
    error('Sampling times of different eyes are not equal')
end

% Parse header
header_values = [header_str{[1,2,4]}];
header_values = splitlines(header_values);
header_values(7) = []; header_values(end-5) = []; header_values(end) = [];
nLines = length(header_values);
header_values_table = cell(nLines, 2);
for iLine = 1:nLines % Remove starting integer at each line
    header_values{iLine} = split(header_values{iLine});
    header_values{iLine} = strjoin(header_values{iLine}(2:end));
    if iLine < 7 % dubbel dot parsing
        firstI = find(header_values{iLine} == ':', 1, 'first');
        header_values_table{iLine, 1} = strrep(header_values{iLine}(1:firstI-1), ' ', '');
        header_values_table{iLine, 1} = strrep(header_values_table{iLine, 1}, ' ', '');
        header_values_table{iLine, 1} = strrep(header_values_table{iLine, 1}, ':', '');
        header_values_table{iLine, 2} = header_values{iLine}(firstI+1:end);
    else % space parsing
        header_values{iLine} = split(header_values{iLine});
        header_values_table(iLine, 1) = header_values{iLine}(1);
        header_values_table{iLine, 1} = strrep(header_values_table{iLine, 1}, ' ', '');
        header_values_table{iLine, 1} = strrep(header_values_table{iLine, 1}, ':', '');
        header_values_table{iLine, 2} = strjoin(header_values{iLine}(2:end));
    end
end
header_dataNotes = splitlines(header_str{3});
header_dataNotes = header_dataNotes(3:end-1);
nLines = length(header_dataNotes);
for iLine = 1:nLines % Remove starting integer at each line
    header_dataNotes{iLine} = header_dataNotes{iLine}(4:end);
end
header_table = [header_values_table; {'DataNotes', header_dataNotes}];

% Add starting and sampling time
time_actualA = str2double(measurement_table(:,1));
time_actualB = str2double(measurement_table(:,12));
time_maxDiff_eyes = max(time_actualA - time_actualB);
if time_maxDiff_eyes < 1
    time_actual = mean([time_actualA, time_actualB], 2);
else
    warning(['Maximum time difference between EyeA and EyeB ' ...
             'is greater than 1 ms.'])
    time_actual = time_actualA;
end
% time_actual = mean(str2double(measurement_table(2:end,[1, 12])), 2);
% time_sampling = mean(diff(time_actual));
time_sampling = 1 / freq_sampling1;
time_estimate = time_sampling * (0:(length(time_actual)-1)) + time_actual(1);
time_maxDiff = max(abs(time_actual - time_estimate')) * 1e3; % max ms time difference
% freq_samplingDiff = abs(freq_sampling1 - 1/time_sampling);
if time_maxDiff > 1
    warning(['Maximum time difference between estimate and time given ' ...
             'by Eyetracker software is greater than 1 ms.'])    
         figure; plot(time_actual(1:end-1), diff(time_actual))
    [~, fn_eyetracker_table, ~] = fileparts(fp_eyetracker_table);
    title(fn_eyetracker_table, 'Interpreter', 'none')
    xlabel('time (s)')
    ylabel('successive time difference (s)')
end

% remove columns which are unnecessary from measurement table
measurement_table_copy = measurement_table;
measurement_colNames_copy = measurement_colNames;
measurement_table(:, [1,2,12,13]) = [];
measurement_colNames([1,2,12,13]) = [];

% Add columns, sampling freq and startime information to header
header_table_met = [{'Columns', measurement_colNames}; header_table];
header_table_met = [{'SamplingFrequency', num2str(freq_sampling1)}; header_table_met];
header_table_met = [{'StartTime', num2str(time_actual(1))}; header_table_met];

% Define output file names
[~, eyetracker_table_name, ~] = fileparts(fp_eyetracker_table);
mea_path = fullfile(fp_output, [eyetracker_table_name, '_eyetracker.tsv']);
mes_path = fullfile(fp_output, [eyetracker_table_name, '_eyetracker-events.tsv']);
met_path = fullfile(fp_output, [eyetracker_table_name, '_eyetracker.json']);

% Write measurement and messages to .tsv and .tsv.gz files
writecell(measurement_table, mea_path, 'FileType', 'text', 'Delimiter', '\t')
% gzip(mea_path); delete(mea_path); mea_path = [mea_path, '.gz'];
gzip(mea_path); mea_path = [mea_path, '.gz'];
writecell(messages_table, mes_path, 'FileType', 'text', 'Delimiter', '\t')
% gzip(mes_path); delete(mes_path);  mes_path = [mes_path, '.gz'];
gzip(mes_path); mes_path = [mes_path, '.gz'];

% Write header to (pretty) .json
header_table_met = cell2table(header_table_met(:,2)', 'VariableNames', header_table_met(:,1));
fid = fopen(met_path, 'w');
header_json = jsonencode(header_table_met);
header_json = strrep(header_json, ',"', sprintf(',\r"'));
header_json = strrep(header_json, '[{', sprintf('[\r{\r'));
header_json = strrep(header_json, '}]', sprintf('\r}\r]'));
fprintf(fid, header_json);
fclose(fid);

% get output struct
output_struct.measurement_table = cell2table(measurement_table_copy, 'VariableNames', matlab.lang.makeValidName(measurement_colNames_copy));
output_struct.measurement_table.time = time_actual;
mea_fields = fields(table2struct(output_struct.measurement_table));
for iMEA = 1:(length(mea_fields)-2)
    output_struct.measurement_table.(mea_fields{iMEA}) = ...
        cellfun(@(x) str2double(x), output_struct.measurement_table.(mea_fields{iMEA}), 'uni', 0);
    output_struct.measurement_table.(mea_fields{iMEA}) = [output_struct.measurement_table.(mea_fields{iMEA}){:}]';
end
output_struct.messages_table = cell2table(messages_table(2:end, :), 'VariableNames', messages_table(1, :));
output_struct.messages_table.TotalTime = cellfun(@(x) str2double(x), output_struct.messages_table.TotalTime, 'uni', 0);
output_struct.messages_table.TotalTime = [output_struct.messages_table.TotalTime{:}]';
output_struct.header_struct = table2struct(header_table_met);

end
