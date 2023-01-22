function [output_signal, output_time] = func_table2signal(input_table, sampling_time, output_time)
% Given a table with start points and durations, get output signal and
% output time, time assumed to be in ms
%
% does it work for [1, 0, 1]? Only 1 difference

if isempty(input_table)
end

% input_table = round(input_table);
if istable(input_table)
    input_table = table2array(input_table);
end

% add sampling time or frequency addition
if ~exist('sampling_time', 'var')
    sampling_time = 1;
end

% sort onsets
[~, sort_ind] = sort(input_table(:, 1));
input_table = input_table(sort_ind, :);

% get table times and values
table_times = zeros(4, size(input_table, 1));
table_times(1, :) = input_table(:, 1)' - 1; % start-1
table_times(2, :) = input_table(:, 1)'; % start
table_times(3, :) = (input_table(:, 1) + input_table(:, 2) - 1)'; % end
table_times(4, :) = (input_table(:, 1) + input_table(:, 2))'; % end + 1
table_times = table_times(:);
table_values = zeros(4, size(input_table, 1));
table_values(2, :) = 1;
table_values(3, :) = 1;
table_values = table_values(:);

% get time vector for interpolation
if ~exist('output_time', 'var')
    output_time = table_times(1) - sampling_time : sampling_time : ...
                  table_times(end) + sampling_time;
end

% use nearest neighbour interpolation to get output signal
output_signal = interp1(table_times, table_values, output_time, 'nearest', 0);

end
