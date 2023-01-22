function [output_table, output_table_all] = func_getFirst_inSequence(input_table)
%func_getFirst_inSequence: Convert a table with start points and durations
% to a table with only the first points of each sequence 
%
% uses median now, which assumes that the number of distance between
% points is larger than 
output_table = input_table;
output_table(:, 3) = [inf; diff(input_table(:, 1))] - input_table(:, 2); % calculate intersignal lengths
output_table(:, 4) = output_table(:,3) > median(output_table(:,3))*3; 

% calculate number of signals associated with first sequence signal
firstSeq_ind = [find(output_table(:, 4)); length(output_table(:, 4))+1];
for i=1:length(firstSeq_ind)-1
    output_table(firstSeq_ind(i), 5) = firstSeq_ind(i+1) - firstSeq_ind(i);
end

% copy before filtering
output_table_all = output_table;

% only keep signals with high intersignals lengths (the first signals in their respective sequences)
output_table = output_table(output_table(:, 4) > 0, [1, 2, 5]);
end

