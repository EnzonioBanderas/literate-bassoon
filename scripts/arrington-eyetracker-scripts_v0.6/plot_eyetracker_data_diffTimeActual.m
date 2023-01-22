% eyetracker_table_path = fullfile('eyetrackdata_test', '20211207_152540_test_100_realtimeHalfwayMovingwindow.txt');
% eyetracker_table_path = fullfile('eyetrackdata_test', '20211207_135237.txt');
eyetracker_table_paths = dir('eyetrack_data\*.txt');

nFile = length(eyetracker_table_paths);
for iFile = 1:nFile
    eyetracker_table_path = fullfile(eyetracker_table_paths(iFile).folder, eyetracker_table_paths(iFile).name);
    [a, b, c, time_actual, time_estimate] = arringtonEyetracker2bids(eyetracker_table_path, 'C:\offline\v0.5\Output');
    % gunzip(a)
    % fileTSV = a(1:end-3);
    % tdfread(fileTSV)
    
    % figure; plot(time_actual, time_estimate)
    figure; plot(time_actual(1:end-1), diff(time_actual))
    title(eyetracker_table_paths(iFile).name)
end
