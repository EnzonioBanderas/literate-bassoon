function SI_syn = func_readSIEMENSphysio(fp_physio, fp_nii)
%func_readSIEMENSphysio: SIEMENS physio load and internal synchronization
%   Load SIEMENS PMU files and synchronize with each other (the start and 
%   end times of the SIEMENS PMU files differ)

%% settings and paths
% time_sampling_PO_SIEMENS = 2.5; % sampling time in ms
ext_used = 'ext';

if isfolder(fp_physio)
    fl_SIEMENS = [dir(fullfile(fp_physio, ['*.', ext_used])); ...
        dir(fullfile(fp_physio, '*.puls')); ...
        dir(fullfile(fp_physio, '*.resp')); ...
        dir(fullfile(fp_physio, '*.ecg'))];
    % fl_SIEMENS = [dir(fullfile(fp_physio, ['*.', ext_used])); ...
    %     dir(fullfile(fp_physio, '*.puls')); ...
    %     dir(fullfile(fp_physio, '*.resp'))];
    fl_SIEMENS = func_dirl2fl(fl_SIEMENS);

    % if optional fp_nii input is given, only load relevant physio files (do not load files of other sessions)
    if exist('fp_nii', 'var')
        fp_nii_sesID = regexp(fp_nii, 'ses-[0-9]*', 'match');
        fp_nii_sesID = fp_nii_sesID{end}; % last ses code is of the filename, most relevant

        fl_SIEMENS_sesID = regexp(fl_SIEMENS, 'ses-[0-9]*', 'match');
        fl_SIEMENS_sesID = vertcat(fl_SIEMENS_sesID{:});

        fl_SIEMENS = fl_SIEMENS(strcmp(fl_SIEMENS_sesID, fp_nii_sesID));
    end
else
    fl_SIEMENS = fp_physio;
end

fprintf('SIEMENS syn structure definition\n')


%% Read in SIEMENS files
% SIEMENS physio file loop
lineData_cell = cell(length(fl_SIEMENS), 1);
logFooter_cell = cell(length(fl_SIEMENS), 1);
linesFooter_cell = cell(length(fl_SIEMENS), 1);
tableData = cell(length(fl_SIEMENS), 1);
startTime = zeros(length(fl_SIEMENS), 1);
startTime_estimate = zeros(length(fl_SIEMENS), 1);
stopTime = zeros(length(fl_SIEMENS), 1);
stopTime_estimate = zeros(length(fl_SIEMENS), 1);
physio_SIEMENS_data = cell(length(fl_SIEMENS), 1);
for i = 1:length(fl_SIEMENS)
    tic

    [~, fn_SIEMENS, ext_SIEMENS] = fileparts(fl_SIEMENS{i});
    
    % get index of file with scantrigger information
    if strcmp(ext_SIEMENS, '.ext')
        ext_index = i;
    end
    
    fprintf(['reading ', fn_SIEMENS, '\n'])
    
    %readPMU, MRIsyn cutout PMU +-20s
    
    
%     % reading method 1
%     [lineData, logFooter, linesFooter] = tapas_physio_read_physlogfiles_siemens_raw(fl_SIEMENS{i});
%     lineData_cell(i) = {lineData};
%     logFooter_cell(i) = {logFooter};
%     linesFooter_cell(i) = {linesFooter};
    
    % reading method 2, compare to 1
    fid             = fopen(fl_SIEMENS{i});
    C               = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    linesFooter = C{1}(2:end);
    logFooter.LogStartTimeSeconds =   str2num(char(regexprep(linesFooter(~cellfun(@isempty,strfind(linesFooter,...
        'LogStartMDHTime'))),'\D',''))) / 1000;
    logFooter.LogStopTimeSeconds =    str2num(char(regexprep(linesFooter(~cellfun(@isempty,strfind(linesFooter,...
        'LogStopMDHTime'))),'\D',''))) / 1000;
    logFooter.ScanStartTimeSeconds = str2num(char(regexprep(linesFooter(~cellfun(@isempty,strfind(linesFooter,...
        'LogStartMPCUTime'))),'\D','')));
    logFooter.ScanStopTimeSeconds = str2num(char(regexprep(linesFooter(~cellfun(@isempty,strfind(linesFooter,...
        'LogStopMPCUTime'))),'\D','')));
    lineData = C{1}{1};
    logFooter_cell(i) = {logFooter};
    linesFooter_cell(i) = {linesFooter};
    
%     %%%%%%%%%%%%%
    % signals start of data logging
    iTrigger = regexpi(lineData, ' 6002 ');

    if ~isempty(iTrigger)
        % crop string after trigger
        lineData = lineData((iTrigger(end)+6):end);
        doCropLater = false;
    else
        % crop first 4 values as in UPenn manual after conversion
        doCropLater = true;
    end
    lineData_cell(i) = {lineData};
% 
%     data = textscan(lineData, '%d', 'Delimiter', ' ', 'MultipleDelimsAsOne',1);
% 
%     if doCropLater
%         % crop first 4 values;
%         data{1} = data{1}(5:end);
%     end
%     %%%%%%%%%%%%%%
    
    if contains(lineData, ': ')
        lineData_split = strsplit(lineData, ': ');
        lineData_split = cellfun(@(x) strsplit(x, ' '), lineData_split, 'uni', 0);
        lineData_split{1} = lineData_split{1}(1:end-1);
        for j=2:(length(lineData_split)-1)
            lineData_split{j} = lineData_split{j}(2:end-1);
        end
        lineData_split{end} = lineData_split{end}(2:end);
        lineData_split = [lineData_split{:}];
    else
        lineData_split = strsplit(lineData, ' ');
    end
%     lineData_uniq = unique(lineData_split);
%     uniqData2{i} = table(lineData_uniq', ...
%         cell2mat(cellfun(@(x) str2double(x), lineData_uniq, 'uni', 0))', ...
%         cell2mat(cellfun(@(x) sum(strcmp(lineData_split, x)), lineData_uniq, 'uni', 0))', ...
%         'VariableNames', {'uniq', 'uniq_double', 'count'});
    
    % remove first 4 integers which contain acquisition information
%     doCropLater = true;
    if doCropLater
        lineData_data = lineData_split(5:end); 
    else
        lineData_data = lineData_split;
    end

    % convert strings to double
    lineData_data = cellfun(@(x) str2double(x), lineData_data, 'uni', 0);
    lineData_data = [lineData_data{:}]';
    lineData_data = lineData_data(~isnan(lineData_data));
    
    % get number of elements before removing triggers
    nSamples_withTriggers = length(lineData_data);
    
    % get markers/triggers evaluated by system
    nDataMarker = length(lineData_data);
    lineData_markers_logical = ismember(lineData_data, [5000, 5002, 5003, 6000, 6002]);
    iMarkers = find(lineData_markers_logical); % indices of markers
    cMarkers = lineData_data(iMarkers); % values of markers
    
%     cellfun(@(x) str2double(x), cell, 'uni', 0)
    
    [~, ~, extension] = fileparts(fl_SIEMENS{i});
    if strcmp(extension, '.ecg')
        
        lineData_markers_ind = find(lineData_markers_logical)-1;
        
        % There can be a marker on the first position, which when 1 is subtracted would be a 0
        % For this exception we can say that the lineData_marker_logical is
        % false, although we could also find another solution for it
        lineData_markers_logical(lineData_markers_ind(lineData_markers_ind <= 0)+1) = false;
        
        nMark = sum(lineData_markers_logical);
        nData = sum(~lineData_markers_logical);
        
        lineData_newMarker = zeros(nDataMarker, 1);
        lineData_newMarker(iMarkers-1) = cMarkers;
        
        lineDataMarker_colIndex = NaN(nDataMarker, 1);
        lineDataMarker_colIndex(~lineData_markers_logical) = 1:nData;
        lineDataMarker_colIndex = mod(lineDataMarker_colIndex, 4);
        lineDataMarker_colIndex(lineDataMarker_colIndex == 0) = 4;
        
        lineDataMarker_index = 1:length(lineData_markers_logical);
        
        % When there is a marker, go to the index before and copy the marker value to a copy of the original linedata with  
        
        % Loop through markers
        zeros(nMark, 1)
        for iMark = 1:nMark
            Marker_Index(iMark) = ...
                lineDataMarker_index(iMarkers(iMark) - 1); %%%%%%%%%%%%%%%%
%             Marker_colIndex(iMark) = ...
%                 lineDataMarker_colIndex(iMarkers(iMark) - 1);
        end
           
        % remove markers from lineData_data
        lineData_data = lineData_data(~lineData_markers_logical);
        lineData_newMarker = lineData_newMarker(~lineData_markers_logical);
        
        iDataStream = 1:nSamples_withTriggers;
        iDataStream(lineData_markers_logical) = [];
        
        nSamples = length(lineData_data);
        nRows = ceil(nSamples/4);
        
        tableData{i} = zeros(nRows,5);
        lineData_newMarker_4CHsplit = zeros(nRows,4);
%         iData_table = zeros(nRows,4);
        
        % split into 4 channels instead
        tableData{i}(1:length(lineData_data(1:4:end)), 1) =  ...
            lineData_data(1:4:end);
        tableData{i}(1:length(lineData_data(2:4:end)), 2) =  ...
            lineData_data(2:4:end);
        tableData{i}(1:length(lineData_data(3:4:end)), 3) =  ...
            lineData_data(3:4:end);
        tableData{i}(1:length(lineData_data(4:4:end)), 4) =  ...
            lineData_data(4:4:end);
        lineData_newMarker_4CHsplit(1:length(lineData_newMarker(1:4:end)), 1) = ...
            lineData_newMarker(1:4:end);
        lineData_newMarker_4CHsplit(1:length(lineData_newMarker(2:4:end)), 2) = ...
            lineData_newMarker(2:4:end);
        lineData_newMarker_4CHsplit(1:length(lineData_newMarker(3:4:end)), 3) = ...
            lineData_newMarker(3:4:end);
        lineData_newMarker_4CHsplit(1:length(lineData_newMarker(4:4:end)), 4) = ...
            lineData_newMarker(4:4:end);
        
        % generalize the markers from 4 channels to 1 (is this valid?)
        tableData{i}(:, 5) = max(lineData_newMarker_4CHsplit, [], 2);
        
        % are the markers really separate across channels? Can they be
        % generalized to a single column, or should they remain as separate
        % columns? Maybe if they can not be generalized, just use max!
        
%         iData_table(1:length(iDataStream(1:4:end)), 1) =  iDataStream(1:4:end);
%         iData_table(1:length(iDataStream(2:4:end)), 2) =  iDataStream(2:4:end);
%         iData_table(1:length(iDataStream(3:4:end)), 3) =  iDataStream(3:4:end);
%         iData_table(1:length(iDataStream(4:4:end), 4) =  iDataStream(4:4:end);
        
%         if mod(nSamples,2) == 1
%             tableData{i}(1:nRows-1,2) = lineData_data(2:2:end);
%             iData_table(1:nRows-1,2) = iDataStream(2:2:end);
%         else
%             tableData{i}(1:nRows,2) = lineData_data(2:2:end);
%             iData_table(1:nRows,2) = iDataStream(2:2:end);
%         end
        
        % now fill up 3rd column with trigger data
        % - for each trigger index in data{1}, check where ECG data with closest
        % smaller index ended up in the data_table ... and put trigger code in
        % same row of that table
%         nTriggers = numel(iMarkers);
        
        % create index without markers to get column position with mod(x,
        % 4) (if 0 then it should be 4)
        
        % loop through each marker
        % the index before
        
        % position - 1 in index without markers is the row index of the
        % undivided (still one line) marker vector
        
%         for iTrigger = 1:nTriggers
%             % find index before trigger event in data stream and
%             % detect it in table
%             % look in 1st column as well whether maybe signal detected there
%             iRow = [];
%             iCol = 1;
%             while isempty(iRow)
%                 iRow = find(iData_table(:,iCol) == iMarkers(iTrigger)-1);
%                 iCol = iCol+1;
%             end
%             tableData{i}(iRow,5) = cMarkers(iTrigger);
%         end
        
    else
    
        % get table data from lineData (currently not ecg compatible!)
        lineData_markerdata = zeros(length(lineData_data), 1);
        lineData_seconddata = zeros(length(lineData_data), 1); %%%%%%%%%%%%%%%%%%%%%%%%
        lineData_markers_ind = find(lineData_markers_logical)-1;
        
        % There can be a marker on the first position, which when 1 is subtracted would be a 0
        % For this exception we can say that the lineData_marker_logical is
        % false, although we could also find another solution for it
        lineData_markers_logical(lineData_markers_ind(lineData_markers_ind <= 0)+1) = false;
        
        lineData_markerdata(lineData_markers_ind(lineData_markers_ind > 0)) = lineData_data(lineData_markers_logical);
        lineData_data = lineData_data(~lineData_markers_logical);
        lineData_markerdata = lineData_markerdata(~lineData_markers_logical);
        lineData_seconddata = lineData_seconddata(~lineData_markers_logical);
        tableData{i} = [lineData_data, ...
            lineData_seconddata, lineData_markerdata];
        
    end
    
    % define sampling time
    time_sampling_PO_SIEMENS = 2.5;
%     % define sampling time depending on file extension
%     switch extension
%         case '.ext'
%             time_sampling_PO_SIEMENS(i) = 2.5;
%         case '.puls'
%             time_sampling_PO_SIEMENS(i) = 2.5;
%         case '.resp'
%             time_sampling_PO_SIEMENS(i) = 2.5; %?
%         case '.ecg'
%             time_sampling_PO_SIEMENS(i) = 2.5;
%     end
    
    % estimate time vector using defined sampling time
    Time = (logFooter_cell{i}.ScanStopTimeSeconds-(size(tableData{i}, 1) - 1)* ...
        time_sampling_PO_SIEMENS: ...
        time_sampling_PO_SIEMENS: ...
        logFooter_cell{i}.ScanStopTimeSeconds)';
    
    % get start and end time from the header and the estimate
    startTime(i) = logFooter_cell{i}.ScanStartTimeSeconds;
    startTime_estimate(i) = Time(1);
    stopTime(i) = logFooter_cell{i}.ScanStopTimeSeconds;
    stopTime_estimate(i) = Time(end);
    
    % add time vector to table data
    tableData{i} = [tableData{i}, Time];

    % elapsed time for reading in file
    toc
end

% figure; plot(tableData{1}(:,end), tableData{1}(:,1)); hold on; plot(tableData{1}(:,end), tableData{1}(:,3)/10000);
% figure; plot(tableData{2}(:,end), tableData{2}(:,1)); hold on; plot(tableData{2}(:,end), tableData{2}(:,3));
% figure; plot(tableData{3}(:,end), tableData{3}(:,1)); hold on; plot(tableData{3}(:,end), tableData{3}(:,2)); plot(tableData{3}(:,end), tableData{3}(:,3)); plot(tableData{3}(:,end), tableData{3}(:,4)); plot(tableData{3}(:,end), tableData{3}(:,5));

%% calculate start and end time estimates for within SIEMENS PMU sync
% startTime_estimate2 = startTime_estimate;
% stopTime_estimate2 = stopTime_estimate;
startTime_estimate = max(startTime_estimate);
stopTime_estimate = min(stopTime_estimate);

% check and synchronize
for i = 1:length(fl_SIEMENS)

    [~, fn_SIEMENS] = fileparts(fl_SIEMENS{i});
    fprintf(['checking ', fn_SIEMENS, '\n'])

    % do some time checks
    clockRatio = (abs(logFooter_cell{i}.LogStartTimeSeconds - logFooter_cell{i}.LogStopTimeSeconds)*1) / ...
        (abs(logFooter_cell{i}.ScanStartTimeSeconds - logFooter_cell{i}.ScanStopTimeSeconds)/1e3);
    time_sampling_pulse_MPCU = (abs(logFooter_cell{i}.ScanStartTimeSeconds - logFooter_cell{i}.ScanStopTimeSeconds)/(size(tableData{i}, 1)-1))/1e3;
    time_sampling_pulse_MDH = abs(logFooter_cell{i}.LogStartTimeSeconds - logFooter_cell{i}.LogStopTimeSeconds)/(size(tableData{i}, 1)-1);
    time_sampling_pulse_compensated = (time_sampling_PO_SIEMENS / 1000) * clockRatio;
    fprintf('time sampling estimated MDH = %4.16f\n', time_sampling_pulse_MDH * 1000)
    fprintf('time sampling estimated MPCU = %4.16f\n', time_sampling_pulse_MPCU * 1000)
    fprintf('time sampling compensated = %4.16f\n', time_sampling_pulse_compensated * 1000)
    fprintf('time sampling defined = %4.16f\n', time_sampling_PO_SIEMENS)
    points_expected = (logFooter_cell{i}.ScanStopTimeSeconds - logFooter_cell{i}.ScanStartTimeSeconds) / time_sampling_PO_SIEMENS;
    points_measured = size(tableData{i}, 1);
    fprintf('number of missing points = %4.16f\n', points_expected-points_measured)
    fprintf(['SIEMENS physio start time is ', ...
        datestr(func_MPCU2datetime(logFooter_cell{i}.ScanStartTimeSeconds)), ...
        '\n'])
    fprintf(['SIEMENS physio stop time is ', ...
        datestr(func_MPCU2datetime(logFooter_cell{i}.ScanStopTimeSeconds)), ...
        '\n'])
    fprintf('difference between end times in seconds = %4.16f\n', ...
        (tableData{i}(end, end)-logFooter_cell{i}.ScanStopTimeSeconds) / 1000)
    fprintf('\n')
    fprintf('difference between start times in seconds = %4.16f\n', ...
        (tableData{i}(1, end)-logFooter_cell{i}.ScanStartTimeSeconds) / 1000)
    fprintf('\n')
    
    % only keep data where all signals are measured
    tableData{i} = tableData{i}(tableData{i}(:, end) >= startTime_estimate & tableData{i}(:, end) <= stopTime_estimate, :);

    [~, ~, extension] = fileparts(fl_SIEMENS{i});
    switch extension
        case '.ext'
            physio_SIEMENS_data{i} = table(tableData{i}(:, 1), ...
                tableData{i}(:, 3), ...
                tableData{i}(:, 4), ...
                'VariableNames', {'ext_signal', 'ext_code', 'time'});
        case '.puls'
            physio_SIEMENS_data{i} = table(tableData{i}(:, 1), ...
                tableData{i}(:, 3), ...
                tableData{i}(:, 4), ...
                'VariableNames', {'puls_signal', 'puls_code', 'time'});
        case '.resp'
            physio_SIEMENS_data{i} = table(tableData{i}(:, 1), ...
                tableData{i}(:, 3), ...
                tableData{i}(:, 4), ...
                'VariableNames', {'resp_signal', 'resp_code', 'time'});
        case '.ecg'
            physio_SIEMENS_data{i} = table(tableData{i}(:, 1), ...
                tableData{i}(:, 2), ...
                tableData{i}(:, 3), ...
                tableData{i}(:, 4), ...
                tableData{i}(:, 5), ...
                tableData{i}(:, 6), ...
                'VariableNames', {'ecg_signal1', 'ecg_signal2', 'ecg_signal3', 'ecg_signal4', 'ecg_code', 'time'});
    end
        

end

% plot
% figure;
% plot(tableData{1}(:,end), tableData{1}(:,3)), hold on
% plot(tableData{2}(:,end), tableData{2}(:,3))
% plot(tableData{3}(:,end), tableData{3}(:,3))
% title('SIEMENS')
% 
% figure;
% plot(tableData{1}(:,end), tableData{1}(:,1)), hold on
% plot(tableData{2}(:,end), tableData{2}(:,1))
% plot(tableData{3}(:,end), tableData{3}(:,1))
% title('SIEMENS')
% 
% figure;
% plot(physio_SIEMENS_data.time, physio_SIEMENS_data.pulse_signal), hold on
% plot(physio_SIEMENS_data.time, physio_SIEMENS_data.pulse_code/5000)
% title('SIEMENS')

% % calculate pulse count
% extCounts = arrayfun(@(x) length(find(physio_SIEMENS_data{1}.ext_code == x)), unique(physio_SIEMENS_data{1}.ext_code), 'uni', 0);
% extCount = extCounts{ismember(unique(physio_SIEMENS_data{1}.ext_code), 5000)};

% split data into ext_table and data_table
scantrigger_table = physio_SIEMENS_data{ext_index}(physio_SIEMENS_data{ext_index}.ext_code==5000, 'time');
scantrigger_table.index = (1:size(scantrigger_table, 1))';
% scantrigger_table.duration = NaN(size(scantrigger_table, 1), 1);
scantrigger_table.duration = ones(size(scantrigger_table, 1), 1) * 30;
% add scantrigger durations?

% shorten data_table so that all modalities have equal amount of points 
physio_SIEMENS_data_length = cellfun(@(x) size(x, 1), physio_SIEMENS_data, 'uni', 0);
physio_SIEMENS_data_length = [physio_SIEMENS_data_length{:}];
physio_SIEMENS_data_length = min(physio_SIEMENS_data_length);
physio_SIEMENS_data = cellfun(@(x) x(1+end-physio_SIEMENS_data_length:end, :), physio_SIEMENS_data, 'uni', 0);

% add to syn structure
SI_syn.fp = fl_SIEMENS;
SI_syn.scantrigger_table = scantrigger_table;
SI_syn.data_table = physio_SIEMENS_data;
SI_syn.data_modality = 'SIEMENS';

% display some basic info
fprintf(['SIEMENS physio start time is ', ...
    datestr(func_MPCU2datetime(physio_SIEMENS_data{1}.time(1))), ...
    '\n'])
fprintf(['SIEMENS physio stop time is ', ...
    datestr(func_MPCU2datetime(physio_SIEMENS_data{1}.time(end))), ...
    '\n'])
fprintf('\n')

% % optional plot
% syn_data_table = [syn.data_table{1}(:, 1:end-1), syn.data_table{2}(:, 1:end-1), syn.data_table{3}(:, 1:end-1)];
% figure, plot(table2array(syn.data_table{1}(:, end)), table2array(syn_data_table(:, [1,2])))
% figure, plot(table2array(syn.data_table{1}(:, end)), table2array(syn_data_table(:, [3,4])))
% figure, plot(table2array(syn.data_table{1}(:, end)), table2array(syn_data_table(:, [5,6,7,8,9])))

end

