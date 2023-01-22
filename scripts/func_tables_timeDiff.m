function tableOut = func_tables_timeDiff(table1, table2)
%func_tables_timeDiff: for each point in starting and duration table 1, 
% find the closest corresponding point in starting and duration table 2
%
% Output: create a merged table with table 1 points and closest
% corresponding table 2 points (with time differences in final column)
    
nRow = size(table1, 1);
tableOut = zeros(nRow, 3);
tableOut(:, 1) = table1(:, 1);
for iM = 1:nRow
    onset_difference = table2(:, 1) - table1(iM, 1);
    [~, minAbsInd_onset_difference] = ...
        min(abs(onset_difference));
    min_onset_difference = onset_difference(minAbsInd_onset_difference);    
    tableOut(iM, 2) = table2(minAbsInd_onset_difference, 1);        
    tableOut(iM, 3) = min_onset_difference;
end

end

