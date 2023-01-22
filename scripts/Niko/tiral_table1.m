strs = char(randi([double('A') double('Z')], 5, 3))           % Create Variable
vTbl = table(strs, 'VariableNames',{'Strings'})                % Create ‘vTbl’
todayDate = datestr(now,'mm/dd/yyyy')                          % Create Date String
datesvct = repmat(todayDate, size(vTbl,1), 1)                  % Create ‘Vector’ Of Date String
vTbl = [ table(datesvct, 'VariableNames', {'VDate'})  vTbl]

 