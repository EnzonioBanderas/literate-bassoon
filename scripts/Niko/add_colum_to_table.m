script_init_study7 
a = dir(fullfile(fp_s,'**','*checklist.xlsx')) %creates a structure with list of all paticipants checklists and their locations
supertab = essbids_readTsv(fullfile(fp_scr,'ppa','test_updatePTab.tsv')) % reads a previous Tsv table
t1 = table
t1.stimulation_threshold=nan(numel(a),1)
t1.participant_id=cell(numel(a),1)
 
stimulation_threshold = nan(size(supertab,1),1)
for i=1:numel(a)
    disp(a(i).name) %dont need this, shows filename of the checklist that is processed by the loop
    [~,~,rawData] = xlsread([a(i).folder filesep a(i).name]) %reads current checklistfile (excel table)
 
    %!!!!!don't realy need next to lines

    b = rawData(:,1) % reads first column of the table from the checklistfile
    c = strcmp(rawData(:,1),'stimulation_threshold') %compares all the values of the first column
    % with string 'stimulation_threshold', gives logical array
   
    idx = find(strcmp(rawData(:,1),'stimulation_threshold'),1) %finds the row number where the value is 1 
    % which means that it is 'stimulation_threshold'
    g = rawData{idx,2} %gives the value in the the second column that corresponds to the first colum idx, which is 7 which is 'stim_treshold'
   if ~ ischar(g)
     t1.stimulation_threshold{i} = (g);% creates string array with threshold values
   end
   tmp = strsplit(a(i).name,"_")
   t1.participant_id{i} = tmp{1};
end
littletab = table(stimulation_threshold) %creates a tab of string array

megatab = [supertab, littletab] %combines to tabs
megatab=megatab(:,[1:3,end,4:end-1])


%DONT UNDERSTAND THIS, just copied P-k-code fore creating something called
%"custom property named JsonSidecar" which is neede for creating TSV table


megatab = addprop(megatab,'JsonSidecar2',{'table'}) %creates a custom property named JsonSidecar 
jinfo = megatab.Properties.CustomProperties.JsonSidecar2 %working variable 

%Define values for each variable
jinfo.participant_id.Description = 'Participant id Description'
jinfo.sfb1280_id.Description = 'Sfb1280 ID description'
jinfo.cs_plus = struct('Description','CS Plus Description','Levels',struct('diamond','shows diamond','square','shows square'))
jinfo.stim_intensity = struct('Description','Stim Description','Units','mA')
jinfo.age = struct('Description','placeholder age','Units','years')
jinfo.sex = struct('Description','sex placeholder','Levels',struct('m','male','f','female','o','other'))

%added next line
jinfo.stimulation_threshold = struct ('Description','Stim Threshold','Units','mA')
%Assign values from working variable to the JsonSidecar2 Property
megatab.Properties.CustomProperties.JsonSidecar2 = jinfo

%Change directory to store created participant table 
cd '/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/Niko/'
essbids_writeTsv('megatab.tsv',megatab) %create tsv file