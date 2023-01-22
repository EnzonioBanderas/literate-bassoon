%%TO-DO

%% Loading relevant paths
script_init_study7 %load in all packages and paths
fp_m = fullfile(fp0,"misc") %load misc folder

%% Get Excel-File from misc folder and only select Z7T7 and SFB-Code
etab = readtable(fullfile(fp_m,"SFB-Z7T7-List.xlsx"),"VariableNamingRule","preserve") %get excel sfb-z7t7 list
clcetab = etab(:,["Z7T7-Code" "SFB-Code"]) %only select participant codes from the excel spreadshet
etab = renamevars(etab,["Z7T7-Code" "SFB-Code"],["participant_id" "sfb1280_id"]) %change variable names

%% Get filelists 
fl_set = dir(fullfile(fp_s,'**','*_settings.csv')) %get file list of all settings (applies to newer participants)
fl_setOld = struct2table(dir(fullfile(fp_s,'**','*_stimsettingratings.txt'))) %get stimsettinglogs
fl_que = dir(fullfile(fp_s,'**','*_sfbquestionnaire.xlsx')) %get file list of all sfbquestionnaires

%% Loop through each stimsettings to aquire CSplus, stimulus intensity and participant ID
% Preallocate space for each variable
array_ParticipantID = strings(size(fl_set,1),1) %creates an empty string array with one column and one row for each file
array_CSplus = strings(size(fl_set,1),1) %see above
array_StimIntensity = nan(size(fl_set,1),1) %creates an array filled with zeros size: see above

% Looping through the file list and extracting values
for idx = 1:size(fl_set,1)  
    fcur = fileread([fl_set(idx).folder filesep fl_set(idx).name]) %open file
    %get current participant id and store in array
    f_participant_id = regexp (fcur,'ID\s*: (?<val>[a-zA-Z0-9\-]+)','names')
    array_ParticipantID(idx,1) = string(f_participant_id.val)
    %get current csplus and store in array
    f_cs_plus = regexp (fcur,'CSplus\s*: (?<val>[a-z]+)','names') 
    array_CSplus(idx,1) = string(f_cs_plus.val)
    %get current us intensity and store in array
    f_stim_intensity = regexp (fcur,'US_mA\s*: (?<val>[0-9\.]+)','names')
    array_StimIntensity(idx,1) = double(string(f_stim_intensity.val))
end

% Create table and 
stab = table(array_ParticipantID,array_CSplus,array_StimIntensity)
stab = renamevars(stab,["array_ParticipantID" "array_CSplus" "array_StimIntensity"],["participant_id" "cs_plus" "stim_intensity"])

%% Loop through each questionnaire excel file to participant id, age and sex
% Preallocate space for each variable
array_ParticipantID = strings(size(fl_que,1),1) %creates an empty string array with one column and one row for each file
array_Age = nan(size(fl_que,1),1) %creates an array filled with zeros size: see above
array_Sex = strings(size(fl_que,1),1) %see array_ParticipantID comment

% Looping through the file list and extracting values
for idx = 1:size(fl_que,1)
    importOpt = detectImportOptions([fl_que(idx).folder filesep fl_que(idx).name]) %get import options for current file
    importOpt.VariableTypes{2} = 'char' %change import option for the second column to character array
    importOpt.VariableNamingRule = "preserve" %change variable naming rule to perserve (avoids warning msg)
    % open and flip the table that is read from the current file
    fcur = rows2vars(readtable([fl_que(idx).folder filesep fl_que(idx).name],importOpt),"VariableNamesSource",'label')
    %get name of the current file and extract the participant id
    fname = fl_que(idx).name
    array_ParticipantID(idx,1) = extractBefore(string(fname),"_sfbquestionnaire.xlsx")
    % get current age and sex and store in respective array
    array_Age(idx,1) = double(string(fcur.age{1,1}))
    array_Sex(idx,1) = string(fcur.sex{1,1})
end

% Create questionnaire table and rename variables
qtab = table(array_ParticipantID,array_Age,array_Sex)
qtab = renamevars(qtab,["array_ParticipantID" "array_Age" "array_Sex"],["participant_id" "age" "sex"])

% %% Old Stimsettings
% %Preallocate space for loop variables
% array_log = false(size(fl_setOld,1),1)
% array_ParticipantID = strings(size(fl_setOld,1),1)
% for idx1 = 1:size(fl_setOld,1) %Looping through filelist
%     curr_participant_id = string(extractBefore(fl_setOld.name(idx1),'_ses')) %get participant id from file name
%     if ismember(curr_participant_id,array_ParticipantID) || ismember(curr_participant_id,stab.participant_id)
%         array_ParticipantID(idx1) = curr_participant_id
%         array_log(idx1,1) = false
%         continue
%     end
%     array_ParticipantID(idx1) = curr_participant_id
%     for idx2 = 1:size(etab,1) %Looping through excel table
%         e_participant_id = string(etab.participant_id(idx2)) %get participants from excel table
%         if isequal(e_participant_id,curr_participant_id)
%             array_log(idx1,1) = true
%         end
%     end
% end
% fl_setOld = fl_setOld(array_log,:)
% %% Aquire data TODO
% for idx1 = 1:size(fl_setOld,1)
%     fcur = fileread([fl_setOld.folder{idx1} filesep fl_setOld.name{idx1}]) %open file
%     f_participant_id = regexp (fcur,'ID\s*: (?<val>[a-zA-Z0-9\-]+)','names')
% 
%     ftab = readtable([fl_setOld.folder{idx1} filesep fl_setOld.name{idx1}])
%     if contains(ftab.Properties.VariableNames,{'Rating'})
%         a = 0
%         for idx2 = 1:size(ftab,1)
%             if string(ftab.Rating{idx2}) == "not sensed"
%                 a = a +1
%                 continue
%             end
%             a = a+1
%             break
%         end
%         f_stim_intensity = regexp (fcur,'final stimulation intensity\s*: (?<val>[0-9\.]+)','names')
%         treshhold_value = double(extractBefore(string(ftab.StimIntensity{a})," mA"))
%     else
%         warning("Stimsettings-File empty for: " + f_participant_id.val)
%     end
% end
% 

%% Join tables using participant id as key
ptab = outerjoin(etab,stab,"Keys","participant_id","MergeKeys",true)
ptab = outerjoin(ptab,qtab,"Keys","participant_id","MergeKeys",true)

%% Create JsonSidecar
ptab = addprop(ptab,'JsonSidecar',{'table'}) %creates a custom property named JsonSidecar 
jinfo = ptab.Properties.CustomProperties.JsonSidecar %working variable 

%Define values for each variable
jinfo.participant_id.Description = 'Participant id Description'
jinfo.sfb1280_id.Description = 'Sfb1280 ID description'
jinfo.cs_plus = struct('Description','CS Plus Description','Levels',struct('diamond','shows diamond','square','shows square'))
jinfo.stim_intensity = struct('Description','Stim Description','Units','mA')
jinfo.age = struct('Description','placeholder age','Units','years')
jinfo.sex = struct('Description','sex placeholder','Levels',struct('m','male','f','female','o','other'))

%Assign values from working variable to the JsonSidecar Property
ptab.Properties.CustomProperties.JsonSidecar = jinfo

%Change directory to store created participant table in ppa folder
cd '/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/ppa/'
essbids_writeTsv('test_updatePTab.tsv',ptab) %create tsv file