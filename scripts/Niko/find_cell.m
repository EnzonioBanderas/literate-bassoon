script_init_study7 
a = dir(fullfile(fp_s,'**','*checklist.xlsx'))
supertab = essbids_readTsv(fullfile(fp_scr,'ppa','test_updatePTab.tsv'))


for i=1:numel(a)
    disp(a(i).name);
    [~,~,rawData] = xlsread([a(i).folder filesep a(i).name])
    b = rawData(:,1)
    c = strcmp(rawData(:,1),'stimulation_threshold')
    idx = find(strcmp(rawData(:,1),'stimulation_threshold'),1)
    g = rawData(idx,2)
    stimulation_threshold(i,1) = string(g)

end
littletab = table(stimulation_threshold)

megatab = [supertab, littletab]

megatab = addprop(megatab,'JsonSidecar2',{'table'}) %creates a custom property named JsonSidecar 
jinfo = megatab.Properties.CustomProperties.JsonSidecar2 %working variable 

%Define values for each variable
jinfo.participant_id.Description = 'Participant id Description'
jinfo.sfb1280_id.Description = 'Sfb1280 ID description'
jinfo.cs_plus = struct('Description','CS Plus Description','Levels',struct('diamond','shows diamond','square','shows square'))
jinfo.stim_intensity = struct('Description','Stim Description','Units','mA')
jinfo.age = struct('Description','placeholder age','Units','years')
jinfo.sex = struct('Description','sex placeholder','Levels',struct('m','male','f','female','o','other'))
jinfo.stimulation_threshold = struct ('Description','Stim Threshold','Units','mA')
%Assign values from working variable to the JsonSidecar2 Property
megatab.Properties.CustomProperties.JsonSidecar2 = jinfo

%Change directory to store created participant table in ppa folder
cd '/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/Niko/'
essbids_writeTsv('megatab.tsv',megatab) %create tsv file