% creates TSV table using the data from
%
% checklist.xlsx
% participants.tsv
% SFB-Z7T7-List.xlsx
% sub-Z7T7xxx_sfbquestionnaire.xlsx
% sub-Z7T7xxx_settings.csv
%
% made by NPE with the help of PPE code from "test_updatePTab_live" and
% assistance of TER. If needed may be replaced with other code or integrated
% to the original code of PPE
script_init_study7
a = dir(fullfile(fp_s,'**','*checklist.xlsx'))
t1 = table
t1.stimulation_threshold=nan(numel(a),1)
t1.participant_id=cell(numel(a),1)
t1.stim_intensity=nan(numel(a),1)
t1.complete_paradigm=nan(numel(a),1)
for i=1:numel(a)
   [~,~,rawData] = xlsread([a(i).folder filesep a(i).name])
   idx = find(strcmp(rawData(:,1),'stimulation_threshold'),1)
   idx2 = find(strcmp(rawData(:,1),'stimulation_final'),1)
   idx3 = find(strcmp(rawData(:,1),'complete paradigm'),1)
   g = rawData{idx,2}
   g2 = rawData{idx2,2}
   g3 = rawData{idx3,2}
      if ~ ischar(g)      
      t1.stimulation_threshold(i) = (g)
      end
      if ~ ischar(g2)
      t1.stim_intensity(i)=(g2)
      end
      if ~ ischar(g3)
      t1.complete_paradigm(i)=(g3)
      end
   tmp = strsplit(a(i).name,"_")
   t1.participant_id{i} = tmp{1}
end
Patricktab = essbids_readTsv(fullfile(fp_scr,'ppa','test_updatePTab.tsv'))
mytab = outerjoin(Patricktab,t1,"Keys","participant_id","MergeKeys",true)
mytab=mytab(:,[1:3,end-1,4,end-2,5,6,end])

mytab = addprop(mytab,'JsonSidecar2',{'table'}) %creates a custom property named JsonSidecar 
jinfo = mytab.Properties.CustomProperties.JsonSidecar2 %working variable 

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
mytab.Properties.CustomProperties.JsonSidecar2 = jinfo

%Change directory to store created participant table 
cd '/media/diskEvaluation/Evaluation/sfb1280a05study7/scripts/Niko/'
essbids_writeTsv('mytab.tsv',mytab) %create tsv file