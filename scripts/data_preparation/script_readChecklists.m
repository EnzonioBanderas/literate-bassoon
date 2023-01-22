fn_ptab = fullfile(fp_d,'participants.tsv');
ptab = readtable(fn_ptab,'filetype','text','delim','tab','treat','n/a');


ptab.protocol         = nan(size(ptab,1),1);
ptab.recall_order     = nan(size(ptab,1),1);
ptab.stim_tacs_group = repmat({'n/a'},size(ptab,1),1);
ptab.stim_tacs_impedance       = nan(size(ptab,1),1);
ptab.stim_us_sensory_threshold = nan(size(ptab,1),1);
ptab.stim_us_pain_threshold    = nan(size(ptab,1),1);
ptab.stim_us_pain_threshold    = nan(size(ptab,1),1);
ptab.experimenter_ses1 = repmat({'n/a'},size(ptab,1),1);
ptab.experimenter_ses2 = repmat({'n/a'},size(ptab,1),1);
ptab.comment = repmat({'n/a'},size(ptab,1),1);
for i=1:size(ptab,1)
  pID = ptab.participant_id{i};
  fn = fullfile(fp_s,pID,'documents',[pID '_checklist.xlsx']);
  if not(exist(fn,'file')==2)
    fprintf('missing : %s\n',fn);
    continue
  end
  t0 = readtable(fn);
  ind = ismember(t0.label,'participant_id');
  if not(contains(pID,t0.value(ind)))
    warning('label error :%s',fn);
    continue
  end
  
  
  ptab.protocol(i)       = str2double(t0.value{ismember(t0.label,'protocol')});
  ptab.recall_order(i)   = str2double(t0.value{ismember(t0.label,'recall_order')});
  ptab.stim_tacs_group{i} = ter_unblindNeuroconnCode(...
    t0.value{ismember(t0.label,'stim_code')});
  ptab.stim_tacs_impedance(i)       = str2double(t0.value{ismember(t0.label,'stim_impedance')});
  ptab.stim_us_sensory_threshold(i) = str2double(t0.value{ismember(t0.label,'stim_sensory_threshold')});
  ptab.stim_us_pain_threshold(i)    = str2double(t0.value{ismember(t0.label,'stim_pain_threshold')});
  
  ptab.experimenter_ses1(i) = t0.value(ismember(t0.label,'experimenter_ses1'));
  ptab.experimenter_ses2(i) = t0.value(ismember(t0.label,'experimenter_ses2'));
  ptab.comment(i) =  t0.value(ismember(t0.label,'comment'));
end

ptab.experimenter_ses1 = strrep(ptab.experimenter_ses1,' ','');
ptab.experimenter_ses2 = strrep(ptab.experimenter_ses2,' ','');