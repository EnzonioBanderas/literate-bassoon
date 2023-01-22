
script_init_study7;
% fp_s    = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/sourcedata2';
fp_bids = fp_d;
% fp_bids = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata';
% fp_de   = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/derivatives';
fl_acq = essbids_listFiles(fullfile(fp_s,'**','sub-*.acq'));

% remove troublesome individuals
fl_acq = fl_acq(not(contains(fl_acq,{'Z7T7095','Z7T7096','Z7T7104',...
  'Z7T7116','Z7T7146','Z7T7146','Z7T7158'}))); %,'Z7T7161','Z7T7162'
% fl_acq = fl_acq(not(contains(fl_acq,{'Z7T7095','Z7T7096','Z7T7104',...
%   'Z7T7116','Z7T7146','Z7T7146','Z7T7158','Z7T7161','Z7T7162'}))); %


%% import acq files
ws = warning('query','essbids:MultiplesDataType');
warning('off','essbids:MultiplesDataType');
error_acq = false(size(fl_acq));
for i=1:numel(fl_acq)
  % check wether output exists and if it does skip
  fn_acq = fl_acq{i};
  bf = essbids_parseLabel(fn_acq);
  bf.type = 'func';
  fp_out = fullfile(fp_bids,...
    fileparts(essbids_buildFileName(bf,'relativepath',true)));
  fl = dir(fullfile(fp_out,sprintf('sub-*_ses-%s*_physio*',bf.ses)));
  if numel(fl)>0
    %fprintf('Output present, skipping %s\n',fn_acq)
    continue;
  end
  % import or complain if it is not working
  try
    essbids_biopac2bids(fn_acq,fp_bids,...
      'TR',1.62,...
      'bidsflags',struct('task','fear','acq','stxtr1620'),...
      'minTrialPerRun',60);
  catch
    fprintf('Trouble importing %s\n',fn_acq)
    error_acq(i)=true;
  end
end
warning(ws.state,'essbids:MultiplesDataType');
fprintf('trouble with these files:\n');
disp(fl_acq(error_acq));



%% corrections
fp_7162 = fullfile(fp_bids,'sub-Z7T7162','ses-2','func');
if exist(fullfile(fp_7162,'sub-Z7T7162_ses-2_task-fear_acq-stxtr1620_run-3_date-20220520_physioevents.tsv'),'file')==2
  try
  delete(fullfile(fp_7162,'sub-*run-2*physio*'));
  movefile(...
    fullfile(fp_7162,'sub-Z7T7162_ses-3_task-fear_acq-stxtr1620_run-3_date-20220520_physioevents.tsv'),...
    fullfile(fp_7162,'sub-Z7T7162_ses-3_task-fear_acq-stxtr1620_run-2_physioevents.tsv'));
  movefile(...
    fullfile(fp_7162,'sub-Z7T7162_ses-3_task-fear_acq-stxtr1620_run-3_date-20220520_physio.tsv.gz'),...
    fullfile(fp_7162,'sub-Z7T7162_ses-3_task-fear_acq-stxtr1620_run-2_physio.tsv.gz'));
  movefile(...
    fullfile(fp_7162,'sub-Z7T7162_ses-3_task-fear_acq-stxtr1620_run-3_date-20220520_physio.json'),...
    fullfile(fp_7162,'sub-Z7T7162_ses-3_task-fear_acq-stxtr1620_run-2_physio.json'));
  catch
  end
end

fp_7172 = fullfile(fp_bids,'sub-Z7T7172','ses-3','func');
if exist(fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-3_physioevents.tsv'),'file')==2
  try
  delete(fullfile(fp_7172,'sub-*run-2*physio*'));
  movefile(...
    fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-3_physioevents.tsv'),...
    fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-2_physioevents.tsv'));
  movefile(...
    fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-3_physio.tsv.gz'),...
    fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-2_physio.tsv.gz'));
  movefile(...
    fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-3_physio.json'),...
    fullfile(fp_7172,'sub-Z7T7172_ses-3_task-fear_acq-stxtr1620_run-2_physio.json'));
  catch
  end
end
fp_7173 = fullfile(fp_bids,'sub-Z7T7173','ses-3','func');
if exist(fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-3_physioevents.tsv'),'file')==2
  try
  delete(fullfile(fp_7173,'sub-*run-2*physio*'));
  movefile(...
    fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-3_physioevents.tsv'),...
    fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-2_physioevents.tsv'));
  movefile(...
    fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-3_physio.tsv.gz'),...
    fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-2_physio.tsv.gz'));
  movefile(...
    fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-3_physio.json'),...
    fullfile(fp_7173,'sub-Z7T7173_ses-3_task-fear_acq-stxtr1620_run-2_physio.json'));
  catch
  end
end


%% check plausible run indices
fl_pe = essbids_listFiles(fullfile(fp_bids,'**',...
  'sub-*acq-stxtr1620*_physioevents.tsv'));
bf = essbids_parseLabel(fl_pe);
ind_err = false(size(bf));
for i=1:numel(fl_pe)
  try
    switch bf{i}.ses
      case '1'
        ind_err(i) = not(ismember(bf{i}.run,{'1','2'}));
      case '2'
        ind_err(i) = not(ismember(bf{i}.run,{'1'}));
      case '3'
        ind_err(i) = not(ismember(bf{i}.run,{'1','2'}));
    end
  catch
    ind_err(i) = true;
  end
end
if sum(ind_err)>0
  fprintf('run index is off here, please check manually:\n')
  disp(fl_pe(ind_err))
  warning('run index off')
end


%% remove date field from file names
fl_wDateField = essbids_listFiles(fullfile(fp_bids,'**','sub-*date-*'));
fl_new = regexprep(fl_wDateField,'_date-[0-9]*','');
for i=1:numel(fl_wDateField)
  if exist(fl_new{i},'file')==2
    fprintf('deleting %s\n',fl_wDateField{i});
    delete(fl_wDateField{i});
  else
    fprintf('renaming %s\n',fl_wDateField{i});
    movefile(fl_wDateField{i},fl_new{i});
  end
end
  

%% merge physioevents with stimlog
fl_pe = essbids_listFiles(fullfile(fp_bids,'**',...
  'sub-*acq-stxtr1620*_physioevents.tsv'));

fl_pe = fl_pe(~contains(fl_pe,{'sub-Z7T7161','sub-Z7T7172','sub-Z7T7173'}));
bf = essbids_parseLabel(fl_pe);
msuccess = false(size(fl_pe));
for i=1:numel(fl_pe)
  fl_sl = essbids_listFiles(fullfile(fp_s,['sub-' bf{i}.sub],'preslog',...
    'sub-*stimlog.txt'));
  fl_sl=fl_sl(contains(fl_sl,['ses-' bf{i}.ses]));
  fl_sl=fl_sl(contains(fl_sl,['run-' bf{i}.run]));
  if numel(fl_sl)==1
    try 
      msuccess(i) = func_mergePhysioEventsAndLog(fl_pe{i},fl_sl{1});
      if msuccess(i)
        delete(fl_pe{i});
      end
    catch
      fprintf('ERROR importing %s\n',fl_pe{i});
    end
  else
    warning('Check multiple preslog files:')
    disp(fl_sl);
  end
end
fprintf('\nTroublesome physioevent files:\n');
disp(fl_pe(not(msuccess)));


%% now prepare for EDA evaluation
func_prepareEda4rub(fp_bids,fp_de);


