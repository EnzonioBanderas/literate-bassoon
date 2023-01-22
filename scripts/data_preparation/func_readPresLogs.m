function func_readPresLogs(fp_s,fp_d)
% reads presentation log files 
% inputs : fp_s  -  path to sourcedata folder
%          fp_d  -  path to rawdata folder

fprintf('\n%s\n%s - reading presentation log files\n',...
  repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))


labelSessID = {'ses-1','ses-2','ses-3'};
nSessFiles  = [2 1 2];  % expected number of session log files 
nQuestFiles = [2 1 2];  % expected number of questionnaire log files 
subdirLogs  = 'stimlogs';

id_runlogs   = '*_stimlog.txt';
id_questlogs = '*_questlog.txt';
id_questresp = '*_questanswers.txt';

if nargin < 2
  dir_scripts = fileparts(fileparts(mfilename('fullpath')));
  fp0  = fileparts(dir_scripts);
  fp_d = fullfile(fp0,'rawdata');
  fp_s = fullfile(fp0,'sourcedata');
end

fn_ptab = fullfile(fp_d,'participants.tsv');
ptab = ter_readptab(fn_ptab);
%ptab = readtable(fn_ptab,'FileType','text','delim','tab');

for i=1:size(ptab,1)
  pID = ptab.participant_id{i};
  sID = struct2cell(dir(fullfile(fp_s,pID,'ses-*')));
  sID = sID(1,:);
  if isempty(sID)
    sID = {''};
  end
  for j=1:numel(sID)
    dirW = fullfile(fp_s,pID,sID{j},subdirLogs);
    fn_mat = strrep(fullfile(dirW,[pID '_' sID{j} '_logdataRuns.mat']),...
      '__','_');
    if exist(fn_mat,'file')==2
      continue % skip existing sets
    end
    % find logs to read
    try
      fl = ter_listFiles(dirW,id_runlogs,0);
    catch
      fl = {};
      warning('Could not list files in %s',dirW);
    end
    if numel(fl) ~= nSessFiles(ismember(labelSessID,sID{j}))
      warning(['Skipping: Number of expected log files not met '...
        'for %s %s (%d instead of %d).\n'...
        'Please handle this by adding to or removing from these files.\n'...
        'Attention: You might need to modify physio acq/mat as well.\n'...
        'Please apply all dataset modifications by adding appropriate '...
        'code to the function "func_correctDataset":'],...
        pID,sID{j},numel(fl),nSessFiles(j));
      fprintf('%d files found for %s %s:\n',numel(fl),pID,sID{j})
      disp(fl);
      fprintf('\n')
      continue
    end
    dirTmp = fullfile(dirW,'temp');
    if isfolder(dirTmp)
      rmdir(dirTmp,'s');
    end
    mkdir(dirTmp);
    cellfun(@(x) copyfile(x,dirTmp),fl);
    which ter_readPresentationLog
    [~,fn_out] = ter_readPresentationLog('source',dirTmp,'target',dirTmp);
    movefile(fn_out,fn_mat);
    rmdir(dirTmp,'s');
  end
end

for i=1:size(ptab,1)
  pID = ptab.participant_id{i};
  sID = struct2cell(dir(fullfile(fp_s,pID,'ses-*')));
  sID = sID(1,:);
  if isempty(sID)
    sID = {''};
  end
  for j=1:numel(sID)
    dirW = fullfile(fp_s,pID,sID{j},subdirLogs);
    fn_mat = strrep(fullfile(dirW,...
      [pID '_' sID{j} '_logdataQuestionnaires.mat']),'__','_');
    if exist(fn_mat,'file')==2
      continue % skip existing sets
    end
    % find logs to read
    try
      fl = ter_listFiles(dirW,id_questlogs,0);
    catch
      fl = {};
      warning('Could not list files in %s',dirW);
    end
    if numel(fl) ~= nQuestFiles(ismember(labelSessID,sID{j}))
      warning(['Skipping: Number of expected questionnaire files not met '...
        'for %s %s (%d instead of %d).\n'...
        'Please handle this by adding to or removing from these files.\n'...
        'Attention: You might need to modify physio acq/mat as well.\n'...
        'Please apply all dataset modifications by adding appropriate '...
        'code to the function "func_correctDataset":'],...
        pID,sID{j},numel(fl),nQuestFiles(j));
      fprintf('%d files found for %s %s:\n',numel(fl),pID,sID{j})
      disp(fl);
      continue
    end
    fl2 = ter_listFiles(dirW,id_questresp,0);
    if numel(fl2) ~= nQuestFiles(j)
      warning(['Skipping: Number of expected questionnaire response '...
        'files not met for %s (%d instead of %d).\n'...
        'Please handle this by adding to or removing from these files.\n'...
        'Attention: You might need to modify physio acq/mat as well.\n'...
        'Please apply all dataset modifications by adding appropriate '...
        'code to the function "func_correctDataset":'],...
        pID,sID{j},numel(fl2),nQuestFiles(j));
      fprintf('%d files found for %s %s:\n',numel(fl2),pID,sID{j})
      disp(fl2);
      continue
    end
    dirTmp = fullfile(dirW,'temp');
    if isfolder(dirTmp)
      rmdir(dirTmp,'s');
    end
    mkdir(dirTmp);
    cellfun(@(x) copyfile(x,dirTmp),[fl;fl2]);
    try
      [~,fn_out] = ter_readPresentationLog('source',dirTmp,'target',dirTmp);
      movefile(fn_out,fn_mat);
    catch
      warning('cloud not work here:%s', dirTmp);
    end
    
    
    rmdir(dirTmp,'s');
  end
end



