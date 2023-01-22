function ptab = func_updatePtab(fp_s,fp_d)
% reads presentation log files 
% inputs : fp_s  -  path to sourcedata folder
%          fp_d  -  path to rawdata folder

fprintf('\n%s\n%s - Checking for ptab update\n',...
  repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))

fn_ptab = fullfile(fp_d,'participants.tsv');
if exist(fn_ptab,'file')==2
  ptab = ter_readptab(fn_ptab);
else
  ptab = table('Size',[0,3],'VariableTypes',{'cell','cell','double'},...
    'VariableNames',{'participant_id','sex','age'});
end

ptab0 = ptab;

dl = struct2table(dir(fullfile(fp_s,'sub-*')));
dl = dl(dl.isdir,:);
dl = dl.name;

pID_list = {};
for i=1:numel(dl)
  tmp = regexp(dl{i},'sub-[0-9a-zA-Z]*$','match');
  if isempty(tmp)
    continue
  end
  pID = dl{i};
  pID_list = [pID_list;pID]; %#ok<AGROW>
  if ismember(pID,ptab.participant_id)
    continue
  end
  fprintf('adding participant %s\n',pID);
  ptab.participant_id{end+1}=pID;
  for j = 2:numel(ptab.Properties.VariableNames)
    vn = ptab.Properties.VariableNames{j};
    if iscell(ptab.(vn)) 
      ptab.(vn){end} = 'n/a';
    elseif isnumeric(ptab.(vn))
      ptab.(vn)(end) = nan;
    else
      warning('neither cell nor numeric column: %s',vn);
    end
  end
end

ind_pID2rm = not(ismember(ptab.participant_id,pID_list));
if sum(ind_pID2rm)>0
  fprintf('removing these %d participants:\n',sum(ind_pID2rm));
end
ptab = ptab(not(ind_pID2rm),:);


%% get some information from demographics questionnaire
fn_demo = fullfile(fp_d,'phenotype','demographics.tsv');
if exist(fn_demo,'file')
  try
    t_demo = ter_readBidsTsv(fn_demo);
    for i=1:size(ptab,1)
      pID = ptab.participant_id{i};
      ind = ismember(t_demo.participant_id,pID);
      if sum(ind==1)
        ptab.age(i) = t_demo.age(ind);
        ptab.sex(i) = t_demo.sex(ind);
      end
    end
  catch
    warning(['Trouble merging particpants.tsv information with'...
      ' demographics.tsv']);
  end
end

%% get information from randomization list
fn_rl = fullfile(fp_s,'documents','randomlist.xlsx');
if exist(fn_rl,'file')==2
  tmp = warning('query','MATLAB:table:ModifiedAndSavedVarnames');
  warning('off','MATLAB:table:ModifiedAndSavedVarnames');
  t0 = readtable(fn_rl);
  warning(tmp.state,'MATLAB:table:ModifiedAndSavedVarnames');
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    ind = ismember(t0.SFBCode,strrep(pID,'sub-',''));
    if sum(ind)>1
      warning('same participant multiple times in random list: %s',pID)
      continue
    elseif sum(ind)==1
      if isnumeric(ptab.code_7t)
        ptab.code_7t = repmat({'n/a'},size(ptab.code_7t));
      end
      ptab.code_7t(i) = t0.ELHCode(ind);
      ptab.drug(i) = t0.Med_Group(ind);
      ptab.cs_sequence(i) = t0.Hab_Ext_Recall(ind);
    else
      try 
        ptab.code_7t{i} = 'n/a';
      catch
        ptab.code_7t(i) = nan;
      end
      ptab.drug(i) = nan;
      ptab.cs_sequence(i) = nan;
      fprintf('Participant not listed in randomlist: %s\n',pID);
    end
  end
end
rl_pIDs = unique(t0.SFBCode);
rl_pIDs = rl_pIDs(not(cellfun(@isempty,rl_pIDs)));
rl_pIDs = cellfun(@(x) ['sub-' x],rl_pIDs,'uni',0);
ind_rl = not(ismember(rl_pIDs,ptab.participant_id));
if sum(ind_rl)>0
  fprintf(['%d participant listed in randomlist but not '...
    'participants.tsv:\n'],sum(ind_rl));
  disp(rl_pIDs(ind_rl))
end

try
  ptab = verify_cs_sequence(fp_d,ptab);
catch
  warning('trouble verifying cs sequence');
end

%% get stim intensity
for i=1:size(ptab,1)
  pID = ptab.participant_id{i};
  fl = ter_listFiles(fullfile(fp_s,pID),...
    [pID '*_stimsettingratings.txt'],2);
  if numel(fl)>1
    warning('Multiple stim settings found, ignoring: %s',pID);
  elseif numel(fl)==1
    try
      t0 = readtable(fl{1},'filetype','text','delim','tab');
      ind = find(ismember(t0.Rating,'not sensed'),1,'last');
      if isempty(ind)
        ind = 0;
      end
      ptab.stim_sensory_threshold(i) = ...
        str2double(strrep(t0.StimIntensity(ind+1),' mA',''));
    catch
      warning('could not read stim sensory threshold for %s',pID')
    end
    try 
      txt = fileread(fl{1});
      tmp = regexp(txt,...
        'final stimulation intensity: (?<si>[0-9\.]*) mA','names');
      if not(isempty(tmp))
        ptab.stim_intensity(i) = str2double(tmp(end).si);
      else
        warning('could not read final stim intensity for %s',pID');
      end
    catch
      warning('could not read final stim intensity for %s',pID');
    end
    try
      ptab.stim_intensity(ptab.stim_intensity==0)=nan;
      ptab.stim_sensory_threshold(ptab.stim_sensory_threshold==0)=nan;
    catch
    end
  end
end


%% get datetimes of each session
for i=1:size(ptab,1)
  pID = ptab.participant_id{i};
  fl = ter_listFiles(fullfile(fp_s,pID),...
    [pID '*_run-1*_stimlog.txt'],2);
  try 
    t1 = struct2table(cellfun(@ter_parseFname,fl));
    if size(t1,1)==1
      t1.ses = {t1.ses};
      t1.datetime = {t1.datetime};
    end
    if numel(unique(t1.ses)) ~= size(t1,1) && size(t1,1)>1
      warning('double stimlogs : %s',pID);
      disp(fl);
      continue
    end
    for j=1:size(t1,1)
      fn = ['datetime_' strrep(t1.ses{j},'-','')];
      if isdatetime(ptab.(fn))
       ptab.(fn)(i) = datetime(strrep(t1.datetime{j},'datetime-',''),...
         'format','yyyy-MM-dd''T''HH:mm:ss');
       ptab.(fn).Format = 'yyyy-MM-dd''T''HH:mm:ss';
      elseif iscell(ptab.(fn))
        dt = strrep(t1.datetime{j},'datetime-','');
        try
          dt = datestr(datetime(dt),'yyyy-mm-ddTHH:MM:ss');
        catch
          dt = [ dt(1:4) '-'  dt(5:6) '-'  dt(7:11) ':' ... 
             dt(11:12) ':'  dt(13:14) ];
        end
        ptab.(fn){i} = dt;
      elseif isnumeric(ptab.(fn))
        ptab.(fn) = NaT(size(ptab.(fn)));
        ptab.(fn).Format = 'yyyy-MM-dd''T''HH:mm:ss';
        ptab.(fn)(i) = datetime(strrep(t1.datetime{j},'datetime-',''),...
          'format','yyyy-MM-dd''T''HH:mm:ss');
      end
    end
  catch
  end
end
try 
  ptab.datetime_ses1(cellfun(@isempty,ptab.datetime_ses1))={'n/a'};
  ptab.datetime_ses2(cellfun(@isempty,ptab.datetime_ses2))={'n/a'};
  ptab.datetime_ses3(cellfun(@isempty,ptab.datetime_ses3))={'n/a'};
catch
end


%% ensure correct order of fields analog to json file
fn_json = strrep(fn_ptab,'.tsv','.json');
if exist(fn_json,'file')==2
  jinfo = jsondecode(fileread(fn_json));
  fn = fieldnames(jinfo);
  if not(ismember('participant_id',fn))
    fn = ['participant_id';fn];
  end
  vn = ptab.Properties.VariableNames;
  ind_vn = nan(size(vn));
  count = 0;
  for i=1:numel(fn)
    tmp = find(ismember(vn,fn(i)));
    if not(isempty(tmp))
      count = count+1;
      ind_vn(count) = tmp;
    end
  end
end
ptab = ptab(:,ind_vn);

ptab = sortrows(ptab);

%% update tsv file if something is new
if not(isequaln(ptab0,ptab))
  fprintf('Updating ptab\n')
  ter_writeptab(ptab,fn_ptab);
else
  fprintf('ptab needed no update\n')
end

end

function ptab = verify_cs_sequence(fp_d,ptab)
  
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    fl = ter_listFiles(fp_d,[pID '*run-1*_events.tsv'],3); 
    % just so happens that habituation extinction and recall are the
    % respective runs 1, all should start either with CS+ (CS_sequence 1)
    % or CS- (CS sequence 2)
    myseq = ptab.cs_sequence(i);
    if isnan(myseq)
      myseq = [];
    end
    for j=1:numel(fl)
      t0 = ter_readBidsTsv(fl{j});
      ind_csp = find(contains(t0.trial_type,'CSplus'),1, 'first');
      ind_csm = find(contains(t0.trial_type,'CSminus'),1,'first');
      if isempty(ind_csp) || isempty(ind_csm)
        % nothing to do
      elseif ind_csm<ind_csp
        myseq = [myseq;2]; %#ok<AGROW>
      elseif ind_csm>ind_csp
        myseq = [myseq;1]; %#ok<AGROW>
      end    
    end
    if numel(unique(myseq))==1
      ptab.cs_sequence(i) = myseq(1);
    elseif numel(unique(myseq))>1
      warning('multiple cs sequences found: %s',pID)
      ptab.cs_sequence(i) = nan;
    else
      ptab.cs_sequence(i) = nan;
    end
  end
end