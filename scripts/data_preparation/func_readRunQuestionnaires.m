function func_readRunQuestionnaires(fp_s,fp_d)

  % if run without inputs try to gather inputs from base workspace
  if nargin<1
    fp_s = evalin('base','fp_s');
  end
  if nargin<2
    fp_d = evalin('base','fp_d');
  end

  fprintf('\n%s\n%s - reading digital questionnaires data to tables \n',...
    repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))

  subdirLogs  = 'stimlogs';
  suffix_mat  = '_logdataQuestionnaires.mat';
  
  protocol_labels = {
    'after acquisition (ses-1_run-2)'
    'after extinction (ses-2_run-1)' 
    'after habituation (ses-1_run-1)'
    'after recall (ses-3_run-1)'      
    'after volatile (ses-3_run-2)'   };
%     'Day 1: After Acquisition (Brain Stimulation)'
%     'Day 1: After Extinction (Brain Stimulation)' 
%     'Day 1: After Habituation (Fear)'             
%     'Day 1: After Acquisition (Fear)'                 
%     'Day 1: After Extinction (Fear)'                  
%     'Day 2: After Recall (Fear)'                  };
   
   
  %protocols_bs = protocol_labels (1:2);
  protocols_fc = protocol_labels;
  protocols_temporalOrder = protocol_labels([3,1,2,4,5]);
  
  %offset_lickertScale = 0; %lickert scale in questionnaires scales from 0 to 8 instead of 1 to 9

  %fn_bs_tsv = fullfile(fp_d,'phenotype','questions_brainstimulation.tsv');
  fn_fc_tsv = fullfile(fp_d,'phenotype','questions_fearconditioning.tsv');
  
  fn_ptab = fullfile(fp_d,'participants.tsv');
  %tsv_param = {'filetype','text','delim','tab','treat','n/a'};
  %ptab = readtable(fn_ptab,tsv_param{:});
  ptab = essbids_readTsv(fn_ptab);


  fn2ignore = {};

  %% first loop over all protocol files and read to table tab0
  tab0 = [];
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    fl = essbids_listFiles(fullfile(fp_s,pID,'**',...
      [pID '*task-fear*_questanswers.txt']));
    try 
      bf = struct2table(cellfun(@essbids_parseLabel, fl));
    catch
      fprintf('some of these filenames do not fit the pattern expected');
      disp(fl);
      continue
    end
    bf = bf(not(ismember(bf.fname,fn2ignore)),:);
    if numel(unique(bf.acq))<size(bf,1)
      fprintf(['some of these filenames do not fit the pattern expected\n'...
      'consider adding filenames to fn2ignore']);
      disp(fl);
      continue
    end
    for j=1:size(bf,1)
      fn1 = fullfile(bf.fpath{j},[bf.fname{j} bf.extension{j}]);
      tmp = readAnswerFile(fn1);
      tmp = tmp(not(ismember(tmp.Result,'Info')),:);
      str2rep = {
        'Ja'      'yes'
        'Nein'    'no'
        'Quadrat' 'square'
        'Raute'   'diamond'};
      for k=1:size(str2rep)
        tmp.Result = strrep(tmp.Result,str2rep{k,:});
      end
      str2add = sprintf('ses%s_run%s_%s_',...
        bf.ses{j},bf.run{j},bf.acq{j}(1:3));
      tmp.Label = strcat(str2add,tmp.Label);
      ind2rm = false(size(tmp,1),1);
      for l=1:size(tmp,1)-1
        if tmp.TrialN(l)>= min(tmp.TrialN(l+1:end))
          ind2rm(l)=true;
        end
      end
      tmp = tmp(~ind2rm,:);
      if isempty(tab0)
        tab0 = tmp;
      else
        tab0 = [tab0;tmp]; %#ok<AGROW>
      end
    end

    ses = struct2cell(dir(fullfile(fp_d,pID,'ses-*')));
    ses = ses(1,:);
    if isempty(ses)
      ses = {''};
    end
    for j=1:numel(ses)
      sID = ses{j};
      fn = fullfile(fp_s,pID,sID,subdirLogs,...
        strrep([pID '_' sID suffix_mat],'__','_'));
      if exist(fn,'file')~=2
        warning('file does not exist: %s',fn)
        break
        %continue
      end
      clearvars logdata;
      load(fn,'logdata');
      fin = fieldnames(logdata);
      for k=1:numel(fin)
        
        if ~isfield(logdata.(fin{k}),'questList')
          continue
        end
        pID2 = logdata.(fin{k}).('participant_id');
        if ~isequal(strfind(pID2,'sub-'),1)
          pID2 = ['sub-' pID2]; %#ok<AGROW>
        end
        if ~isequal(pID,pID2)
          warning('participant ID mismatch: %s %s %s',pID,pID2,fn);
          continue
        end
        tmp = logdata.(fin{k}).('questList');
        tmp = tmp(~ismember(lower(tmp.Result),{'info','skipped',''}),:);
        tmp = tmp(~contains(lower(tmp.Label),...
          {'comments','trainingtrial'}),:);
        tmp.participant_id = repmat({pID},size(tmp,1),1);
        tmp.session = repmat({sID},size(tmp,1),1);
        tmp.time =  repmat({logdata.(fin{k}).('datetime')},size(tmp,1),1);
        tmp.protocol = repmat({logdata.(fin{k}).('protocol')},...
          size(tmp,1),1);
        
        tmp = tmp(:,[end-3:end,1:end-4]);
        if iscell(tmp.TrialN) 
          tmp.TrialN= str2double(tmp.TrialN);
        end
        %removefirst answers to repeated questions and questions skipped
        %after step back
        ind2rm = false(size(tmp,1),1);
        for l=1:size(tmp,1)-1
          if tmp.TrialN(l)>= min(tmp.TrialN(l+1:end))
            ind2rm(l)=true;
          end
        end
        tmp = tmp(~ind2rm,:);
        if isempty(tab0)
          tab0 = tmp;
        else
          tab0 = [tab0;tmp]; %#ok<AGROW>
        end
      end
    end
  end
  tab0.protocol = strrep(tab0.protocol,'ses-3_run1','ses-3_run-1');
  str2rep = {
    'Ja'      'yes'
    'Nein'    'no'
    'Quadrat' 'square'
    'Raute'   'diamond'};
  for i=1:size(str2rep,1)
    tab0.Result = strrep(tab0.Result,str2rep{i,:});
  end
  for i=1:size(tab0,1)
    tmp = regexp(tab0.protocol{i},'run-[0-9]*','match');
    try
      tab0.run{i} = tmp{1};
    catch
      tab0.run{i} = '';
    end
  end
  CS = {'Diamond','Square'};
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    csp = ptab.cs_plus{i};
    csp(1) = upper(csp(1));
    csm = CS{not(ismember(CS,csp))};
    ind = ismember(tab0.participant_id,pID) & contains(tab0.Label,csp);
    tab0.Label(ind) = strrep(tab0.Label(ind),csp,'_CSplus');
    ind = ismember(tab0.participant_id,pID) & contains(tab0.Label,csm);
    tab0.Label(ind) = strrep(tab0.Label(ind),csm,'_CSminus');
  end
  for i=1:size(tab0,1)
    try
      tab0.Label{i} = strrep(sprintf('%s_%s_%s_%s',tab0.session{i},...
        tab0.run{i},tab0.protocol{i}(7:9),tab0.Label{i}),'-','');
    catch
      tab0.Label{i} = strrep(sprintf('%s_%s_%s',tab0.session{i},...
        tab0.run{i},tab0.Label{i}),'-','');
    end
  end
  t0 = table;
  t0.participant_id = ptab.participant_id;
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    ind = find(ismember(tab0.participant_id,pID));
    for j=1:numel(ind)
      t0.(tab0.Label{ind(j)})(i) = tab0.Result(ind(j));
    end
  end
  vn = t0.Properties.VariableNames;
  ind0 = find(contains(vn,'ForcedChoiceContin'));
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    ind = find(ismember(t0.participant_id,pID));
    if strcmpi(ptab.cs_plus{i},'diamond')
      for j=1:numel(ind0)
        if strcmpi(t0{ind,ind0(j)},'diamond')
          t0{ind,ind0(j)} = {'CSplus'};
        elseif strcmpi(t0{ind,ind0(j)},'square')
          t0{ind,ind0(j)} = {'CSminus'};
        end
      end
    elseif strcmpi(ptab.cs_plus{i},'square')
      for j=1:numel(ind0)
        if strcmpi(t0{ind,ind0(j)},'diamond')
          t0{ind,ind0(j)} = {'CSminus'};
        elseif strcmpi(t0{ind,ind0(j)},'square')
          t0{ind,ind0(j)} = {'CSplus'};
        end
      end
    end
  end
      
  % reorder 
  sID = unique(tab0.session);
  for i=1:numel(sID)
    rID = unique(tab0.run(ismember(tab0.session,sID(i))));
    for j=1:numel(rID)
      label = strrep(sprintf('%s_%s_',sID{i},rID{j}),'-','');
      vn = t0.Properties.VariableNames;
      ql = {'Valence','Arousal','Fear','Expectation'};
      ind = find(contains(vn,label) & contains(vn,ql));
      tmp = sortrows(table(vn(ind)',ind'));
      tmp2 = t0(:,tmp.Var2);
      t0(:,ind)=[];
      t0 = [t0(:,1:(min(ind)-1)) tmp2 t0(:,min(ind):end)];
    end
  end
  
  
  
  
  
  fn_temp = fullfile(fileparts(fn_fc_tsv),'temp.tsv');
  ter_writeBidsTsv(t0,fn_temp);
  tab_fc = ter_readBidsTsv(fn_temp);
  delete(fn_temp);
  
  
%   for i=1:numel(tab0.Result)
%     if numel(regexp(tab0.Result{i},'\d*')) == numel(tab0.Result{i})
%       tab0.Result{i} = str2double(tab0.Result{i});
%     end
%   end
  
  %% double check questionniare presence
  pIDs = unique(tab0.participant_id);
  for i=1:numel(pIDs)
    pID = pIDs{i};
    ind = ismember(tab0.participant_id,pID);
    tmp = unique(tab0(ind,1:4));
    %tmp = unique(tab0(ind,1:3));
    ind2 = ismember(protocol_labels,tmp.protocol);
    if prod(ind2)==0
      warning('protocols missing for %s :',pID)
      disp(protocol_labels(~ind2));
      disp(tmp)
    elseif ~isequal(tmp.protocol,protocols_temporalOrder)
      warning('temporal order messed up for %s :',pID)
      tmp2 = table(tmp.protocol,protocols_temporalOrder);
      tmp2.Properties.VariableNames = {'detected_order','should_be'};
      disp(tmp2);
      disp(tmp)
    end
  end
  
  
  
  
  
  
%     
%   %% now do the same for fear conditioning questionnaires
%   for i=1:size(tab0,1)
%     tmp = regexp(tab0.protocol{i},'run-[0-9]*','match');
%     try
%       tab0.run{i} = tmp{1};
%     catch
%       tab0.run{i} = '';
%     end
%   end
%   for i=1:size(tab0,1)
%     mylabels{i} = regexp(tab0.protocol{i},'run-[0-9]*','match');
%     try
%       tab0.run{i} = tmp{1};
%     catch
%     end
%   end
%   
%   
%   if iscell(tab0.run{1})
%     tab0.run = cellfun(@(x) x{1},tab0.run,'uni',0);
%   end
%   for i=1:numel(protocols_fc)
%     ind =ismember( tab0.protocol,protocols_fc(i));
%     mylabels = unique(cellfun(@(x,y) sprintf('%02.0f_%s',x,y),...
%       num2cell(tab0.TrialN(ind)),tab0.Label(ind),'uni',0));
%     mylabels0 = cellfun(@(x) x(4:end), mylabels,'uni',0);
%     for j=1:numel(mylabels0)
%       mylabels{j} = sprintf('ses%d_q%02.0f_%s',i,j,mylabels0{j});
%       ind2 = ismember(tab0.Label,mylabels0{j});
%       tab0.question(ind & ind2) = mylabels(j);
%     end    
%   end
%   ind = ismember( tab0.protocol,protocols_fc);
%   tab1 = tab0(ind,:);
%   pIDs = unique(tab1.participant_id);
%   qIDs = unique(tab1.question);
%   tab_fc = [table(pIDs),array2table(nan(numel(pIDs),numel(qIDs)))];
%   tab_fc.Properties.VariableNames = ['participant_id'; qIDs];
%   for i=1:numel(pIDs)
%     pID = pIDs{i};
%     for j=1:numel(qIDs)
%       qID = qIDs{j};
%       ind1 = find(ismember(tab1.participant_id,pID) & ...
%         ismember(tab1.question,qID));
%       % questions can be answered twice or more, use the very last answer
%       % only
%       if numel(ind1)==0
%         continue
%       end
%       ind1 = ind1(end);
%       if numel(ind1)==1
%         if contains(qID,{'StimReceived' 'StimNumber' 'PercStimAfter' 'Contingency'})
%           tab_fc.(qID)(i) = tab1.Result(ind1);
%         else
%           tab_fc.(qID)(i) = tab1.Result(ind1) + offset_lickertScale ;
%         end
%       else
%         warning('mismatch')
%       end
%     end
%   end
 
  % now check wether information is different from previous table, if it is
  % write to phenotype folder
  if exist(fn_fc_tsv,'file')==2
    %tab_fc0 = readtable(fn_fc_tsv,tsv_param{:});
    tab_fc0 = ter_readBidsTsv(fn_fc_tsv);
  else
    tab_fc0 = [];
  end
  if ~isequaln(tab_fc,tab_fc0)
    if ~isfolder(fileparts(fn_fc_tsv))
      mkdir(fileparts(fn_fc_tsv));
    end
    ter_writeBidsTsv(tab_fc,fn_fc_tsv);
%     writetable(tab_fc,fn_fc_tsv,tsv_param{1:4});
%     % rename NaN to n/a to comply with BIDS standard
%     txt = fileread(fn_fc_tsv);
%     txt2 = strrep(txt,'NaN','n/a');
%     if ~isequal(txt,txt2)
%       fid = fopen(fn_fc_tsv,'w');
%       fprintf(fid,txt2);
%       fclose(fid);
%     end
    fprintf('Updated file: %s\n',fn_fc_tsv);
  else
    fprintf('Questionnaire responses are up to date.\n');
  end
  
  fn_j = strrep(fn_fc_tsv,'.tsv','.json');
  if exist(fn_j,'file')~=2
    vn = tab_fc.Properties.VariableNames;
    clearvars jinfo;
    for j=2:numel(vn)
      jinfo.(vn{j}).Description = vn{j};
    end
    ter_savePrettyJson(fn_j,jinfo);
    warning('Writting dummy json file, please edit manually.');
  end
  
  
end

    
function t0 = readAnswerFile(fn_answers)

  txt = fileread(fn_answers);
  txt = strrep(txt,'\r','');
  txt2 = txt(strfind(txt,'TrialN'):end);
  tmp = strsplit(txt2,{'\t','\n','\r'});
  tmp = tmp(1:floor(numel(tmp)/3)*3);
  tmp = strrep(tmp,'\t','');
  tmp = reshape(tmp,[3,numel(tmp)/3])';
  t0 = cell2table(tmp(2:end,:),'VariableNames',tmp(1,:));
  t0 = addprop(t0,{'participant_id','time','protocol'},...
    {'table','table','table'});
  tmp = regexp(txt,'ParticipantID *: *(?<pID>[a-zA-Z0-9/\.\\\-_]*)',...
    'names');
  t0.Properties.CustomProperties.participant_id = tmp.pID;
  tmp = regexp(txt,...
    'Protocol filename *: *(?<protocol>[a-zA-Z0-9/\.\\\-_]*)','names');
  t0.Properties.CustomProperties.protocol = tmp.protocol;
  tmp = regexp(txt,'time *: *(?<time>[a-zA-Z0-9/\.\\\-_\ \:]*)','names');
  t0.Properties.CustomProperties.time = tmp.time;
  t0.TrialN = str2double(t0.TrialN);
end

