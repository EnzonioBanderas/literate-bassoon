function out = ter_bidsPhysio2cardioEval(dir_bids,dir_out)
  
  %% checking input 
  if iscell(dir_bids)
    if nargin < 2
      out = cellfun(@ter_bids2cardioEval,dir_bids);
    elseif iscell(dir_out)
      out = cellfun(@ter_bids2cardioEval,dir_bids,dir_out);
    else
      out = cellfun(@(x) ter_bids2cardioEval(x,dir_out),dir_bids);
    end
    return
  end
  if ~isfolder(dir_bids)
    disp(dir_bids)
    error('bids folder innvalid')
  end
  if nargin < 2
    dir_out = '';
  elseif iscell(dir_out)
    if numel(dir_out)==1
      out = ter_bids2cardioEval(dir_bids,dir_out{1});
      return
    else
      error('dir_out is a cell array of inappropriate size');
    end
  end
  
  
  %% ćreate target folder if neccessary
  if isempty(dir_out)
    dir_out = fullfile(fileparts(dir_bids),'derivatives','cardioEval');
  end
  if ~isfolder(dir_out) 
    if isfolder(fileparts(dir_out))
      mkdir(dir_out);
    else
      error('Target folder does not exist: %s\n',dir_out);
    end
  end
    
  
  fl_gz  = ter_listFiles(dir_bids,'sub-*_physio.tsv.gz');
  fl_tsv = ter_listFiles(dir_bids,'sub-*_physio.tsv');
  fl_gz = fl_gz(~contains(fl_gz,'task-rest'));
  fl_tsv = fl_gz(~contains(fl_tsv,'task-rest'));
  fl_tmp = cellfun(@(x) [x '.gz'],fl_tsv,'uni',0);
 
  
  fl2work = [fl_gz(~ismember(fl_gz,fl_tmp));fl_tsv];
  
  for i=1:numel(fl2work)
    fn_tsv = fl2work{i};
    [~,fn,fe] = fileparts(fn_tsv);
    
    pID = regexp(fn,'sub-[a-zA-Z0-9]*','match');
    ses = regexp(fn,'ses-[a-zA-Z0-9]*','match');
    pID = pID{1};
    if isempty(ses)
      ses = '';
    else
      ses = ses{1};
    end
    fp_eval = fullfile(dir_out,pID,ses);
    if ~isfolder(fp_eval)
      mkdir(fp_eval)
    end
    if strcmpi(fe,'.gz')
      fn_cardio = fullfile(fp_eval,strrep(fn,'physio.tsv','cardio.mat'));
      fn_json   = strrep(fn_tsv,'.tsv.gz','.json');
      fn_events = strrep(fn_tsv,'physio.tsv.gz','events.tsv');
      %gunzip(fn_tsv,fp_eval)
      %phystab = readtable(fullfile(fp_eval,fn),'filetype','text',...
      %  'delim','tab');
      %delete(fullfile(fp_eval,fn));
    else
      fn_cardio = fullfile(fp_eval,strrep(fn,'physio','cardio.mat'));
      fn_json   = strrep(fn_tsv,'.tsv','.json');
      fn_events = strrep(fn_tsv,'physio.tsv','events.tsv');
      %phystab = readtable(fn_tsv,'filetype','text',...
      %  'delim','tab');
    end
    fn_et = strrep(fn_cardio, '_cardio.mat','_eventTable.mat');
    if exist(fn_et,'file')==2 && exist(fn_cardio,'file')==2
      continue
    end
    if exist(fn_json,'file')~=2
      error('json file does not exist: %s',fn_json);
    elseif exist(fn_events,'file')~=2
      error('event file does not exist: %s',fn_events);
    end
    % now load physio data
    if strcmpi(fe,'.gz')
      gunzip(fn_tsv,fp_eval)
      phystab = readtable(fullfile(fp_eval,fn),'filetype','text',...
        'delim','tab');
      delete(fullfile(fp_eval,fn));
    else
      phystab = readtable(fn_tsv,'filetype','text',...
        'delim','tab');
    end
    evtab = readtable(fn_events,'filetype','text','delim','tab');
    jinfo = jsondecode(fileread(fn_json));
    phystab.Properties.VariableNames = jinfo.Columns;
    for j=1:numel(jinfo.Columns)
      colID = phystab.Properties.VariableNames{j} ;
      if isfield(jinfo,colID)
        try 
          phystab.Properties.VariableUnits{j} = jinfo.(colID).Units;
        catch
        end
      end
    end
    
        
    % select only the cardio relevant colums
    ind = contains(lower(jinfo.Columns),{'cardiac','pulseox','ecg','ekg'});
    phystab = phystab(:,ind);
    phystab.Time = 1000 * (jinfo.StartTime + ...
      (0:size(phystab)-1)/jinfo.SamplingFrequency)';
    phystab = phystab(:,[end,1:end-1]);
    phystab.Properties.VariableUnits{1} = 'ms';
    
    
    vn = phystab.Properties.VariableNames;
    phystab.Properties.VariableNames{ismember(vn,'cardiac')} = 'ECG';
    phystab.Properties.VariableNames{ismember(vn,'pulseoximeter')} = 'PPG';
    cardiodata = phystab;
    save(fn_cardio,'cardiodata','-v7.3');
    
    %% now we are getting very specific for the exüperiment in question
    durTIR = 4; %seconds
    
    event_label     = {};
    trial_index     = [];
    edaeval_index   = 0;
    evalBlock_start = [];
    evalBlock_stop  = [];
    event_onsets     = [];
    
    %#ok<*AGROW>
    for j=1:size(evtab,1)
      label0 = evtab.trial_type{j};
      if 1<evtab.duration(j) && evtab.duration(j)<6
        % treat prolonged, but still short single events as FIR and TIR only
        fulldur = round(evtab.duration(j),3);
        %FIR
        event_onsets(end+1)    = evtab.onset(j) ;
        evalBlock_start(end+1) = evtab.onset(j) ;
        evalBlock_stop(end+1)  = evtab.onset(j) + fulldur;
        trial_index(end+1)      = j;
        edaeval_index(end+1)    = edaeval_index(end) + 1;
        event_label{end+1}      = ['FIR_' label0];
        %TIR
        event_onsets(end+1)    = evtab.onset(j) ;
        evalBlock_start(end+1) = evtab.onset(j) + fulldur;
        evalBlock_stop(end+1)  = evtab.onset(j) + fulldur + durTIR;
        trial_index(end+1)      = j;
        edaeval_index(end+1)    = edaeval_index(end) + 1;
        event_label{end+1}      = ['TIR_' label0];
      elseif 6<evtab.duration(j) && evtab.duration(j)<20
        % long events cut into first, second and thirs interval response
        fulldur = round(evtab.duration(j),3);
        halfdur = round(evtab.duration(j)/2,3);
        % FIR
        event_onsets(end+1)    = evtab.onset(j) ;
        evalBlock_start(end+1) = evtab.onset(j);
        evalBlock_stop(end+1)  = evtab.onset(j) + halfdur;
        trial_index(end+1)     = j;
        edaeval_index(end+1)   = edaeval_index(end) + 1;
        event_label{end+1}     = ['FIR_' label0];
        % SIR
        event_onsets(end+1)    = evtab.onset(j) ;
        evalBlock_start(end+1) = evtab.onset(j) + halfdur ;
        evalBlock_stop(end+1)  = evtab.onset(j) + fulldur ;
        trial_index(end+1)     = j;
        edaeval_index(end+1)   = edaeval_index(end) + 1;
        event_label{end+1}     = ['SIR_' label0];
        % TIR
        event_onsets(end+1)    = evtab.onset(j) ;
        evalBlock_start(end+1) = evtab.onset(j) + fulldur;
        evalBlock_stop(end+1)  = evtab.onset(j) + fulldur + durTIR;
        trial_index(end+1)     = j;
        edaeval_index(end+1)   = edaeval_index(end) + 1;
        % next few lines are specific for the frea conditioning case here,
        % markup TIR depending on apearance or omission of US
        if j == size(evtab,1)
          event_label{end+1} = ['TIR-noUS_' label0];
        elseif strcmpi(evtab.trial_type{j+1}(1:2),'US')
          event_label{end+1} = ['TIR-US_' label0];   
        else
          event_label{end+1} = ['TIR-noUS_' label0]; 
        end
      else
        % treating single events as TIR equivalent, as they usually are
        % shocks
        addTIR = false;
        if j==1 
          addTIR = true;
        elseif diff(evtab.onset([j-1,j]))+diff(evtab.duration([j-1,j]))>.01 
          addTIR = true;
        end
        if addTIR
          event_onsets(end+1)    = evtab.onset(j);
          evalBlock_start(end+1) = evtab.onset(j);
          evalBlock_stop(end+1)  = evtab.onset(j) + durTIR;
          trial_index(end+1)     = j;
          edaeval_index(end+1)   = edaeval_index(end) + 1;
          event_label{end+1}     = ['TIR_' label0];
        end
      end
    end
    event_onsets              = event_onsets';
    trial_index               = trial_index';
    edaeval_index             = edaeval_index(2:end)';
    event_label               = event_label';
    evalBlock_start           = evalBlock_start';
    evalBlock_stop            = evalBlock_stop'; 
    eventtable = table(event_onsets,trial_index,edaeval_index,...
      event_label,evalBlock_start,evalBlock_stop);

    save(fn_et,'eventtable','-v7.3');
    
  end
    
  
  
  
  
end
  