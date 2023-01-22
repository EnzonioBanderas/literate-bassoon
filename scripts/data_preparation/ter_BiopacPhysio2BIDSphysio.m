function [fn_out, tables_physio] = ter_BiopacPhysio2BIDSphysio(fn_biopacmat,fp_out,TR,addInfo,minDur,minRunSpacing)
% TR in seconds
% if multiple datasets are found they will be saved as run-1 .. run-n,
% independent of run information in addInfo
% will generate events table from all digital channels using channel labels
%
%
% fMRI experiments
%   Use TR=nan to automatically detect TR from trigger channel, use numeric
%   value in seconds for a fixed specific TR. Value "minDur" represents 
%   the minimum number of triggers in any given run. Runs shorter will be
%   trimed. minRunSpacing is not utilized and ignored.
%
% Purely behavioral experiments
%   Use TR=0 for purely behavioral experiments with no scan triggers. In 
%   this case and, if present, the "Trial" channel will be used to trim the
%   dataset, with no "trial" channel present everything will be presented 
%   in one block. 
%   minDur will represent the minimal number of trials per run. 
%   minRunSpacing represents the minimal temporal distance between 
%   "trials" in seconds to be considered a cutoff between runs. E.g. for
%   minRunSpacing=100 the first two trial trails more than 100 s apart will 
%   be considered the last trial of run-1 and the first trial of run-2 and
%   so on.  
%
%   2021-05-23 : run seperation for purely behavioral experiments
%                introduced
%
%   2020-12-18 : digital channels 2 events.tsv
%                use SamplingFrequency instead of SamplingRate
%

%   fn_out = {};
%   tables_physio= {};
   
  if nargin < 2
    tmp = strsplit(which(fn_biopacmat),'sourcedata');
    fp_out = fullfile(tmp{1},'rawdata');
  end
  if nargin < 3
    TR = nan;
  end
  if nargin <4
    addInfo = '';
  end
  if nargin <5
    minDur = 1;
  end
  if nargin <5
    minRunSpacing = 60; % 60 seconds
  end
  
  if ~isfolder(fp_out)
    error('invalid output folder specified')
  end
  
  tolerance = 5;  % trigger tolerance in ms
  timePre   = 20; % time before mri scan start in seconds 
  timePost  = 20; % time after mri scan stop in seconds 
  
  % assemble BIDs flags
  bids_flags = {'sub','ses','task','acq','dir','run'}; % run has to be the last
  
  for i=1:numel(bids_flags)
    tmp = regexp(addInfo,[bids_flags{i} '-[0-9a-zA-Z]*'],'match');
    if isempty(tmp)
      tmp = regexp(fn_biopacmat,[bids_flags{i} '-[0-9a-zA-Z]*'],'match');
    end
    if isempty(tmp)
      binfo.(bids_flags{i}) = '';
    else
      binfo.(bids_flags{i}) = tmp{end};
    end
  end
 
  
  % load and process biopac physio data
  ds = load(fn_biopacmat);
  
  labels = cellstr(ds.labels);
  units = cellstr(ds.units);
  bids_labels = labels;
  
  if strcmpi(ds.isi_units,'ms')
    jinfo.SamplingFrequency = round(1e3/ds.isi);
  else
    jinfo.SamplingFrequency = round(1/ds.isi,3);
  end
  
  labels2change = {
    'ecg',       'cardiac';
    'trigger',   'trigger';
    'ppg'        'pulseoximeter';
    'resp'       'respiratory'
    'eda'        'skinconductance'};
    
  labels_order = {'trigger','cardiac','pulseoximeter','respiratory',...
    'skinconductance'};

  for i=1:size(labels2change,1)
    bids_labels(contains(lower(bids_labels),labels2change{i,1})) = ...
      labels2change(i,2);
  end
  
  if ~ismember('trigger',bids_labels) && TR~=0
    error('cannot identify fmri trigger channel')
  end
  
  ind4physio = [];
  for i=1:numel(labels_order)
    ind4physio = [ind4physio,...
      find(ismember(bids_labels,labels_order(i)),1,'first')]; %#ok<AGROW>
  end
  t1 = array2table(ds.data(:,ind4physio));
  t2_labels = strrep(strrep(strrep(strrep(cellstr(ds.labels),' ',''),...
    '-','_'),',','_'),'.','_');
  t2 = array2table(ds.data,'VariableNames',t2_labels);
  t1.Properties.VariableNames = bids_labels(ind4physio);
  t1.Properties.VariableUnits = units(ind4physio);
  
  %identify digital channel 2 use for events table
  labels_other = bids_labels(~ismember(bids_labels,bids_labels(ind4physio)));
  labels_other = strrep(strrep(strrep(strrep(labels_other,' ',''),...
    '-','_'),',','_'),'.','_');
  labels4events = {};
  for i=1:numel(labels_other)
    if prod(ismember(t2.(labels_other{i}),[0,1,5]))==1
      labels4events = [labels4events;labels_other{i}]; %#ok<AGROW>
    end
  end
    
  jinfo.StartTime = nan;
  jinfo.Columns = t1.Properties.VariableNames;
  for i=2:numel(jinfo.Columns)
    jinfo.(jinfo.Columns{i}).Units = t1.Properties.VariableUnits{i};
  end
  
  % add some more info to json  
  jinfo.Manufacturer = 'BIOPAC';
  fn_biopacacq = strcat(fn_biopacmat(1:end-3),'acq');
  if exist(fn_biopacacq,'file')==2
    try 
      txt = fileread(fn_biopacacq);
      station = unique(regexp(txt,'MP[A-Za-Z0-9\-]*','match'));
      devices = unique(regexp(txt,' - [A-Za-Z0-9\-]*','match'));
      devices = strrep(devices,' - ','');
      mmn = station{1};
      for j=1:numel(devices)
        mmn = strcat(mmn,',',devices{j});
      end
      jinfo.ManufacturerModelName = mmn;
    catch
      fprintf('Cannot read acq file: %s\n',fn_biopacacq);
    end
  end
  
  
  % differtntiate further steps wether fmri data nor purely behavioral
  % dataset 
  if TR~=0
    fprintf('fMRI experiment \n');
    % now calculate blocks of trigger events
    t1.trigger = round(t1.trigger/5);
    t1.trigger = [false;diff(t1.trigger)==1];
    tdiff = diff(find(t1.trigger));

    if isnan(TR)
      TR_samples =  median(tdiff);
      TR = TR_samples/jinfo.SamplingFrequency;
      fprintf('auto detected TR to be %.3f seconds, please verify\n',TR)
    else
      TR_samples = TR*jinfo.SamplingFrequency;  
    end
    tol = round(tolerance/1e3*jinfo.SamplingFrequency);

    % clear away all non-fitting trigger inputs
    st = find(t1.trigger);
    st0 = st;
    for j=2:numel(st)-1
      tr1 = diff(st(j-1:j));
      tr2 = diff(st(j:j+1));
      if (tr1 < TR_samples-tol || tr1 > TR_samples+tol) &&  ...
          (tr2<TR_samples-tol || tr2>TR_samples+tol)
        t1.trigger(st(j))=false;
        %disp(st(j))
      end
    end
    st = find(t1.trigger);
    tr1 = diff(st(1:2));
    tr2 = diff(st(end-1:end));
    if (tr1 < TR_samples-tol || tr1 > TR_samples+tol) 
      t1.trigger(st(1))=false;
    end
    if (tr2 < TR_samples-tol || tr2 > TR_samples+tol) 
      t1.trigger(st(end))=false;
    end
    if numel(st0)<numel(st)
      fprintf('Cleared %d non-fitting triggers\n',numel(st)-numel(st0));
    end


    ind_last = [find([diff(st)>TR_samples+tol;false]);sum(t1.trigger)];
    ind_first = [1;ind_last(1:end-1)+1];
    ntrig = [ind_last(1);diff(ind_last)];

    if sum(ntrig < minDur)>0
      fprintf('ignoring %d scans with a total of %d trigger pulses\n',...
        sum(ntrig < minDur), sum(ntrig(ntrig < minDur)));
    end

    ind_first = ind_first(ntrig >= minDur);
    ind_last  = ind_last(ntrig >= minDur);
  else
    fprintf('Behavioral experiment \n');
    onset_first = size(t2,1);
    stop_last   = 1;
    if ismember("Trial",labels4events)
      % use trial channel to detect start and end points and to split into
      % separate runs
      dat = t2.Trial;
      dat = dat./max([1,max(dat)]); % normalize to 1, ;
      onset  = find([0;diff(dat)]==1);
      stops  = find([0;diff(dat)]==-1);
      if onset(1)>stops(1)
        stops = stops(2:end);
      end
      if onset(end)>stops(end)
        stops(end+1) = size(t2,1); %#ok<NASGU>
      end
      iti = (onset(2:end) - stops(1:end-1) )./jinfo.SamplingFrequency ;
      ind_runstarts = [1;1+find(iti >= minRunSpacing)];
      ind_runends   = [find(iti >= minRunSpacing); numel(stops)];
      nTrials = ind_runends - ind_runstarts;
      ind_first = onset(ind_runstarts(nTrials >= minDur));
      ind_last  = stops(ind_runends(nTrials >= minDur));
      indRun = 0;
      for i=1:numel(nTrials)
        if nTrials < minDur
          fprintf('ignoring brief run with only %d trials\n',nTrials(i));
        else
          indRun = indRun +1 ;
          fprintf('detected run-%d with %d trials\n',indRun,nTrials(i))
        end
      end
    else
      % no Trial channel, hence use all events to setect start end end of
      % one single large block
      for j=1:numel(labels4events)
        lbl = labels4events{j};
        dat = t2.(lbl);
        dat = dat./max([1,max(dat)]); % normalize to 1, ;
        onset  = [0;diff(dat)]==1;
        stops  = [0;diff(dat)]==-1;
        if numel(onset)<numel(stops)
          onset = [1;onset]; %#ok<AGROW>
        end
        if numel(stops)~=0 
          if stops(1) < onset(1)
            stops(1) = [];
          end
        end
        if numel(stops)<numel(onset)
          stops(end+1) = size(t2,1); %#ok<AGROW>
        end
        onset_first = min([onset_first,find(onset,1,'first')]);
        stop_last   = max([stop_last,find(stops,1,'last')]);  
      end
      ind_first = max([1,onset_first]);
      ind_last  = min([size(t2,1),stop_last]);
    end
  end
  
  dout = cell(size(ind_first));
  fn_out        = cell(size(dout));
  tables_physio = cell(size(dout));
  
  for i=1:numel(ind_first)
    if TR == 0
      ind_0 = max(1,ind_first(i)-timePre*jinfo.SamplingFrequency);
      ind_1 = min(size(t1,1),ind_last(i)+timePost*jinfo.SamplingFrequency);
      dout{i}.StartTime = (ind_0-ind_first(i))/jinfo.SamplingFrequency;
      subfolder_out = 'beh';
    else
      ind_0 = max(1,st(ind_first(i))-timePre*jinfo.SamplingFrequency);
      ind_1 = min(size(t1,1),st(ind_last(i))+timePost*jinfo.SamplingFrequency);
      dout{i}.StartTime = (ind_0-st(ind_first(i)))/jinfo.SamplingFrequency;
      subfolder_out = 'func';
    end
    dout{i}.data = t1(ind_0:ind_1,:);
    %suppress trigger leading up to acquisition (prescans, opposed phase,
    % SBref, etc)
    ind_pre = 1:(-1*dout{i}.StartTime*jinfo.SamplingFrequency);
    if isfield(dout{i}.data,'trigger')
      dout{i}.data.trigger(ind_pre) = 0;
    end
    tmp = (1:size(dout{i}.data,1))'/jinfo.SamplingFrequency+dout{i}.StartTime;
    tables_physio{i} = [ ...
      table(tmp,'VariableName',{'time'}) ...
      t1(ind_0:ind_1,1) t2(ind_0:ind_1,:)];
    if numel(dout) == 1
      dout{1}.fname = fullfile(fp_out,binfo.sub,binfo.ses,subfolder_out,...
        [binfo.sub '_' binfo.ses '_' binfo.task  '_' ...
        binfo.acq '_' binfo.dir '_' binfo.run '_physio']);
    else
      runinfo = sprintf('run-%d',i);
      dout{i}.fname = fullfile(fp_out,binfo.sub,binfo.ses,subfolder_out,...
        [binfo.sub '_' binfo.ses '_' binfo.task '_' ...
        binfo.acq '_' binfo.dir '_' runinfo '_physio']);
    end
    while contains(dout{i}.fname,'__')
      dout{i}.fname = strrep(dout{i}.fname,'__','_');
    end
  end
  
    
  for i=1:numel(dout)
    while contains(dout{i}.fname,'__')
      dout{i}.fname = strrep(dout{i}.fname,'__','_');
    end
  end
    
  for i=1:numel(dout)
    jinfo.StartTime = dout{i}.StartTime;
    if TR ~= 0
      fprintf('Writting %s\n with %d scans\n',dout{i}.fname,...
        sum(dout{i}.data.trigger));
    else
      fprintf('Writting %s\n',dout{i}.fname)
    end
    if ~isfolder(fileparts(dout{i}.fname))
      mkdir(fileparts(dout{i}.fname));
    end
    tmp = regexp(dout{i}.fname,'task-(?<task>[0-0a-zA-Z]*)','names');
    if not(isempty(tmp))
      jinfo.TaskName = tmp(end).task;
    else
        if isfield(jinfo, 'TaskName')
            jinfo = rmfield(jinfo,'TaskName');
        end
    end
    ter_savePrettyJson([dout{i}.fname '.json'],jinfo);
    
    
    writetable(dout{i}.data,[dout{i}.fname '.tsv'],'fileType','text',...
      'delim','tab','writeVar',0);
    txt = fileread([dout{i}.fname '.tsv']);
    txt2 = strrep(txt,'NaN','n/a');
    if ~isequal(txt,txt2)
      fid = fopen([dout{i}.fname '.tsv'],'w');
      fprintf(fid,'%s',txt2);
      fclose(fid);
    end
    system(sprintf('pigz %s',strrep([dout{i}.fname '.tsv'],' ','\ ')));
    fn_out{i} = dout{i}.fname;
  end
  
  
  for i=1:numel(tables_physio)
    fn_events = strrep(fn_out{i},'_physio','_events.tsv');
    tme = tables_physio{i}.time;
    et = table('size',[0,3],'VariableType',{'double','double','cell'},...
      'VariableNames',{'onset','duration','trial_type'});
    for j=1:numel(labels4events)
      lbl = labels4events{j};
      dat = tables_physio{i}.(lbl);
      dat = dat./max([1,max(dat)]); % normalize to 1, ;
      onset   = tme([0;diff(dat)]==1);
      stops    = tme([0;diff(dat)]==-1);
      if numel(onset)<numel(stops)
        onset = [tme(1);onset]; %#ok<AGROW>
      end
      if numel(stops)~=0 
        if stops(1) < onset(1)
          stops(1) = [];
        end
      end
      if numel(stops)<numel(onset)
        stops(end+1) = tme(end); %#ok<AGROW>
      end
      duration = stops-onset;
      trial_type = repmat({lbl},size(duration));
      et = [et;table(onset,duration,trial_type)]; %#ok<AGROW>
    end
    et = sortrows(et);
    if isempty(et)
      fprintf('No events recorded\n')
    else
      fprintf('%d events saved to %s\n',size(et,1),fn_events);
      digits2use =  ceil(log(jinfo.SamplingFrequency)/log(10));
      myExpr = sprintf('%%.%df',digits2use);
      et.onset    = strrep(cellstr(num2str(et.onset,myExpr)),' ','');
      et.duration = strrep(cellstr(num2str(et.duration,myExpr)),' ','');
      writetable(et,fn_events,'filetype','text','delim','tab');
    end
  end
  
  
    
end
  