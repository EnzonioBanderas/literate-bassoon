function ter_skinconductanceBids2RUBeval(fl_physiobids, dir_eval)
  % output sampling rate fixed to 1 kHz
  % 2021-11-04 : 

  pIDs = cell(size(fl_physiobids));
  sess = cell(size(fl_physiobids));
  runs = cell(size(fl_physiobids));
  for i=1:numel(fl_physiobids)
    [~,fn,~] = fileparts(fl_physiobids{i});
    pID = regexp(fn,'sub-[a-zA-Z0-9]*','match');
    pIDs{i} = pID{1};
    ses = regexp(fn,'ses-[a-zA-Z0-9]*','match');
    if not(isempty(ses))
      sess{i} = ses{1};
    else
      sess{i} = '';
    end
    run = regexp(fn,'run-[0-9]*','match');
    if not(isempty(run))
      runs{i} = run{1};
    else
      runs{i} = '';
    end
  end
  
  for i=1:numel(fl_physiobids)
    pID = pIDs{i};
    ses = sess{i};
    
    fp_eval = fullfile(dir_eval,pID,ses);
    if ~isfolder(fp_eval)
      mkdir(fp_eval)
    end
    [~,fn,fe] = fileparts(fl_physiobids{i});
    if strcmpi(fe,'.gz')
      fn_eda    = fullfile(fp_eval,strrep(fn,'physio.tsv','EDA.mat'));
      fn_json   = strrep(fl_physiobids{i},'.tsv.gz','.json');
      fn_events = strrep(fl_physiobids{i},'physio.tsv.gz','events.tsv');
    else
      fn_eda    = fullfile(fp_eval,strrep(fn,'physio','EDA.mat'));
      fn_json   = strrep(fl_physiobids{i},'.tsv','.json');
      fn_events = strrep(fl_physiobids{i},'physio.tsv','events.tsv');
    end
    if exist(fn_json,'file')~=2
      error('json file does not exist: %s',fn_json);
    elseif exist(fn_events,'file')~=2
      error('event file does not exist: %s',fn_events);
    end
    
    evtab = readtable(fn_events,'filetype','text','delim','tab');
    jinfo = jsondecode(fileread(fn_json));
    
    
    %% filter data and write EDA file to output
    if exist(fn_eda,'file')~=2
      %% read physio data file
      %phystab = ter_readBidsTsv(fl_physiobids{i});
      phystab = essbids_readTsv(fl_physiobids{i});
      phystab.skinconductance = func_rmEdaArtefacts(fl_physiobids{i});
      %pause
%       if strcmpi(fe,'.gz')
%         gunzip(fl_physiobids{i},fp_eval)
%         phystab = readtable(fullfile(fp_eval,fn),'filetype','text',...
%           'delim','tab');
%         delete(fullfile(fp_eval,fn));
%       else
%         phystab = readtable(fl_physiobids{i},'filetype','text',...
%           'delim','tab');
%       end
%       phystab.Properties.VariableNames = jinfo.Columns;
      %% prepare EDA data file
      if jinfo.SamplingFrequency ~= 1000
        %error('sampling rate has to be 1000, transform for other sampling rates not yet implemented')
        %2do allow for any sampling rate
        fprintf('resampling from %d kHz to 1 kHz\n',...
          jinfo.SamplingFrequency/1000);
        if ismember(jinfo.SamplingFrequency/1000,1:20)
          myFactor = jinfo.SamplingFrequency/1000;
          tmp1 = phystab.skinconductance;
          tmp1 = [tmp1;nan(ceil(numel(tmp1)/myFactor)*myFactor-numel(tmp1),1)];
          tmp2 = tmp1(1:myFactor:end);
          for j=2:myFactor
            tmp2(:,j) = tmp1(j:myFactor:end);
          end
          skinconductance = nanmean(tmp2,2);
        else
          error('');
          %skinconductance = resample(phystab.skinconductance,...
          %  1000,jinfo.SamplingFrequency);
        end
      else
        skinconductance = phystab.skinconductance;
      end
      
      ind_nan = isnan(skinconductance);
      
      if sum(ind_nan)>0
        skinconductance = ter_interpolateNan(skinconductance);
      end
      
      %fprintf('%s - applying bandpass to : %s\n',...
      %  datestr(now,'yyyy-mm-dd HH:MM:SS'),fn);
      %data_filtered = ter_bandpassAnalogBiopac(skinconductance,1000);
      %data          = data_filtered;
      data = skinconductance;
      %data(ind_nan) = nan;
      isi           = 1;
      units         = jinfo.skinconductance.Units;
      labels        = 'GSR100C'; 
      start_sample  = 0;      
      isi_units     = 'ms';
      save(fn_eda,'data','isi','units','labels','start_sample',...
        'isi_units','-v7.3');
    end
    
    
    %% prepare trial definition file
    fn_td = strrep(fn_eda, '_EDA.mat','_trialDefinition.mat');
    fn_et = strrep(fn_eda, '_EDA.mat','_eventTable.mat');
    if exist(fn_td,'file')==2 && exist(fn_et,'file')==2
      continue
    end
    
    % ignore certain events
    events2exclude = {'instructions'};
    evtab = evtab(not(ismember(lower(evtab.trial_type),events2exclude)),:);
    % also ignore all events way to long, i.e. 20 seconds or longer
    evtab = evtab(evtab.duration < 20,:);
    evtab.onset = evtab.onset - jinfo.StartTime;
    
    % seconds before and after the stimulus
    t_before   =  4;
    t_after    =  8;
    offsetFIR  =  1.0;  % offset first  interval response
    offsetSIR  =  1.0;  % offset second interval response
    offsetTIR  =  0.5;  % offset third  interval response
    durTIR     =  5;
    threshold_FIR_SIR = 7.5;
    threshold_longEvents = 8;
     
    
    %#ok<*AGROW>
    event_label     = {};
    trial_index     = [];
    edaeval_index   = 0;
    evalBlock_start = [];
    evalBlock_stop  = [];
    event_onsets     = [];
    interval_boundaries_lower = [];
    interval_boundaries_upper = [];
    count_trialindex = 0;
    evtab = evtab(ismember(evtab.trial_type,{'CSplus','CSminus'}),:);
    numTrials = size(evtab,1);
    trialStart = evtab.onset-2;
    trialStop = trialStart + 18;
    stimOnset = evtab.onset;
    recordingSystem = 'GSR100C';
    samplRate = 1000;
    sPos{1}.time = evtab.onset+1;
    sPos{1}.name =  'EIR';
    sPos{1}.color =  'c';
    sPos{2}.time = evtab.onset+6+5+.5;
    sPos{2}.name =  'end';
    sPos{2}.color =  'c';
    
    sPos{3}.time  = -1*ones(size(evtab,1),1);
    sPos{3}.name  =  'CSplus';
    sPos{3}.color =  'b';
    
    sPos{4}.time  = -1*ones(size(evtab,1),1);
    sPos{4}.name  =  'CSminus';
    sPos{4}.color =  'b';
    
    sPos{5}.time  = -1*ones(size(evtab,1),1);
    sPos{5}.name  =  'TIR-noUS';
    sPos{5}.color =  'b';
    
    sPos{6}.time  = -1*ones(size(evtab,1),1);
    sPos{6}.name  =  'TIR-US';
    sPos{6}.color =  'b';
    
    for j=1:size(evtab,1)
      if ismember(evtab.trial_type(j),'CSplus')
        sPos{3}.time(j) = evtab.onset(j);
      elseif ismember(evtab.trial_type(j),'CSminus')
        sPos{4}.time(j) = evtab.onset(j);
      end
      if isequal(evtab.reinforced(j),1) || isequal(evtab.reinforced(j),{'1'}) 
        sPos{6}.time(j) = evtab.onset(j)+6.4;
      else
        sPos{5}.time(j) = evtab.onset(j)+6.4;
      end
    end
    trialSep=[];
    save(fn_td,'numTrials','recordingSystem','sPos','samplRate',...
      'stimOnset','trialSep','trialStart','trialStop')
%       end
%         
%       label0 = evtab.trial_type{j};
%       if 1<evtab.duration(j) && evtab.duration(j)< threshold_FIR_SIR 
%         % treat prolonged, but still short single events as FIR and TIR only
%         fulldur = round(evtab.duration(j),3);
%         count_trialindex = count_trialindex+1;
%         %FIR
%         event_onsets(end+1)    = evtab.onset(j) ;
%         evalBlock_start(end+1) = evtab.onset(j) - t_before;
%         evalBlock_stop(end+1)  = evtab.onset(j) + t_after;
%         interval_boundaries_lower(end+1) = evtab.onset(j) + offsetFIR;
%         interval_boundaries_upper(end+1) = evtab.onset(j) + offsetTIR + ...
%           +  fulldur;
%         trial_index(end+1)      = count_trialindex;
%         edaeval_index(end+1)    = edaeval_index(end) + 1;
%         event_label{end+1}      = ['EIR_' label0];
%         %TIR
%         event_onsets(end+1)    = evtab.onset(j) + fulldur;
%         evalBlock_start(end+1) = evtab.onset(j) + fulldur - t_before;
%         evalBlock_stop(end+1)  = evtab.onset(j) + fulldur + t_after;
%         interval_boundaries_lower(end+1) = evtab.onset(j) + fulldur + ...
%           offsetTIR;
%         interval_boundaries_upper(end+1) = evtab.onset(j) + fulldur + ...
%           offsetTIR + durTIR;
%         trial_index(end+1)      = count_trialindex;
%         edaeval_index(end+1)    = edaeval_index(end) + 1;
%         event_label{end+1}      = ['TIR_' label0];
%         
%       elseif threshold_FIR_SIR < evtab.duration(j) && ...
%           evtab.duration(j)< threshold_longEvents
%         fulldur = round(evtab.duration(j),3);
%         halfdur = round(evtab.duration(j)/2,3);
%         count_trialindex = count_trialindex+1;
%         % FIR
%         event_onsets(end+1)     = evtab.onset(j) ;
%         evalBlock_start(end+1) = evtab.onset(j) - t_before;
%         evalBlock_stop(end+1)  = evtab.onset(j) + t_after;
%         interval_boundaries_lower(end+1) = evtab.onset(j) + offsetFIR;
%         interval_boundaries_upper(end+1) = evtab.onset(j) + offsetSIR + ...
%           +  halfdur;
%         trial_index(end+1)      = count_trialindex;
%         edaeval_index(end+1)   = edaeval_index(end) + 1;
%         event_label{end+1}       = ['FIR_' label0];
%         % SIR
%         event_onsets(end+1)    = evtab.onset(j) + halfdur;
%         evalBlock_start(end+1) = evtab.onset(j) + halfdur - t_before;
%         evalBlock_stop(end+1)  = evtab.onset(j) + halfdur + t_after;
%         interval_boundaries_lower(end+1) = evtab.onset(j) + offsetSIR + ...
%           + halfdur;
%         interval_boundaries_upper(end+1) = evtab.onset(j) + offsetTIR + ...
%           + fulldur;
%         trial_index(end+1)     = count_trialindex ;
%         edaeval_index(end+1)   = edaeval_index(end) + 1;
%         event_label{end+1}     = ['SIR_' label0];
%         % TIR
%         event_onsets(end+1)    = evtab.onset(j) + fulldur;
%         evalBlock_start(end+1) = evtab.onset(j) + fulldur - t_before;
%         evalBlock_stop(end+1)  = evtab.onset(j) + fulldur + t_after;
%         interval_boundaries_lower(end+1) = evtab.onset(j) + fulldur + ...
%           offsetTIR;
%         interval_boundaries_upper(end+1) = evtab.onset(j) + fulldur + ...
%           offsetTIR + durTIR;
%         trial_index(end+1)     = count_trialindex ;
%         edaeval_index(end+1)   = edaeval_index(end) + 1;
%         % next few lines are specific for the frea conditioning case here,
%         % markup TIR depending on apearance or omission of US
%         if j == size(evtab,1)
%           event_label{end+1} = ['TIR-noUS_' label0];
%         elseif strcmpi(evtab.trial_type{j+1}(1:2),'US')
%           event_label{end+1} = ['TIR-US_' label0];   
%         else
%           event_label{end+1} = ['TIR-noUS_' label0]; 
%         end
%       else
%         % treating single events as TIR equivalent, as they usually are
%         % shocks
%         addTIR = false;
%         if j==1 
%           addTIR = true;
%         elseif diff(evtab.onset([j-1,j]))+diff(evtab.duration([j-1,j]))>.01 
%           addTIR = true;
%         end
%         if addTIR
%           count_trialindex = count_trialindex+1;
%           event_onsets(end+1)    = evtab.onset(j);
%           evalBlock_start(end+1) = evtab.onset(j) - t_before;
%           evalBlock_stop(end+1)  = evtab.onset(j) + t_after;
%           interval_boundaries_lower(end+1) = evtab.onset(j) + offsetTIR;
%           interval_boundaries_upper(end+1) = evtab.onset(j) + offsetTIR ...
%             + durTIR;
%           trial_index(end+1)     = count_trialindex;
%           edaeval_index(end+1)   = edaeval_index(end) + 1;
%           event_label{end+1}     = ['TIR_' label0];
%         end
%       end
%     end
%     event_onsets              = event_onsets';
%     trial_index               = trial_index';
%     edaeval_index             = edaeval_index(2:end)';
%     event_label               = event_label';
%     evalBlock_start           = evalBlock_start';
%     evalBlock_stop            = evalBlock_stop'; 
%     interval_boundaries_lower = interval_boundaries_lower';
%     interval_boundaries_upper = interval_boundaries_upper';
%     eventtable = table(event_onsets,trial_index,edaeval_index,...
%       event_label,evalBlock_start,evalBlock_stop,...
%       interval_boundaries_lower,interval_boundaries_upper);
% 
%     
%     %% now write the trialDefinition file
%     s = sprintf([...
%       '%% automatically created trial definition file for EDA eval\n'...
%       'numTrials       = %d;\n'...
%       'trialSep        = [];\n'...
%       'trialStart      = [%s];\n'...
%       'trialStop       = [%s];\n'...
%       'stimOnset       = [%s];\n'],...
%       numel(evalBlock_start),...
%       ter_formatedOutput(evalBlock_start),...
%       ter_formatedOutput(evalBlock_stop ),...
%       ter_formatedOutput(event_onsets    ));
%     
%     
%     %% markers for the lower interval boundaries
%     countMarkers = 0;
%     if numel(interval_boundaries_lower)>0
%       countMarkers = countMarkers + 1;
%       markerlist = interval_boundaries_lower';
%       for k = (numel(markerlist)-1):-1:1
%         if markerlist(k) < markerlist(k+1) && markerlist(k) == 0
%           markerlist(k) = markerlist(k+1);
%         end
%       end
%       s = sprintf(['%s\n'...
%         'sPos{%d}.time    = [%s];\n'...
%         'sPos{%d}.color   = ''c'';\n'...
%         'sPos{%d}.name    = ''intervalBoundary_l_o_w'';\n'],...  
%       s, ...
%       countMarkers,...
%       ter_formatedOutput(markerlist),...
%       countMarkers,...
%       countMarkers);
%     end
%     
%     
%     %% markers for the upper boundaries
%     if numel(interval_boundaries_upper)>0
%       countMarkers = countMarkers + 1;
%       markerlist = interval_boundaries_upper';
%       for k=2:numel(markerlist)
%         if markerlist(k) == 0
%           markerlist(k) = markerlist(k-1);
%         end
%       end
%       s = sprintf(['%s\n'...
%         'sPos{%d}.time    = [%s];\n'...
%         'sPos{%d}.color   = ''g'';\n'...
%         'sPos{%d}.name    = ''intervalBoundary_u_p'';\n'],...
%         s, ...
%         countMarkers,...
%         ter_formatedOutput(markerlist),...
%         countMarkers,...
%         countMarkers);
%     end
%     
%     
%     %% finalize the trialDefinition file
%     s = sprintf(['%s\n'...
%       'samplRate       = %d;\n'...
%       'recordingSystem = ''EDA100C-MRI'';\n\n\n'],...
%       s,1000); % now sampling rate fixed to 1 kHz
%       %jinfo.SamplingFrequency);
%     for k=1:numel(event_label)
%       s = sprintf('%s%% %2.0f. %s\n',s,k,event_label{k});
%     end
%     fid      = fopen(fn_td,'w');
%     fprintf(fid,'%s',s);
%     fclose(fid);
%     % save event table
%     save(fn_et,'eventtable','-v7.3');
%     
    
  end
end



function mystring = ter_formatedOutput(numArray)
  s = ''; 
  for i=1:numel(numArray)
    s = sprintf('%s %9.4f',s,numArray(i));
    if mod(i,5) == 0 && i ~= numel(numArray)
      s = sprintf('%s ...\n  ',s);
    end
  end
  mystring = s;
end