function func_study7_skinconductanceBids2RUBeval(fl_physiobids, dir_eval)
  % output sampling rate fixed to 1 kHz
  % 2022-06-29: put out unfiltered data as well

  device_label = 'GSR100C';
  bids_label = cellfun(@essbids_parseLabel,fl_physiobids, 'uni', 0);
  bids_label_include = true(length(bids_label), 1);
  for i=1:length(bids_label)
      if ~ismember('acq', fieldnames(bids_label{i}))
            fprintf('%s', fullfile(bids_label{i}.fpath, bids_label{i}.fname))
            bids_label_include(i) = false;
            delete(fullfile(bids_label{i}.fpath, bids_label{i}.fname))
      end
  end
  bids_label = bids_label(bids_label_include);
  fl_physiobids = fl_physiobids(bids_label_include);
  bids_label = vertcat(bids_label{:});

  if strcmp(dir_eval(end),filesep) && numel(dir_eval)>1
    dir_eval = dir_eval(1:end-1);
  end
  dir_eval0 = [dir_eval '_unfiltered'];

  for i=1:numel(fl_physiobids)
    fn_phys = fl_physiobids{i};
    fn_json = essbids_findJsonSidecar(fn_phys);
    
    bf = bids_label(i);
    bf = rmfield(bf,'type');
    bf.suffix = 'EDA';
    bf.extension = '.mat';
    fn_eda  = essbids_buildFileName(bf,'addpath',dir_eval);
    fn_eda0 = essbids_buildFileName(bf,'addpath',dir_eval0); % for unfiltered data
    
    
    bf = bids_label(i);
    bf.suffix = 'events';
    bf.extension = '.tsv';
    fn_events = essbids_buildFileName(bf,'addfilepath',true);
    
    if exist(fn_phys,'file')~=2
      error('physio file does not exist: %s',fn_phys);
    elseif exist(fn_json,'file')~=2
      error('json file does not exist: %s',fn_json);
    elseif exist(fn_events,'file')~=2
      error('event file does not exist: %s',fn_events);
    end
    fp_eval = fileparts(fn_eda);
    if ~isfolder(fp_eval)
      mkdir(fp_eval);
    end
    if ~isfolder(fileparts(fn_eda0))
      mkdir(fileparts(fn_eda0));
    end


    %% filter data and write EDA file to output
    if exist(fn_eda,'file')~=2
      fprintf('Processing %s\n',fn_eda);
      
      % read events data file
      evtab = essbids_readTsv(fn_events);
      
      % read physio data file
      phystab = essbids_readTsv(fl_physiobids{i});
      jinfo = phystab.Properties.CustomProperties.JsonSidecar;
      
      % remove artifacts and apply lowpass
      skinconductance0 = phystab.skinconductance;
      phystab.skinconductance = ...
        func_edaInterpolatePulsedArtifacts(fl_physiobids{i});
        


      % downsample to 1 kHz
      if jinfo.SamplingFrequency ~= 1000
        % 2do allow for any sampling rate different to full mulitples of
        % 1000
        fprintf('resampling from %d kHz to 1 kHz\n',...
          jinfo.SamplingFrequency/1000);
        if ismember(jinfo.SamplingFrequency/1000,1:20)
          myFactor = jinfo.SamplingFrequency/1000;
          tmp1 = phystab.skinconductance;
          tmp1 = [tmp1;nan(ceil(numel(tmp1)/myFactor)*myFactor-numel(tmp1),1)]; %#ok<AGROW> 
          tmp2 = tmp1(1:myFactor:end);
          for j=2:myFactor
            tmp2(:,j) = tmp1(j:myFactor:end);
          end
          skinconductance = nanmean(tmp2,2); %#ok<NANMEAN> 
          tmp1 = skinconductance0;
          tmp1 = [tmp1;nan(ceil(numel(tmp1)/myFactor)*myFactor-numel(tmp1),1)]; %#ok<AGROW> 
          tmp2 = tmp1(1:myFactor:end);
          for j=2:myFactor
            tmp2(:,j) = tmp1(j:myFactor:end);
          end
          skinconductance0 = nanmean(tmp2,2); %#ok<NANMEAN> 
          


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
      
      data          = skinconductance;
      isi           = 1;
      units         = jinfo.skinconductance.Units;
      labels        = device_label; 
      start_sample  = 0;      
      isi_units     = 'ms';
      save(fn_eda,'data','isi','units','labels','start_sample',...
        'isi_units','-v7.3');
    
      data          = skinconductance0;
      save(fn_eda0,'data','isi','units','labels','start_sample',...
        'isi_units','-v7.3');
    
    
      %% prepare trial definition file
      fn_td = strrep(fn_eda, '_EDA.mat','_trialDefinition.mat');
      fn_td0 = strrep(fn_eda0, '_EDA.mat','_trialDefinition.mat');
      fn_et = strrep(fn_eda, '_EDA.mat','_eventTable.mat');
      if exist(fn_td,'file')==2 && exist(fn_et,'file')==2
        continue
      end
      
      % ignore certain events
      events2exclude = {'instructions'};
      evtab = evtab(not(ismember(lower(evtab.trial_type),events2exclude)),:);
      % also ignore all events way to long, i.e. 20 seconds or longer
      evtab = evtab(evtab.duration < 20,:);
      % for this experiment select only CSplus and CSminus events, get
      % further info from reinforce column
      evtab = evtab(ismember(evtab.trial_type,{'CSplus','CSminus'}),:);
      evtab.onset = evtab.onset - jinfo.StartTime;
      
      % seconds before and after the stimulus
      t_before   =  0;
      dur_trial  = 20;
         
      numTrials       = size(evtab,1);
      trialStart      = evtab.onset - t_before;
      trialStop       = trialStart  + dur_trial;
      stimOnset       = evtab.onset;
      recordingSystem = device_label;
      samplRate       = 1000;
      trialSep        = []; % leave empty, trial splitting performed by 
                            % trialStart/Sstop
  
      offsetEIR  =  1.0;  % offset first  interval response
      offsetTIR  =  6.5;  % offset third  interval response, after CS onset
      durTIR     =  5;
          
      % build trial seperators
      sPos{1}.time  = stimOnset+offsetEIR;
      sPos{1}.name  = 'EIR';
      sPos{1}.color = 'c';
      
      sPos{2}.time  = stimOnset+offsetTIR+durTIR;
      sPos{2}.name  = 'end';
      sPos{2}.color = 'c';
      
      % markers to identify CS and US types (sPos 3 to 6), set time to -1 
      % if not present. Use loop over events table to fill in proper valies
      sPos{3}.time  = -1*ones(size(evtab,1),1);
      sPos{3}.name  =  'CSplus';
      sPos{3}.color =  'b';
      
      sPos{4}.time  = -1*ones(size(evtab,1),1);
      sPos{4}.name  =  'CSminus';
      sPos{4}.color =  'b';
      
      sPos{5}.time  = -1*ones(size(evtab,1),1);
      sPos{5}.name  =  'noUS';
      sPos{5}.color =  'b';
      
      sPos{6}.time  = -1*ones(size(evtab,1),1);
      sPos{6}.name  =  'US';
      sPos{6}.color =  'b';
      
      for j=1:size(evtab,1)
        if ismember(evtab.trial_type(j),'CSplus')
          sPos{3}.time(j) = evtab.onset(j);
        elseif ismember(evtab.trial_type(j),'CSminus')
          sPos{4}.time(j) = evtab.onset(j);
        end
        if isequal(evtab.reinforced(j),1) ||...
            isequal(evtab.reinforced(j),{'1'}) 
          sPos{6}.time(j) = evtab.onset(j)+5.9;
        else
          sPos{5}.time(j) = evtab.onset(j)+5.9;
        end
      end
      
      % save data
      save(fn_td,'numTrials','recordingSystem','sPos','samplRate',...
        'stimOnset','trialSep','trialStart','trialStop')
      copyfile(fn_td,fn_td0);
    end

  end

end


