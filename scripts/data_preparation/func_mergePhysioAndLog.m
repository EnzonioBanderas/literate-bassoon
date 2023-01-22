function func_mergePhysioAndLog(fp_d,fp_s)


  fprintf('\n%s\n%s : Merging log and physio data\n',...
    repmat('=',1,72),datestr(now,'yyyy-mm-dd hh:MM:ss'));
  
  tolerance_ev = .02;%05; % sx
  rf = 3;  % rounding value
  code_scantrig = '1';
  %type_scantrig = 'Pulse';
  phaseIDtext   = 'Experiment part : ';
  trialtypes2ignore = {'Shock' 'US'};% 'Trial' 
  spacingWithinUS = 33; % ms
  id_CS = 'CS';
  id_US = {'ShockOut'};
  id_instructions = {'InstructionText'};
   
  fn_ptab = fullfile(fp_d,'participants.tsv');
  tsv_param = {'filetype','text','delim','tab','treat','n/a'};
  %ptab = readtable(fn_ptab,tsv_param{:});
  ptab = ter_readptab(fn_ptab);
  
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    fl_events = ter_listFiles(fullfile(fp_d,pID),[pID '*_events.tsv']);
    [~,sess] = cellfun(@(x) fileparts(fileparts(fileparts(x))),...
      fl_events,'uni',0);
    
    
    uSess = unique(sess);
    logs = cell(size(sess));
    skip_subject = false;
    for j=1:numel(uSess)
      clearvars logdata
      sID = uSess{j};
      fn_log = fullfile(fp_s,pID,sID,'stimlogs',...
        strrep([pID '_' sID '_logdataRuns.mat'],'__','_'));
      if exist(fn_log,'file')==2
        load(fn_log,'logdata');
        ind = find(ismember(sess,sID));
        myFields = fieldnames(logdata);
        if numel(ind) ~= numel(myFields)
          continue
        end
        for k=1:numel(myFields)
          logs{ind(k)} = logdata.(myFields{k});
        end
      else
        warning('Missing log data: %s',fn_log)
        skip_subject = true;
      end
    end
    if skip_subject
      fprintf('incomplete log data for this participant %s, skipping\n',...
        pID);
      continue
    end
    
    for j=1:numel(sess)
      log = logs{j};
      if isempty(log)
        continue
      end
      %isok = true;
      %read in the session event table 
      et = readtable(fl_events{j},'FileType','text','delim','tab');
      % check for unidentified CS, skip process if none present, i.e. already
      % identified
      sID = sess{j};
      rID = regexp(fl_events{j},'run-[0-9]*','match');
      %disp(fl_events{j})
      try
        rID = rID{end};
      catch
        warning('run inidentified for: %s',fl_events{j})
      end
      if ~ismember(id_CS,et.trial_type) % 
        %fprintf('Log and phys events already synced, skipping: %s, %s\n',...
        %    pID,sess{j});
        continue
      else
        fprintf('Syncing log and phys events for : %s, %s, %s\n',...
            pID,sID,rID);   
      end
      et(ismember(et.trial_type,trialtypes2ignore),:)=[];
      
      
    
      %% identify the CS und US events in the event table
      ind_us = ismember(et.trial_type,id_US);
      indUS = find(ind_us);
      durUS = repmat(0.001,size(indUS));
      rmUS  = false(size(durUS));
      for l = numel(indUS):-1:2
        if et.onset(indUS(l)) - et.onset(indUS(l-1)) < ...
            spacingWithinUS*1e-3 + tolerance_ev
          rmUS(l) = true; 
          durUS(l-1) = durUS(l)+et.onset(indUS(l))-et.onset(indUS(l-1));
        end
      end
      et.duration(indUS(~rmUS)) = durUS(~rmUS);
      et.trial_type(indUS(~rmUS)) = {'US'};
      et(indUS(rmUS),:) = [];
      ind_cs = ismember(et.trial_type,id_CS); % logical vector
      %ind_cs2 = find(ind_cs);                % numbered vector 
      ind_us = ismember(et.trial_type,'US');
      %indUS = find(ind_us);
      
      
      ind_stim = find(ind_cs |ind_us);
      % calculate time differences between the events of the event table
      %timediff_phys_cs = diff(et.onset(ind_cs)) ; 
      %timediff_phys_us = diff(et.onset(ind_us)) ;
      timediff_phys_stim = diff(et.onset(ind_stim)) ; %#ok<FNDSB>

      %% now the logs
      %for k=1:numel(log)
        % select the log info for each run, 
      try
        el = log.eventList;
      catch
        continue
      end
      % remove scantriggers to make reading easier
      el = el(not(ismember(el.Code,code_scantrig)),:); 
      
      el.Time = el.Time / 1e4; % transfer time to seconds
      
      el = el(not(contains(el.Code,'CS+ : ')),:);
      el = el(not(contains(el.Code,'Context : ')),:);
      
      % handle pulse trains of US, adding up to one events, calculate
      % duration, 
      indUS = find(ismember(el.Code,'electricalUS'));
      durUS = repmat(0.001,size(indUS));
      rmUS  = false(size(indUS));
      for l = numel(indUS):-1:2
        if el.Time(indUS(l)) - el.Time(indUS(l-1)) < .04
          rmUS(l) = true; 
          durUS(l-1) = durUS(l)+el.Time(indUS(l))-el.Time(indUS(l-1));
        end
      end
      durUS = durUS(not(rmUS));
      el(indUS(rmUS),:) = [];
      indUS = find(ismember(el.Code,'electricalUS'));
      
      % identify phase 
      indPhaseID  = find(contains(el.Code,phaseIDtext),1,'last');
      phaseID     = strrep(strrep(el.Code{indPhaseID},phaseIDtext,''),...
        ' ','');
      
      % find CS events 
      indCSons    = not(cellfun(@isempty,regexp(el.Code,'^CS[\+\-]')));
      %ismember(el.Code,{'CS+','CS-'});
      indCSend    = false(size(indCSons));
      indCSend_   = find(ismember(el.Code,'fixationImage'));
      indCSons_   = find(indCSons);
      for l=1:numel(indCSons_)
        indCSend(indCSend_(find(indCSend_>indCSons_(l),1,'first'))) = true;
      end
      dur_cs      = el.Time(indCSend) - el.Time(indCSons);

      %timediff_logs_cs = diff(el.Time(indCSons));
        
      % prepare relevant labels for CS events
      stims       = cellfun(@(x) ['stim-' x], ...
        el.Code([indCSons(2:end);false]),'uni',0);
      phases      = repmat({['phase-' phaseID]},size(stims));
%       contexts    = strrep(strrep(...
%         el.Code([indCSons(3:end);false;false]),' ','-'),'Cont','cont');
%       cs_labels   = cellfun(@(x,y,z,a) [x '_' y '_' z '_' a],...
%         strrep(strrep( el.Code(indCSons),'-','minus'),'+','plus'),...
%         phases,contexts,stims,'uni',0);
      cs_labels   = cellfun(@(x,y,z,a) [x '_' y '_' z],...
        strrep(strrep( el.Code(indCSons),'CS-','CSminus'),...
        'CS+','CSplus'),phases,stims,'uni',0);
  
      % prepare relevant labels for US events
      %contexts = cell(size(indUS));
      prevStim = cell(size(indUS));
      for l=1:numel(indUS)
        %look at logs up to 40 sec before US
        ind_40sec = el.Time < el.Time(indUS(l)) & ...
          el.Time(indUS(l))-el.Time < 40;
        ev_before = el.Code(ind_40sec);
%         ind_context = find(contains(ev_before,'Context '),1,'last');
%         if isempty(ind_context)
%           contexts{l} = 'context-none';
%         else
%           context = ev_before{ind_context};
%           context = ['context-' context(end)];
%           contexts{l} = context;
%         end
        %ind_stimUS = find(ismember(ev_before,{'CS+','CS-'}),...
        ind_stimUS = find(not(...
          cellfun(@isempty,regexp(ev_before,'^CS[\+\-]'))),1,'last');
        if isempty(ind_stimUS)
          prevStim{l} = 'prevstim-none';
        else
          tmp = strrep(ev_before{ind_stimUS}(1:3),'+','plus');
          tmp = strrep(tmp,'-','minus');
          prevStim{l} = ['prevstim-' tmp];
        end
      end
      us_labels = cellfun(@(x,y) ['US_phase-' phaseID '_' x],...
        prevStim,'uni',0);

      % create a table of US and CS stimuli in log file
      clearvars stimtab
      stimtab.indStim = [find(indCSons);indUS];
      stimtab.durStim = [dur_cs;durUS];
      stimtab.labels  = [cs_labels;us_labels];
      stimtab.time    = el.Time(stimtab.indStim);
      stimtab = sortrows(struct2table(stimtab));

      timediff_logs_stim = diff(stimtab.time);

      if numel(timediff_logs_stim) ~= numel(timediff_phys_stim)
        warning('different number of events in log and event table');
        continue
      end

      td_stim = abs(round(timediff_logs_stim - ...
         timediff_phys_stim,3));
      if max(td_stim) > tolerance_ev 
        warning('misfit: %s, %s, %d, max diff %d ms',pID,sess{j},k,...
         round(max(td_stim*1000)))
        isok=false;
      else
        ind_CSUS = ismember(et.trial_type,{'CS','US'});
        et.trial_type(ind_CSUS) = stimtab.labels;
        et.duration(ind_CSUS)   = stimtab.durStim;
        isok = true;
        %logsUsed = [logsUsed;fn_logs{k}]; %#ok<AGROW>
        %break
      end
      
      instr = el(find(ismember(el.Code,id_instructions),1,'last'),:);
      if not(isempty(instr))
        instr.Duration = instr.Duration ./1e4;
        try
          td = mean(et.onset(ismember(et.trial_type,stimtab.labels),:)-...
            stimtab.time);
          et.onset(end+1) = round(instr.Time+td,3);
          et.duration(end) = instr.Duration;
          et.trial_type{end} = 'Instructions';
        catch
          warning('trouble with the instructions ');
        end
      end
      
      % if information is sufficient write the updated event.tsv file 
      if isok
        temp_evtable = sortrows(et); % sort in the added context info
        temp_evtable.duration = arrayfun(@(x) eval(...
          sprintf('sprintf(''%%.%df'',%d)',rf,x)),...
          temp_evtable.duration,'uni',0);
        temp_evtable.onset = arrayfun(@(x) eval(...
          sprintf('sprintf(''%%.%df'',%d)',rf,x)),...
          temp_evtable.onset,'uni',0);
        writetable(temp_evtable,fl_events{j},'filetype','text',...
          'delim','tab');
      end
    end
  end
  
  
  %%update ptab
  ptab0 = ptab;
  if ~ismember('cs_plus',ptab.Properties.VariableNames)
    ptab.cs_plus = repmat({'n/a'},size(ptab,1),1);
  end
  for i=1:size(ptab,1)
    if ismember(ptab.cs_plus(i),'n/a')
      pID = ptab.participant_id{i};
      fl = ter_listFiles(fullfile(fp_d,pID),[pID '*events.tsv'],4);
      if isempty(fl)
        continue
      end
      tt = {};
      for j=1:numel(fl)
        t1 = readtable(fl{j},'filetype','text','delim','tab');
        tt = unique([tt;t1.trial_type]);
      end
      tt = tt(contains(tt,'CSplus_'));
      picCSplus = unique(cellfun(@(x) regexp(x,'-Pic[a-zA-Z]*','match'),tt));
      if numel(picCSplus)>1
        warning('CSplus image is not consistent: %s',pID)
      else
        if isempty(picCSplus)
          fprintf('No CS plus detected for subject %s\n', pID);
        else
          ptab.cs_plus{i} =strrep(lower(picCSplus{1}),'-pic','');
        end
      end
    end
  end
  ptab = sortrows(ptab);
  if not(isequaln(ptab,ptab0))
    writetable(ptab,fn_ptab,tsv_param{1:4});
    txt  = fileread(fn_ptab);
    txt2 = strrep(txt, 'NaN','n/a');
    txt2 = strrep(txt2,'NaT','n/a');
    if not(isequal(txt,txt2))
      fid = fopen(fn_ptab,'w');
      fprintf(fid,txt2);
      fclose(fid);
    end
  end




end