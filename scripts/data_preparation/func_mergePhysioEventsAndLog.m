function [bSuccess,fn_events,time0] = func_mergePhysioEventsAndLog(fn_physioevents,fn_log)
% function merges temporary physio events table with NeuroBS presentation
% log file

  label_TR_phys = 'trigger';
  label_TR_log  = 'Pulse';
  
  log_codes2ignore = {'65','97'};
  
  id_CS = 'visualCS';
  id_US = {'electricStimOut'};
  
  tolerance_ev = .02;%05; % sx
  rf = 3;  % rounding value
  spacingWithinUS = 33; % ms
  
  
  log = essbids_readNeuroBSPresentationStimLog(fn_log,'tab');
  et = essbids_readTsv(fn_physioevents);
  

  if contains(fn_physioevents,{'7179','7180','7161','7162','7172','7173','7175'})
    et.trial_type(ismember(et.trial_type,'electricStim')) = {'electricStimOut'};
    et.trial_type(ismember(et.trial_type,'visualCS')) = {'electricStim'};
    et.trial_type(ismember(et.trial_type,'trial')) = {'visualCS'};
  end  
    
  %TR_et = median(round(diff(...
  %  et.onset(ismember(et.trial_type,label_TR_phys))),2));
  
  
  % identify timepoint T=0;
  %ind_trig = find(ismember(el.EventType,label_TR_log));
  %bSuccess  = false;
  fn_events = '';
  time0     = nan;
  
  

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
  ind_us = ismember(et.trial_type,'US');
  
  %% calculate time differences between the events of the event table
  ind_stim = find(ind_cs |ind_us);
  %timediff_phys_cs = diff(et.onset(ind_cs)) ; 
  %timediff_phys_us = diff(et.onset(ind_us)) ;
  timediff_phys_stim = diff(et.onset(ind_stim)) ; %#ok<FNDSB>
  
  %% identify the CS und US events in the presentation log table
  el = log.GeneralEventList;
  el = el(not(ismember(string(el.Code),log_codes2ignore)),:);
  el.Time = el.Time/1e4; 
  el = el(not(contains(string(el.Code),'CS+ : ')),:);
  el = el(not(contains(string(el.Code),'Context : ')),:);
  el = el(not(ismember(string(el.EventType),label_TR_log)),:);
  
  % handle pulse trains of US, adding up to one event, calculate duration
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

  % find CS events 
  indCSons    = not(cellfun(@isempty,regexp(string(el.Code),'^CS[\+\-]')));
  %ismember(el.Code,{'CS+','CS-'});
  indCSend    = false(size(indCSons));
  indCSend_   = find(ismember(string(el.Code),'fixationImage'));
  indCSons_   = find(indCSons);
  for l=1:numel(indCSons_)
    indCSend(indCSend_(find(indCSend_>indCSons_(l),1,'first'))) = true;
  end
  dur_cs      = el.Time(indCSend) - el.Time(indCSons);
  
  labels_cs=string(el.Code(indCSons));
  labels_us=repmat({'US'},numel(indUS),1);
  
  

  % create a table of US and CS stimuli in log file
  clearvars stimtab
  stimtab.indStim = [find(indCSons);indUS];
  stimtab.durStim = [dur_cs;durUS];
  stimtab.time    = el.Time(stimtab.indStim);
  stimtab.labels  = [labels_cs;labels_us];
  stimtab = sortrows(struct2table(stimtab));

  timediff_logs_stim = diff(stimtab.time);
  
  if numel(timediff_logs_stim) == numel(timediff_phys_stim)
    td_stim = abs(round(timediff_logs_stim - timediff_phys_stim,4));
    bSuccess = max(td_stim)<tolerance_ev;
  else
    minTrig = min(numel(timediff_logs_stim),numel(timediff_phys_stim));
    if minTrig>4
      fprintf(['Number of trigger pulses mismatch, '...
        'comparing only first %d triggers:\n'...
        '   phys    : %d\n   stimlog : %d\n'],...
        minTrig,numel(timediff_phys_stim),numel(timediff_logs_stim));
      td_stim = abs(round(timediff_logs_stim(1:minTrig) - ...
        timediff_phys_stim(1:minTrig),4));
      bSuccess = max(td_stim)<tolerance_ev;
    else
      bSuccess=false;
    end
  end
  if not(bSuccess)
    warning('sync:misfit','Logs do not fit:\n  %s\nand\n  %s',...
      fn_log,fn_physioevents)
    return
  end
  
  relativeTime = min(stimtab.time)-min(et.onset(ind_stim));
  stimtab.time = stimtab.time -relativeTime;
  
 
  
  el.Time = el.Time - relativeTime;
  
  eventtab = table();
  eventtab.onset      = round(stimtab.time,rf);
  eventtab.duration   = round(stimtab.durStim,rf);
  eventtab.trial_type  = cell(size(eventtab,1),1);
  eventtab.trial_index = nan(size(eventtab,1),1);
  eventtab.reinforced  = nan(size(eventtab,1),1);
  eventtab.stim_file  = cell(size(eventtab,1),1);
   
  for i=1:size(stimtab)
    label = stimtab.labels{i};
    tmp = strsplit(label,'_');
    if strcmpi(tmp(1),'CS+')
      eventtab.trial_type(i) = {'CSplus'};
    elseif strcmpi(tmp(1),'CS-')
      eventtab.trial_type(i) = {'CSminus'};
    else
      eventtab.trial_type(i) = tmp(1);
    end
    tmp = regexp(label,'trial-(?<ti>[0-9]+)','names');
    if not(isempty(tmp))
      eventtab.trial_index(i) = str2double(tmp(1).ti);
    elseif i>1
      eventtab.trial_index(i) = eventtab.trial_index(i-1);
    end
    if ismember(eventtab.trial_type(i),{'CSplus','CSminus'})
      eventtab.reinforced(i) = 1.0* contains(label,'Paired');
    end
    if ismember(eventtab.trial_type(i),{'CSplus','CSminus'})
      ind_pic = find(ismember(el.Trial,el.Trial(stimtab.indStim(i))) & ...
        ismember(el.EventType,'Picture'));
      eventtab.stim_file{i}  = strcat('images/',...
        strrep(lower(char(el.Code(ind_pic(1)))),'pic',''),'.png');
    else
      eventtab.stim_file{i} = 'n/a';
    end
  end
  disp(eventtab);
  
  fn_events = strrep(GetFullPath(fn_physioevents),'_physioevents','_events');
  %delete(fn_physioevents);
  essbids_writeTsv(fn_events,eventtab);
  
  time0 = relativeTime;
end