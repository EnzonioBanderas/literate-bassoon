function hfig = ter_plotBidsPhysioTsv(fn_tsv)
  
  hfig = figure('units','pixel','outerposition',[0 0 1000 1414]);
  [data,t] = ter_readBidsTsv(fn_tsv);
  data.time = t;
  data.Properties.VariableUnits(size(data,2)) = {'s'};
  ind_numeric = false(1,size(data,2));
  for i=1:size(data,2)
    ind_numeric(i) = isnumeric(data.(data.Properties.VariableNames{i}));
  end
  data = data(:,ind_numeric);
 
  %pinfo = ter_parseFname(fn_tsv);
  %pID = pinfo.pID;
  tmp = ter_parseFname(fn_tsv);
  fn_events = strrep(fn_tsv,['_' tmp.suffix] ,'_events');
  try
    if strcmpi(fn_events(end-2:end),'.gz')
      fn_events = fn_events(1:end-3);
    end
  catch
  end
  if exist(fn_events,'file')==2
    et = readtable(fn_events,'filetype','text','delim','tab');
  else
    et = table('size',[0,3],...
      'VariableNames',{'onset','duration','trial_type'},...
      'VariableTypes',{'double','double','cell'});
  end
  
  tt = et.trial_type;
  for i=1:numel(tt)
    tmp = regexp(tt{i},'_trial-[0-9]*','match');
    if not(isempty(tmp))
      tt{i} = strrep(tt{i},tmp{end},'');
    end
  end
  events = unique(tt);
  
  if numel(events)==0
    heightEv = 0;
  elseif numel(events) < 6
    heightEv = 1;
  else
    heightEv = 2;
  end
  heightPlots = size(data,2)-1 + heightEv;
      
  time_min = min([data.time;et.onset]);
  time_max = max([data.time;(et.onset+et.duration)]);
      
  clf
  myXLimits = [time_min-0.05*(time_max-time_min), ...
    time_max+.05*(time_max-time_min)]./60;
  for k=1:(size(data,2)-1)
    subplot(heightPlots,1,k);
    myLabel = data.Properties.VariableNames{k};
    myUnit  = data.Properties.VariableUnits{k};
    plot(data.time/60,data.(myLabel))
    if isempty(myUnit)
      ylabel(myLabel,'interpreter','none');
    else
      ylabel(sprintf('%s\n[%s]',myLabel,myUnit),'interpreter','none')
    end
    
    xlabel('t [min]')
    xlim(myXLimits)
    if min(data.(myLabel))==0 && max(data.(myLabel))==1
      ylim([-.1,1.1]);
    end
  end
      
  if numel(events)>0
    subplot(heightPlots,1,(heightPlots+1-heightEv):heightPlots);
    
    for k=1:numel(events)
      starts = et.onset(ismember(tt,events{k}));
      stops  = et.duration(ismember(tt,events{k})) + starts;
      
      tmp_time = sortrows([starts;starts;stops;stops]);
      tmp_line = repmat([0;1;1;0],numel(starts),1);
      tmp_time = [time_min;tmp_time;time_max]; %#ok<AGROW>
      tmp_line = [0;tmp_line;0];%#ok<AGROW>
      tmp_line = tmp_line.*1/numel(events)*.9/2 + (k-1)/numel(events);
      plot(tmp_time./60,tmp_line);
      
      text(tmp_time(1)/60,tmp_line(1)+3/4/numel(events),events{k},...
        'FontWeight','normal','interpreter','none','FontSize',8)
      hold on
    end
    ylabel('events');
    xlabel('t [min]');
    ylim([-0.05,1.05]);
    set(gca,'YTick',[])
    xlim(myXLimits);
  end
  [~,fn] = fileparts(fn_tsv);
  [~,fn] = fileparts(fn);
  suptitle(strrep(fn,'_','\_'));
  
end

