function ter_plotBIDSPhysData(fp_d)
 
  ptab = readtable(fullfile(fp_d,'participants.tsv'),'filetype','text',...
    'delim','tab');
    
  % test wether export_fig is available
  tmp = which('export_fig');
  if isempty(tmp)
    error('Cannot find export_fig toolbox. Please add toolbox to path')
  end
  
%   try
%     spm('version');
%   catch
%     startspm 12;
%     spm('quit');
%   end
%   hFig = spm_figure;
%   warning('off','MATLAB:legend:IgnoringExtraEntries');

%   fl_p1 = {};
%   fl_p2 = {};
%   fn_out1 = fullfile(fileparts(fp_data),'plotPhys_all.pdf');
%   fn_out2 = fullfile(fileparts(fp_data),'plotPhys_byPhase.pdf');
%   
%   fl = cellstr(spm_select('FPListRec',fp_data,'plotPhys.*\.pdf'));
%   [~,pdir,~] = cellfun(@(x) fileparts(fileparts(x)),fl,'uni',0);
%   fl = fl(ismember(pdir,'behavior'));
%   fl1_old =fl(contains(fl, '_all.pdf'));
%   fl2_old =fl(contains(fl,'_byPhase.pdf'));
  fn_out = fullfile(fileparts(fp_d),'physio_plots.pdf');
  if exist(fn_out,'file')==2
    count = 1;
    while exist(strrep(fn_out,'.pdf',sprintf('_%d.pdf',count)),'file')==2
      count = count+1;
    end
    fn_out = strrep(fn_out,'.pdf',sprintf('_%d.pdf',count));
  end
  
  hfig = figure('units','pixel','outerposition',[0 0 1000 1414]);
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    if ~isfolder(fullfile(fp_d,pID))
      warning('ther is no rawdata folder for %s',pID)
      continue
    end
    fprintf('Plotting participant %d of %d : %s\n',i,size(ptab,1),pID);
    fl = ter_listFiles(fullfile(fp_d,pID),'*_physio.json');
    for j=1:numel(fl)
      [~,fn] = fileparts(fl{j});
%       sID = regexp(fl{j},'ses-[a-zA-Z0-9]*','match');
%       if isempty(sID)
%         sID = '';
%       else
%         sID = sID{end};
%       end
%       tID = regexp(fl{j},'task-[a-zA-Z0-9]*','match');
%       if isempty(tID)
%         tID = '';
%       else
%         tID = tID{end};
%       end
%       rID = regexp(fl{j},'run-[a-zA-Z0-9]*','match');
%       if isempty(rID)
%         rID = '';
%       else
%         rID = rID{end};
%       end
      jinfo = jsondecode(fileread(fl{j}));
      fn_tsv = strrep(fl{j},'.json','.tsv');
      if exist(fn_tsv,'file')==2
        data = readtable(fn_tsv,'filetype','text','delim','tab');
      elseif exist([fn_tsv '.gz'],'file')==2
        gunzip([fn_tsv '.gz']);
        data = readtable(fn_tsv,'filetype','text','delim','tab');
        delete(fn_tsv);
      end
      fn_events = strrep(fn_tsv,'physio.tsv','events.tsv');
      if exist(fn_events,'file')==2
        et = readtable(fn_events,'filetype','text','delim','tab');
      else
        et = table('size',[0,3],...
          'VariableNames',{'onset','duration','trial_type'},...
          'VariableTypes',{'double','double','cell'});
      end
      data.Properties.VariableNames(1:numel(jinfo.Columns)) = jinfo.Columns';
      myUnits = repmat({''},size(jinfo.Columns));
      for k=1:numel(jinfo.Columns)
        try 
          myUnits{k} = jinfo.(jinfo.Columns{k}).Units;
        catch
        end
      end
      
      data.time = (0:size(data,1)-1)' ./ jinfo.SamplingFrequency + jinfo.StartTime;
      myUnits =[ myUnits;'s']; %#ok<AGROW>
      data.Properties.VariableUnits = myUnits;
      events = unique(et.trial_type);
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
          ylabel(myLabel);
        else
          ylabel(sprintf('%s\n[%s]',myLabel,myUnit))
        end
        
        xlabel('t [min]')
        xlim(myXLimits)
      end
      
      if numel(events)>0
        subplot(heightPlots,1,(heightPlots+1-heightEv):heightPlots);
        
        for k=1:numel(events)
          starts = et.onset(ismember(et.trial_type,events{k}));
          stops  = et.duration(ismember(et.trial_type,events{k})) + starts;
          
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
      suptitle(strrep(fn,'_','\_'))
      
      export_fig(fn_out,hfig,'-append')
      
      
    end
  end
  fprintf('All figures saved to:\n  %s\n',fn_out)
  close(hfig);
  
end
