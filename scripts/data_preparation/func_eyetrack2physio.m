function func_eyetrack2physio(fp_s,fp_d)

  fprintf('\n%s\n%s : Converting eyetracking data and add to physio data\n',...
    repmat('=',1,72),datestr(now,'yyyy-mm-dd hh:MM:ss'));
  
  fn_ptab = fullfile(fp_d,'participants.tsv');
  ptab = ter_readptab(fn_ptab);
  
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    fl = ter_listFiles(fullfile(fp_s,pID),[pID '_ses-*_run-*_eyetrack.asc']);
    finfo = cellfun(@ter_parseFname,fl,'uni',0);
    for j=1:numel(fl)
      fn_asc = fl{j};
      sID = finfo{j}.ses;
      rID = finfo{j}.run;
      fn_events = ter_listFiles(fullfile(fp_d,pID,sID),...
        [pID '*_' sID '*_' rID '*_events.tsv']);
      if numel(fn_events)>1
        error('multiple events.tsv possible : %s',fn_asc)
      elseif numel(fn_events)==0
        warning('no events table found: %s %s %s',pID,sID,rID)
        continue
      end
      fn_events = fn_events{1};
      fn_new = strrep(fn_events,'_events.tsv','_eyetrack');
      fn_tsv = strcat(fn_new,'.tsv');
      fn_json = strcat(fn_new,'.json');
      
      if exist(fn_tsv,'file')==2 || exist(strcat(fn_tsv,'.gz'),'file')==2 
        continue
      end
      fprintf(' %s : Processing %s\n',datestr(now,'hh:MM:ss'),fn_asc);
      txtET = fileread(fn_asc);
  
      % example: "SAMPLES	GAZE	LEFT	RIGHT	RATE	 500.00"
      tmp = regexp(txtET,'SAMPLES\t\s*GAZE[A-Z\s]*(?<sf>[0-9\.]*)',...
        'names');
      if numel(tmp)>1
        warning('Multiple lines of sampling frequency detected, please check wether there is a problem')
      end
      et_info.SamplingFrequency = str2double(tmp(end).sf);
      et_info.StartTime = nan;
      et_info.Manufacturer = "SR-Research";
      et_info.ManufacturerModelName = "Eyelink 1000";
      et_info.Columns = {
        'position_horizontal_left';
        'position_vertical_left';
        'pupilsize_left';
        'position_horizontal_right';
        'position_vertical_right';
        'pupilsize_right';
        'flags'
        'fixation_left'
        'fixation_right'
        'sacchade_left'
        'sacchade_right'
        'blink_left'
        'blink_right'
        };

      tmp = regexp(txtET,'START\t(?<starttime>[0-9]*)','names');
      if numel(tmp)>1
        warning('Multiple lines of start time detected, please check wether there is a problem')
      end
      starttime = str2double(tmp(end).starttime);
      tmp = regexp(txtET,['\n'...
        '(?<time>[0-9]*)\t\s*'...
        '(?<xpl>[\-0-9\.]*)\t\s*'...
        '(?<ypl>[\-0-9\.]*)\t\s*'...
        '(?<psl>[\-0-9\.]*)\t\s*'...
        '(?<xpr>[\-0-9\.]*)\t\s*'...
        '(?<ypr>[\-0-9\.]*)\t\s*'...
        '(?<psr>[\-0-9\.]*)\t\s*'...
        '(?<flags>[I\.][C\.][R\.][C\.][R\.])'],'names');
      ettab = struct2table(tmp);
      f2d = {'time','xpl','ypl','psl','xpr','ypr','psr'};
      for k=1:numel(f2d)
        ettab.(f2d{k}) = str2double(ettab.(f2d{k}));
      end
  
      %a=regexp(txt0,'SFIX');
      %b=regexp(txt0,'EFIX');
      %if numel(a)~= numel(b)    
      %  warning('number of fixation starts and ends not equal')
      %end
  
      ev2query = {
        'fixL' 'EFIX\s*L\s*(?<tstart>[\-0-9\.]*)\s*(?<tstop>[\-0-9\.]*)';
        'fixR' 'EFIX\s*R\s*(?<tstart>[\-0-9\.]*)\s*(?<tstop>[\-0-9\.]*)';
        'sacL' 'ESACC\s*L\s*(?<tstart>[\-0-9\.]*)\s*(?<tstop>[\-0-9\.]*)';
        'sacR' 'ESACC\s*R\s*(?<tstart>[\-0-9\.]*)\s*(?<tstop>[\-0-9\.]*)';
        'blkL' 'EBLINK\s*L\s*(?<tstart>[\-0-9\.]*)\s*(?<tstop>[\-0-9\.]*)';
        'blkR' 'EBLINK\s*R\s*(?<tstart>[\-0-9\.]*)\s*(?<tstop>[\-0-9\.]*)';
      };
      for k=1:size(ev2query,1)
        tmp = struct2table(regexp(txtET,ev2query{k,2},'names'));
        tmp.tstart = str2double(tmp.tstart); 
        tmp.tstop  = str2double(tmp.tstop);
        ettab.(ev2query{k,1}) = zeros(size(ettab,1),1);
        for l=1:numel(tmp.tstart)
          ind = ettab.time >= tmp.tstart(l) & ettab.time <= tmp.tstop(l);
          ettab.(ev2query{k,1})(ind)=1;
        end
      end
      %ettab.time = ettab.time./1e3; % transfer time in seconds
  
      % gather the events
      tmp = regexp(txtET,...
        '\nMSG	(?<onset>[0-9]*) TRIALID TRIAL (?<trialindex>[0-9]*)',...
        'names');
      events = struct2table(tmp);
      events.onset = str2double(events.onset)/1e3;
      tmp = regexp(txtET,'\nMSG	(?<trialendtime>[0-9]*) SYNCTIME_End',...
        'names');
      tmp = struct2table(tmp);
      tmp.trialendtime = str2double(tmp.trialendtime)/1e3;  
      events.duration = (tmp.trialendtime-events.onset);

      et = readtable(fn_events,'filetype','text','delim','tab',...
        'treat','n/a');
      et1 = et(ismember(et.trial_type,'Trial'),:);
  
      
      % now calculate relative time in between the events time and the et
      % time
      relTimes = events.onset-starttime/1e3-et1.onset;
      resTemp = 1/et_info.SamplingFrequency; %temporal resolution in seconds
      relTime = mean( relTimes((relTimes>= median(relTimes)-5* resTemp) & ...
        (relTimes<= median(relTimes)+5* resTemp)));
      relTime = round(relTime/resTemp)*resTemp;
      et_info.StartTime = -relTime;
      fprintf('writting file : %s\n',fn_json);
      ter_savePrettyJson(fn_json,et_info);
       
      fprintf('writting file : %s.gz\n',fn_tsv);
      writetable(ettab(:,2:end),fn_tsv,'fileType','text',...
        'delim','tab','writeVar',0);
      txt1 = fileread(fn_tsv);
      txt2 = strrep(txt1,'NaN','n/a');
      if not(isequal(txt1,txt2))
        fid=fopen(fn_tsv,'w');
        fprintf(fid,txt2);
        fclose(fid);
      end
      % finally gzip
      gzip(fn_tsv);
      delete(fn_tsv);
      
    end
  end
end