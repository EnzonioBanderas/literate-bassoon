function log = essbids_readNeuroBSPresentationStimLog(fn_log,delim)
%% ESSBIDS_READNEUROBSPRESENTATIONSTIMLOG
%
%
% Reads one  NeuroBS Presentation stimulation log file and splits it in its
% up to 5 separate elements:
%
%  - GeneralEventList (table)
%  - StimulusList (table)
%  - VideoList (table) (to Do, currently no log avaiable with this list)
%  - warnings
%  - metadata (json decoded)
%
% Output: workspace variable
%
% Inputs: 
% fn_log : file name of stimulus log file (windows formated, i.e. lines 
%          end with CRLF
% delim  : optional delimiter input, if left empty delimiter will be
%          automatically detected
%
% Writes temporary files of selected content and utilized readtable 
% to create tables.
%
% Example:
%  log = essbids_readNeuroBSPresentationStimLog(fn_log)
%  >> Autodetection identified NeuroBS Presentation log delimiter: \t
%  >> log = 
%  >>   struct with fields:
%  >>             metadata: [1×1 struct]
%  >>     GeneralEventList: [2901×13 table]
%  >>
%
% Changelog: 
%   2022-06-21 : modified readtable options
%   2022-06-17 : first version
%

  % NeuroBS presentation allows these eight delimiter characters
  potential_delims = {'\t','!','"','#','$','%','&',''''};
  
  % parse delimiter input, if empty autodetect delimiter later
  if nargin < 2
    delim = '';
  elseif isempty(delim)
    delim = '';
  else
    if ~ischar(delim) 
      error('essbids:InvalidParameter',...
        'Delimiter must be a single character');
    elseif strcmpi(delim,'tab')
      delim = '\t';
    end
    if ~ismember(delim,potential_delims)
      fprintf('Delimiter must be one of these characters:\n%s',...
        evalc('disp(potential_delims'')'))
      error('essbids:InvalidParameter',...
        'Delimiter parameter invalid: %s',delim);
    end
  end
  
  % read content of file, split of metadata section (json formated)
  txt = fileread(fn_log);
  txt = essbids_substituteChars(txt);
  ind_meta = regexp(txt,'\r\n{\r\n');
  if isempty(ind_meta) 
    ind_meta = numel(txt);
  else
    log.metadata   = jsondecode(txt(ind_meta(1)+2:end));
  end
  txt_main = txt(1:ind_meta(1));
  txt_split = strsplit(txt_main, '\n');
  txt_split{1} = strsplit(txt_split{1}, ' - ');
  txt_split{2} = strsplit(txt_split{2}, ' - ');
  log.metadata.(matlab.lang.makeValidName(txt_split{1}{1})) = txt_split{1}{2};
  log.metadata.(matlab.lang.makeValidName(txt_split{2}{1})) =datetime(txt_split{2}{2});
  
  % auto-detect
  if isempty(delim)
    count = nan(size(potential_delims));
    for i=1:numel(potential_delims)
      count(i) = sum(ismember(txt_main,potential_delims{i}));
    end
    ind = find(count == max(count));
    if numel(ind)>1
      error(['Multiple potential delimiter characters found.\n'...
        'Please specify which one to use'])
    end
    delim=potential_delims{ind};
    delimOut = delim;
    if regexp(delim,'\t')
      delimOut = '\\t';
    end
    fprintf(['Autodetection identified NeuroBS Presentation '...
      'log delimiter: %s\n'],delimOut)
  end

  
  % split blocks dependend on empty lines
  ind_emptylines = regexp(txt_main,'\r\n\r\n');
  % proper blocks are setup in this order :empty line, column titles,
  % empty line, table, empty line
  for i=1:3:numel(ind_emptylines)-2
    h1 = txt_main(ind_emptylines(i)+4:ind_emptylines(i+1)-1);
    h1 = strrep(h1,' ','');
    h1 = regexprep(h1,['TTime' delim 'Uncertainty'],...
      ['TTime' delim 'TimeUncertainty']);
    h1 = regexprep(h1,['Duration' delim	'Uncertainty'],...
      ['Duration' delim 'DurUncertainty']);
    c1 = txt_main(ind_emptylines(i+1)+4:ind_emptylines(i+2)-1);
    fp_tmp = essbids_makeTempDir;
    fn_tmp = fullfile(fp_tmp,'tmp.csv'); 
    fid = fopen(fn_tmp,'w');
    %fprintf(fid,'%s\r\n%s',h1,c1);
    fprintf(fid,c1);
    fclose(fid);
    
    opts = detectImportOptions(fn_tmp,'delim',delim);
    opts.DataLines = [1,Inf];
    opts.VariableNames=strsplit(h1,delim);
    ind_code = ismember(opts.VariableNames,'Code');
    ind_char = ismember(opts.VariableTypes,'char');
    opts.VariableTypes(ind_char | ind_code)={'categorical'};
    t1 = readtable(fn_tmp,opts);
    rmdir(fp_tmp,'s');
    %t1.Properties.VariableNames=strsplit(h1,delim);
    if ismember(t1.Properties.VariableNames{1},{'Subject','Trial'})
      log.GeneralEventList = t1;
    elseif ismember(t1.Properties.VariableNames{1},{'EventType'})
      log.StimulusList = t1;
    else
      log.VideoList = t1;
    end
  end

  if mod(numel(ind_emptylines),3)>0
    log.Warnings = txt_main(ind_emptylines(end-1)+4:ind_emptylines(end)-1);
    log.Warnings = regexprep(log.Warnings,'\r','');
  end
  
  if isfield(log,'Warnings')
    warning('essbids:PresentationLogWarnings',...
      'In file\n  %s\n%s',fn_log,log.Warnings);
  end
  
