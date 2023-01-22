function  [logdata,fname_out_mat]  = ter_readPresentationLog(varargin)
% TER_READPRESENTATIONLOG
%
%  Reads all the potential presentation logs (.log, .txt, .csv, .dat) in 
%  one folder and all subfolder levels below. 
%  Returns a structured dataset "logdata" and saves it in a .mat-file.
%  structure fields will be named as filenames of logs 
%
%  Will handle events tables, stimulation tables and questionnaire response
%  files (formated as Giorgi's NeuroBS Presentation questionnaire response
%  txt files) 
%
%   @author:      Thomas M. Ernst
%   @institute:   University Clinic Essen, Motor Control Lab
%   @date:        2020-12-29
%   @version:     0.4
%
%  Changelog: 
%
%    2020-12-29 : added questionnaire readout of datetime, participant and
%                 protocol, added "fname_out_mat" to outputs
%


  default_source = '';
  default_target = '';
  default_select = '*';

  skipExisting = ismember('skipexisting',lower(varargin));
  beQuiet  = ismember('quiet',lower(varargin));
  varargin = varargin(~ismember(lower(varargin),{'skipexisting','quiet'}));
  p = inputParser;
  addParameter(p, 'source', default_source, @isfolder);
  addParameter(p, 'target', default_target, @isfolder);
  addParameter(p, 'select', default_select, @ischar);

  parse(p,varargin{:});

  source          = p.Results.source;
  target          = p.Results.target;
  select          = p.Results.select;

  count = 0;
  while ~isfolder(source) && count < 3
    source = uigetdir(pwd,'source');
    count = count + 1;
  end
  count = 0;
  while isfolder(source) && ~isfolder(target) && count < 3
    target = uigetdir(source,'target');
    count = count + 1;
  end
  
  % NeuroBS presentation identifiers
  scenarioID     = 'Scenario - ';
  logdatetime    = 'Logfile written - ';
  % questionnaire identifiers
  qID_sub      = 'ParticipantID         : ';
  qID_start    = 'Protocol start time   : ';
  qID_protocol = 'Protocol phase        : ';

  questListStart = 'TrialN';
  eventListStart = {'Trial'; sprintf('Subject\tTrial')};
  stimListStart  = 'Event Type';
  logfileExt     = {'.csv','.dat','.log','.txt'};
 
  selectFilter   = cellfun(@(x) [select x],logfileExt,...
    'UniformOutput',false);
  forbidden_char = {
    '\' ''  ; '/' ''  ; '%' '_' ; '$' '_'; '&' '_'; 
    ' ' ''  ; '?' '_' ; '#' '_' ; '-' '_'; '@' '_at_';
    '?' 'ae'; '?' 'ue'; '?' 'oe';
    '?' 'Ae'; '?' 'Ue'; '?' 'Oe'; };
  
  logfilelist = ter_listFiles(source,selectFilter,0);
  %logfilelist = ter_listAllExtension(source,logfileExt,0,0);
  logfilelist = logfilelist(~cellfun(@isempty,logfilelist));
%   dirname = strsplit(source,filesep);
%   dirname = dirname{end};
  if ~beQuiet
    fprintf(['\nNumber of potential presentation log-files found in\n' ...
      '''%s'':\n  %d\n'],source,numel(logfilelist));
  end
  
  fname = strsplit(source,filesep);
  fname_out = fname{end};
  fname_out = ter_replaceForbiddenChar(fname_out,forbidden_char(:,1),...
    forbidden_char(:,2));
  % tabext = '.xlsx';
  % fname_out_tab = fullfile(target,[fname_out tabext]);
  fname_out_mat = fullfile(target,'mat_files',[fname_out '.mat']);
  if skipExisting && exist(fname_out_mat,'file')
    if ~beQuiet
      fprintf('mat-file already exists - skipping preslog2mat conversion\n');
    end
    %logdata = [];
    load(fname_out_mat,'logdata');
    return
  end
  % fname_out_tab = ter_checkFilenameExistance(fname_out_tab);
  fname_out_mat = ter_checkFilenameExistance(fname_out_mat);
  datasetname = cell(numel(logfilelist),1);
  for i=1:numel(logfilelist)
    % extract dataset name from log file name
    [~,logfname0,~] = fileparts(logfilelist{i});
    logfname = logfname0;
    count = 1;
    %avoid doubles
    while ismember(logfname,datasetname(1:(i-1)))
      logfname = sprintf('%s_%1.0f',logfname0,count);
      count    = count+1;
    end
    datasetname{i} = ter_replaceForbiddenChar(logfname,...
      forbidden_char(:,1),forbidden_char(:,2));
    
    if ~isnan(str2double(datasetname{i}(1)))
      datasetname{i} = ['file_' datasetname{i}];
    end
    whichList = 'none';
    
    disp(logfilelist{i})
    
    % now read dataset line by line
    fid   = fopen(logfilelist{i},'r');
    tline = fgetl(fid);
    while ~isequal(tline, -1)
      if isempty(tline)
        % skip it
      elseif strncmpi(scenarioID,tline,numel(scenarioID))
        res.(datasetname{i}).scenarioID = tline(numel(scenarioID):end );
      elseif strncmpi(logdatetime,tline,numel(logdatetime))
        try 
          res.(datasetname{i}).datetime   = datestr(datevec(...
            tline(numel(logdatetime):end)),'yyyy-mm-dd HH:MM:SS');
        catch
          res.(datasetname{i}).datetime   = tline(numel(logdatetime):end);
        end
      elseif strncmpi(questListStart,tline,numel(questListStart))
        my_delimiter      = tline(numel(questListStart)+1);
        my_questListName  = strsplit(tline,my_delimiter);
        for j=1:numel(my_questListName)
          for k=1:size(forbidden_char,1)
            my_questListName{j} = strrep(my_questListName{j},...
              forbidden_char{k,1},forbidden_char{k,2});
          end
        end
        my_questListName = ter_renameDoubleVarNames(my_questListName);
        res.(datasetname{i}).questList = cell(0,numel(my_questListName));
        whichList = 'quest';
      elseif strncmpi(eventListStart{1},tline,numel(eventListStart{1}))
        my_delimiter = tline(numel(eventListStart{1})+1);
        my_eventListName  = strsplit(tline,my_delimiter);
        for j=1:numel(my_eventListName)
          for k=1:size(forbidden_char,1)
            my_eventListName{j} = strrep(my_eventListName{j},...
              forbidden_char{k,1},forbidden_char{k,2});
          end
        end
        my_eventListName = ter_renameDoubleVarNames(my_eventListName);
        res.(datasetname{i}).eventList = cell(0,numel(my_eventListName));
        whichList = 'event';
      elseif strncmpi(eventListStart{2},tline,numel(eventListStart{2}))
        my_delimiter = tline(numel(eventListStart{2})+1);
        my_eventListName  = strsplit(tline,my_delimiter);
        for j=1:numel(my_eventListName)
          for k=1:size(forbidden_char,1)
            my_eventListName{j} = strrep(my_eventListName{j},...
              forbidden_char{k,1},forbidden_char{k,2});
          end
        end
        my_eventListName = ter_renameDoubleVarNames(my_eventListName);
        res.(datasetname{i}).eventList = cell(0,numel(my_eventListName));
        whichList = 'event';
      elseif strncmpi(stimListStart,tline,numel(stimListStart))
        my_delimiter    = tline(numel(stimListStart)+1);
        my_stimListName = strsplit(tline,my_delimiter);
        for j=1:numel(my_stimListName)
          for k=1:size(forbidden_char,1)
            my_stimListName{j} = strrep(my_stimListName{j},...
              forbidden_char{k,1},forbidden_char{k,2});
          end
        end
        my_stimListName = ter_renameDoubleVarNames(my_stimListName);
        res.(datasetname{i}).stimList  = cell(0,numel(my_stimListName));
        whichList = 'stim';
      elseif strcmpi(whichList,'quest')
        mynewline = cell(1,numel(my_questListName));
        mynewcontent = strsplit(tline,my_delimiter);
        mynewline(1:numel(mynewcontent)) = mynewcontent;
        res.(datasetname{i}).questList = [res.(datasetname{i}).questList;
          mynewline];
      elseif strcmpi(whichList,'event')
        mynewline = cell(1,numel(my_eventListName));
        mynewcontent = strsplit(tline,my_delimiter);
        mynewline(1:numel(mynewcontent)) = mynewcontent;
        res.(datasetname{i}).eventList = [res.(datasetname{i}).eventList;
          mynewline];
      elseif strcmpi(whichList,'stim')
        mynewline = cell(1,numel(my_stimListName));
        mynewcontent = strsplit(tline,my_delimiter);    
        mynewline(1:numel(mynewcontent)) = mynewcontent;
        res.(datasetname{i}).stimList  = [res.(datasetname{i}).stimList;
          mynewline];
      end
      tline = fgetl(fid);
    end
    fclose(fid);
    if ~isfield(res,datasetname{i})
      fprintf('Not a NEUROBS Presentation log file:\n   %s\n',...
        logfilelist{i});
      continue
    end
    myfields = fieldnames(res.(datasetname{i}));
    if ismember('questList',myfields)
      % get further information
      txt = fileread(logfilelist{i});
      tmp = regexp(txt,[qID_sub '(?<pID>sub-[0-9a-zA-Z]*)'],'names');
      res.(datasetname{i}).participant_id = tmp.pID;
      tmp=regexp(txt,[qID_start '(?<datetime>[0-9a-zA-Z\ \-\:]*)'],'names');
      try
        my_datetime = datetime(tmp.datetime,'Locale','en_us');
      catch
        try
          my_datetime = datetime(tmp.datetime,'Locale','de_de');
        catch
          my_datetime = tmp.datetime;
        end
      end
      try 
        my_datetime = datestr(my_datetime,'yyyy-mm-dd HH:MM:ss');
      catch
      end
      res.(datasetname{i}).datetime = my_datetime;
      tmp = regexp(txt,[qID_protocol '(?<protocol>[\w\ \-\.\,\:\;\(\)\[\]]*)'],'names');
      res.(datasetname{i}).protocol = tmp.protocol;
      
      res.(datasetname{i}).questList = cell2table(...
        res.(datasetname{i}).questList);
      res.(datasetname{i}).questList.Properties.VariableNames = ...
        my_questListName;
      writetable(res.(datasetname{i}).questList, 'temptable.txt');
      res.(datasetname{i}).questList = readtable('temptable.txt');
      delete('temptable.txt');
      
    end
    if ismember('eventList',myfields)
      res.(datasetname{i}).eventList = cell2table(...
        res.(datasetname{i}).eventList);
      res.(datasetname{i}).eventList.Properties.VariableNames = ...
        my_eventListName;
      writetable(res.(datasetname{i}).eventList, 'temptable.txt');
      res.(datasetname{i}).eventList = readtable('temptable.txt');
      delete('temptable.txt');
    end
    if ismember('stimList',myfields)
      res.(datasetname{i}).stimList = cell2table(...
        res.(datasetname{i}).stimList);
      res.(datasetname{i}).stimList.Properties.VariableNames  = ...
        my_stimListName;
      writetable(res.(datasetname{i}).stimList, 'temptable.txt');
      res.(datasetname{i}).stimList = readtable('temptable.txt');
      delete('temptable.txt');
    end
  end
  if exist('res','var')
    logdata = res;
    if ~isfolder(fileparts(fname_out_mat))
      mkdir(fileparts(fname_out_mat))
    end
    save(fname_out_mat,'logdata');
  else
    logdata = [];
  end
end


function outStr = ter_replaceForbiddenChar(inStr,forbiddenChar,replaceChar)
  outStr = inStr;
  if iscell(inStr)
    for j=1:numel(outStr)
      for i=1:numel(forbiddenChar)
        outStr{j} = strrep(outStr{j},forbiddenChar{i},replaceChar{i});
      end
    end
  else
    for i=1:numel(forbiddenChar)
      outStr = strrep(outStr,forbiddenChar{i},replaceChar{i});
    end
  end  
end

function outVarName = ter_renameDoubleVarNames(inVarName)
  outVarName = inVarName;
  uniqueVarName = unique(inVarName);
  for i=1:numel(uniqueVarName)
    myIndex = ismember(inVarName,uniqueVarName{i});
    if sum(myIndex)>1
      myIndex = find(myIndex);
      for j=1:numel(myIndex)
        if myIndex(1)==1
          outVarName{myIndex(j)} = sprintf('%s%d',...
            outVarName{myIndex(j)},j);
        else
          outVarName{myIndex(j)} = sprintf('%s%s',...
            outVarName{myIndex(j)-1},outVarName{myIndex(j)});  
        end
      end
    end
  end
end




function list=ter_listAllExtension(parent,extension,subDirs,csens)
%TER_LISTALLEXTENSION
%   Generate a cell array of strings containing all filenames including 
%   paths. 
%   
%   @author:      Thomas M. Ernst
%   @institute:   University Clinic Essen, Motor Control Lab
%   @date:        2014-09-11
%   @version:     1.0
% 
%   last modified:    2014-09-11
% 
% 
%   LIST = TER_LISTALLEXTENSION(PARENT,EXTENSION,SUBDIRS,CSENS)
%     returns a cell array of strings containing the finename including
%     paths of all files in the PARENT folder matching the EXTENSION
%     statement. Subfolders up to the level specified in SUBDIRS are
%     searched. Search can be performed case sensitive or insensitive 
%     (CSENS).
% 
%   PARAMETERS:
%     parent      Parent directory to be searched.
%                 Accepts strings and 1-dim. cell arrays of strings, e.g.
%                 '/folder/subf' and {'/f1/f2/f3' '/f1/f4' '/f5'}
%                 NOTE: Double paths or partially overlapping directories
%                 will result in duplicate results.
%     extension   Extension to be searched, e.g. '.img'.
%                 Does accept 1-dim. cell arrays, too, e.g. {'.img' '.hdr'}
%                 Please remember to start each extension with a ".".
%     subDirs     Specify if subdirectories should be searched.
%                 Accepts numeric input (i.e. subfolder level) or 'all'
%                 Negative input is inverted to positive values,  
%                 decimals are discarded.
%     csens       specify if search should be performed case sensitively
%                 (i.e. e is not equal E, csens=1) or not (csens=0)
% 
%   
%--------------------------------------------------------------------------
%

list={}; %#ok<NASGU>

% Boring parameter checking and rectifing ...dumdidum...
stopit=0;
if not(iscell(parent))
  parent={parent};
end
if not(iscell(extension))
  extension={extension};
end
for i=1:length(parent)
  if not(exist(parent{i},'dir'))
    stopit=1;
    fprintf(['\nFolder missing: %s\n'...
      'Please enter only valid folder paths as first parameter.\n'],...
      char(parent(i))); 
  end
end
if isnumeric(subDirs)
  searchSubLvl=abs(subDirs);
  % fprintf('\nSearching in subfolders up to level: %d\n',searchSubLvl);
elseif strcmp(subDirs,'all') 
  searchSubLvl=Inf;
  % fprintf('\nSearching in all subfolders.\n');
else
  disp('Please enter positive interger or "all" as third parameter');
  stopit=1;
end
if isnumeric(csens)
  if not(csens==0 || csens==1)
    stopit=1;
    disp('Please enter 0 or 1 as fourth parameter (case sensitivity)');
  end
else
  disp('Please enter 0 or 1 as fourth parameter (case sensitivity).');
  stopit=1;
end
if stopit
  error('\nParameter Error. Stop function execution.\n');
end

% initialize result arrays
arrDir      = parent;
arrDirLvl   = zeros(length(arrDir),1);
arrFiles    = cell(1,1);

% Count number of directories und files seperately
countD=1;
countF=0;

% At this point the young padavan has learned that for-loops do not
% reevaluate the condition statement but while-loops do.
i=0;
while i<length(arrDir)
  i=i+1;
  p=arrDir{i};
  if isempty(p)
    continue
  end
  % remember current directories and subfolder level
  currDir=dir(char(p));
  currLvl=arrDirLvl(i);
  for j=1:length(currDir)
    filename = currDir(j).name;
    % discarding '.' and '..' 
    if  strcmp( filename,'.') || strcmp( filename,'..')
      continue;
    end
    % handle files
    if currDir(j).isdir == false
      [~, ~, ext] = fileparts(filename);
      for k=1:length(extension)
        testExtension=extension(k);
        if csens==1
          matched=strcmp(ext, testExtension);
        else
          matched=strcmpi(ext, testExtension);
        end
        if matched
          countF=countF+1;
          arrFiles(countF) = {fullfile(char(p),filename)};
          break
        end
      end
    % handle subfolders, i.e. ignore them or add them to the list of 
    % folders to search
    elseif searchSubLvl>arrDirLvl(countD)
      arrDir(countD+1)    = {fullfile(char(p),filename)};
      arrDirLvl(countD+1,1) = currLvl+1;
      countD=countD+1;
    end
  end
end

% tiny report to shell:
% fprintf([...
%   '\nSearched: %6d folders\nFound:    %6d files\n\n'],...
%   countD,countF);

% return the list of files
list = arrFiles;
return


end



function fnameNew = ter_checkFilenameExistance(fname0,varargin)
%TER_CHECKFILENAMEEXISTANCE
% Checks wether a file name exists and genrerates a new non-existing 
% filename. 
%      
%   @author:      Thomas M. Ernst
%   @institute:   University Clinic Essen, Motor Control Lab
%   @date:        2015-05-05
%   @version:     1.1
% 
%   last modified:     2016-01-06
% 
%   History:       1.1 2016-01-06
%                      allow input of additional search paths, maximum+1
%                      index will be returned
%                  1.0 2015-05-05
% Examples:
%  
%  (a) ter_checkFilenameExistance('C:/MyPath/myFile.txt')
%      returns
%        if the file does not exist: C:/MyPath/myFile.txt 
%        if the file exists:         C:/MyPath/myFile_01.txt 
%        if the file and '_01'exist: C:/MyPath/myFile_02.txt
%        ...
%  (b) ter_checkFilenameExistance('C:/MyPath/myFile_05.txt')
%      returns
%        if the file does not exist: C:/MyPath/myFile_05.txt 
%        if the file exists:         C:/MyPath/myFile_05_01.txt
%        ...
% =========================================================================
  
  additionalDirs = '';
  p = inputParser;
  addParameter(p, 'dirs',additionalDirs, ...
    @(x) iscell(x) && prod(cellfun(@ischar,x)));
  parse(p,varargin{:});
  additionalDirs = p.Results.dirs;
  
  if isempty(additionalDirs)
    fnameNew           = fname0;
    [fpath,fname,fext] = fileparts(fname0);
    counter = 1;
    while exist(fnameNew,'file')
      fnameNew = fullfile(fpath,sprintf('%s_%02.0f%s',fname,counter,fext));
      counter  = counter + 1;
    end
  else

    [fpath,fname,fext] = fileparts(fname0);
    fname              = [fname fext];
    additionalDirs     = [fpath additionalDirs];
    fnameList          = cell(numel(additionalDirs),1);
    for i=1:numel(additionalDirs)
      fnameList{i} = ter_checkFilenameExistance(...
        fullfile(additionalDirs{i},fname));
      [~,fnametmp,fexttmp] = fileparts(fnameList{i});
      fnameList{i} = [fnametmp fexttmp];
    end
    fnameLength = cellfun(@numel,fnameList);
    fnameList   = fnameList(fnameLength == max(fnameLength)); 
    if numel(fnameList)==1
      fnameNew = fullfile(fpath, fnameList{1});
      return;
    elseif isequal(fnameList{:})
      fnameNew = fullfile(fpath, fnameList{1});
      return;
    end
    fname_index = cell(numel(fnameList),1);
    for i=1:numel(fnameList)
      temp_index = strfind(fnameList{i},'_');
      fname_index{i} = ...
        str2double(fnameList{i}(temp_index(end)+1:end-numel(fext)));
      if isempty(fname_index{i})
        fname_index{i} = 0;
      end
    end
    fname_index = cell2mat(fname_index);
    final_index = find(fname_index == max(fname_index));
    fnameNew = fullfile(fpath, fnameList{final_index(1)});
    return;
  end
end
