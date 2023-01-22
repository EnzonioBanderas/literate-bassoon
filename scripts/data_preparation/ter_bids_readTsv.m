function [data,timevector] = ter_bids_readTsv(fn_tsv)
  %% TER_BIDS_READTSV
  %  
  % Read TSV file to matlab table object. If the json sidecars contains
  % column names these will be used. If the json sidecar contains a
  % StartTime and SamplingFrequency field a time vector will be generated 
  % for each point in time, else timevector output will be empty. Does 
  % allow for gziped input files.
  %
  % Examples:
  %  
  %  [datatable,timevector] = ter_bids_readTsv('path/to/myfile.tsv');
  %  [datatable,timevector] = ter_bids_readTsv('path/to/myfile.tsv.gz');
  %
  % Changes:
  %  2021-12-02 Renamed file, allow cell inputs, utilize %USER% path for
  %             temporary storage of unzipped tsv files in a uniquely named
  %             temporary folder 
  %  2021-11-09 Unzip into scripts folder, not in same directory as source
  %             file
  %  2021-07-01 Added units to output table
  %
  
  if iscell(fn_tsv)
    data = cell(size(fn_tsv));
    timevector = cell(size(fn_tsv));
    for i=1:numel(fn_tsv)
      [tmp1,tmp2] = ter_bids_readTsv(fn_tsv{i});
      data{i}=tmp1;
      timevector{i}=tmp2;
    end
    return;
  elseif ~ischar(fn_tsv)
    disp('Invalid Input:');
    disp(fn_tsv);
    error('input must be string file name or cell array of file names')
  elseif exist(fn_tsv,'file')~=2
    error('file not found: %s',fn_tsv)
  end
  
  [fp,fn,fe] = fileparts(fn_tsv);
  if strcmpi(fe,'.gz')
    [~,fn,fe2] = fileparts(fn);
    fe = strcat(fe2,fe);
  end
  
  fn_json = fullfile(fp,strcat(fn,'.json'));
  if exist(fn_json,'file')==2
    jinfo = jsondecode(fileread(fn_json));
  else
    jinfo = struct([]);
  end
  
  tsv_readparam = {'filetype','text','delim','tab','treat','n/a'};
  if strcmpi(fe,'.tsv.gz')
    fn2 = fullfile(fp,strcat(fn,fe(1:4)));
    if exist(fn2,'file')==2
      data = readtable(fn2,tsv_readparam{:});
    else
      % create temporary folder in user directory (you should always have
      % writing rights there) and delete it after unzip
      if ispc()
        userDir = getenv('USERPROFILE');
      else
        userDir= getenv('HOME');
        %char(java.lang.System.getProperty('user.home'));
      end
      % create (virtually) uniquely named temporary directory 
      fp_tmp = fullfile(userDir,sprintf('temp_%s_%09.0f',...
        datestr(now,'yyyy-mm-dd-HH-MM-SS.FFF'),floor(rand*(1e9-1))));
      mkdir(fp_tmp);
      %fp_tmp = fileparts(mfilename('fullpath'));
      fn_tmp = gunzip(fn_tsv,fp_tmp);
      fn_tmp = fn_tmp{1};
      data = readtable(fn_tmp,tsv_readparam{:});
      delete(fn_tmp);
      rmdir(fp_tmp);
    end 
  else
    data = readtable(fn_tsv,tsv_readparam{:});
  end
  
  % add column names
  if isfield(jinfo,'Columns')
    data.Properties.VariableNames = jinfo.Columns;
  end
  
  % add units if available
  if isfield(jinfo,'Columns')
    myUnits = repmat({''},size(jinfo.Columns));
    for k=1:numel(jinfo.Columns)
      try 
        myUnits{k} = jinfo.(jinfo.Columns{k}).Units;
      catch
      end
    end
  else
    vn = data.Properties.VariableNames;
    myUnits = repmat({''},size(vn));
    for k=1:numel(vn)
      try 
        myUnits{k} = jinfo.(vn{k}).Units;
      catch
      end
    end
  end
  data.Properties.VariableUnits = myUnits;
  
  % generate time vector
  if isfield(jinfo,'StartTime') && isfield(jinfo,'SamplingFrequency')
    timevector = (1:size(data,1))'./jinfo.SamplingFrequency +  ...
      jinfo.StartTime;
  else
    timevector = [];
  end
  
end