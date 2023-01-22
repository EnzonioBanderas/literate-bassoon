function [data,timevector] = ter_readBidsTsv(fn_tsv)
  %% 
  % 2021-11-09 Unzip into scripts folder, not in same directory as source
  % file
  % 2021-07-01 Added units to output table
  if ~ischar(fn_tsv)
    error('input must be string file name')
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
      fp1 = fileparts(mfilename('fullpath'));
      fn_tmp = gunzip(fn_tsv,fp1);
      fn_tmp = fn_tmp{1};
      data = readtable(fn_tmp,tsv_readparam{:});
      delete(fn_tmp);
    end 
  else
    data = readtable(fn_tsv,tsv_readparam{:});
  end
  
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
  
  if isfield(jinfo,'Columns')
    data.Properties.VariableNames = jinfo.Columns;
  end
  if isfield(jinfo,'StartTime')
    timevector = (1:size(data,1))'./jinfo.SamplingFrequency +  ...
      jinfo.StartTime;
  else
    timevector = [];
  end
  
end