function ter_savePrettyJson(jsonFilename,struct2convert)
% TER_SAVEPRETTYJSON
%   nicely formats an existing json file, or converts a structure to json 
%   text and saves that nicely formated to a file
%
%   calls a python function: python3 -m json.tool
%

  % Umlaut etc. conversion
  str2rep = {
    'ä' 'ae'
    'Ä' 'Ae'
    'ö' 'oe'
    'Ö' 'Oe'
    'ü' 'ue'
    'Ü' 'Ue'
    'ß' 'ss'
    '?' ''''
    };

  if exist(jsonFilename,'file')==2 && nargin < 2
    % nothing to do really
  else
    strTemp = jsonencode(struct2convert);
    % replace generic placeholder x of numeric fieldnames ,
    % e.g. rename field "x1" to "1" in json text
    for i=1:size(str2rep,1)
      strTemp = strrep(strTemp,str2rep{i,1},str2rep{i,2});
    end
    fn2repl = unique(struct2cell(...
      regexp(strTemp,'(?<fn>"x_[0-9]+")','names')));
    for i=1:numel(fn2repl)
      strTemp = strrep(strTemp,fn2repl{i},...
        [fn2repl{i}(1),'-',fn2repl{i}(4:end)]);
    end
    fn2repl = unique(struct2cell(...
      regexp(strTemp,'(?<fn>"x[0-9]+")','names')));
    for i=1:numel(fn2repl)
      strTemp = strrep(strTemp,fn2repl{i},fn2repl{i}([1,3:end]));
    end
    %strTemp = strrep(strTemp,'%','%%'); % fprintf stumbles on single '%'
    fid     = fopen(jsonFilename,'w');
    fprintf(fid,'%s',strTemp);
    fclose(fid);
  end
  
  if ispc
    pyver = evalc('!py --version');
  else
    pyver = evalc('!python3 --version');
  end
  if isempty(regexp(pyver,'Python 3\.', 'once'))
    warning('python3 missing, please install')
    return
  end
  
  if isunix
    myCmd = sprintf('!python3 -m json.tool < %s',...
        strrep(jsonFilename,' ','\ ')); 
    % strrep should take care of spaces in file name or path
  elseif ispc
    myCmd = sprintf('!py -m json.tool < %s',...
      strrep(jsonFilename,' ','\ ')); 
  end
  
  strTemp = evalc(myCmd);
  %strTemp = strrep(strTemp,'%','%%');
  
  fid     = fopen(jsonFilename,'w');
  %fprintf(fid,strTemp);
  fprintf(fid,'%s',strTemp);
  fclose(fid);

end