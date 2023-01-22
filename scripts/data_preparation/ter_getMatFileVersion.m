function matFileVersion = ter_getMatFileVersion(matfilename)
  mymsg = getMatComment(matfilename);
  
  mymsg = strsplit(mymsg,' ');
  
  if numel(mymsg) < 2 
    error('function type cannot identify file');
  elseif isempty(strfind(mymsg{1}, 'MATLAB'))
    error('function type output weird');
  end
    
  matFileVersion = str2double(mymsg{2});
end

function txt=getMatComment(x)
  fid= fopen(x);
  txt= fread(fid,[1,140],'*char');
  txt= [txt,0];
  txt= txt(1:find(txt==0,1,'first')-1);
  fclose(fid);
end

% alternative approach

% function matFileVersion = ter_getMatFileVersion(matfilename)
%   if ~exist(matfilename,'file')
%     error('file does not exist')
%   end
%   
%   mymsg = evalc(sprintf('type %s',matfilename));
%   
%   % expected start of output "MATLAB 7.3 MAT-file" or similar
%   
%   % number of chars at the start of the file that might be of interest
%   numCharOfInterest = 20;
%   
%   if numel(mymsg)>numCharOfInterest
%     mymsg = mymsg(1:numCharOfInterest);
%   end
%   mymsg = strsplit(mymsg,' ');
%   
%   if numel(mymsg) < 2 
%     error('function type cannot identify file');
%   elseif isempty(strfind(mymsg{1}, 'MATLAB'))
%     error('function type output weird');
%   end
%     
%   matFileVersion = str2double(mymsg{2});
% 
% end
% 
