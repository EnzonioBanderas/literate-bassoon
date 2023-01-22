function out = func_checkForPython3Availability

%% check for availability of python3 
% just needed for nicely formated json files, nothing more right now
if ispc
  pyver = evalc('!py --version');
else
  pyver = evalc('!python3 --version');
end
if isempty(regexp(pyver,'Python 3\.', 'once'))
  warning('python3 missing, please install fot nicely formated json files')
  out = false;
end

pause(3);
