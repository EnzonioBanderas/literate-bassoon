function PackagesFolder = func_getPackagesFolder()
%func_getPackagesFolder Get packages folder by searching for startspm package
%   might be better if it just actually searches for packages folders, but
%   this should work as long as startspm is added to the path
    fl_path = strsplit(path, ':');
    fl_logical = cellfun(@(y) strcmp(y{end}, 'startspm'), cellfun(@(x) strsplit(x, filesep), fl_path, 'UniformOutput', false), 'UniformOutput', false);
    PackagesFolder = fileparts(fl_path{find([fl_logical{:}], 1)});
end

