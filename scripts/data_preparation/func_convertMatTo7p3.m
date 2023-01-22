function func_convertMatTo7p3(fp_s)

fprintf('\n%s\n%s - Converting .mat files to Version 7.3\n',...
  repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))

fl_mat = ter_listFiles(fp_s,'*.mat');
ind0   = ~(cellfun(@ter_getMatFileVersion,fl_mat)>7);

fl_nonZip = fl_mat(ind0);


fn_temp = fullfile(fp_s,'temp_convertFunction.m');
if exist(fn_temp,'file')==2
  delete(fn_temp)
end

try
  tmp=cellfun(@dir,fl_nonZip);
  size_before = sum(arrayfun(@(x) x.bytes,tmp));
catch
  size_before = nan;
end
  

% run a matlab function loading and saving it's own little workspace
for i=1:numel(fl_nonZip)
  fprintf('Converting to matlab file version v7.3: %d of %d : %s\n',...
    i,numel(fl_nonZip),fl_nonZip{i});
  fid = fopen(fn_temp,'w');
  fprintf(fid,['function temp_convertFunction\n'...
    'load(''%s'');\n'...
    'save(''%s'',''-v7.3'');\n'], fl_nonZip{i},fl_nonZip{i});
  fclose(fid);
  run(fn_temp);
  delete(fn_temp)
end

try
  tmp=cellfun(@dir,fl_nonZip);
  size_after = sum(arrayfun(@(x) x.bytes,tmp));
catch
  size_after = nan;
end

if ~isnan(size_before) && ~isnan(size_after)
  fprintf(['%d mat files compressed, average compression %2.0f %%,'...
    ' %d MB saved\n'],...
    numel(fl_nonZip),size_after/size_before*100,...
    round((size_before-size_after)/2^20));
end