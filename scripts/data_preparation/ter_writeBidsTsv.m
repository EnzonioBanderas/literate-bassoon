function ter_writeBidsTsv(myTable,fn_table)
% 2021-10-21 : allow for gz extension

fn_out = fn_table;
[fp,fn,fe] = fileparts(fn_table);
if strcmpi(fe,'.gz')
  fn_out = fullfile(fp,fn);
end

% replace empty cells with string "n/a"
vn = myTable.Properties.VariableNames;
for i=1:numel(vn)
  if iscell(myTable.(vn{i}))
    ind = cellfun(@isempty,myTable.(vn{i}));
    myTable.(vn{i})(ind) = {'n/a'};
  end
end

writetable(myTable,fn_out,'filetype','text','delim','tab');
txt0 = fileread(fn_out);
txt = strrep(txt0,'NaN','n/a');
if not(isequal(txt0,txt))
  fid = fopen(fn_out,'w');
  fprintf(fid,txt);
  fclose(fid);
end

if strcmpi(fe,'.gz')
  gzip(fn_out);
  delete(fn_out);
end