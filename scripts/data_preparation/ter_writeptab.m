function ter_writeptab(ptab,fn_ptab)

writetable(ptab,fn_ptab,'filetype','text','delim','tab');
txt0 = fileread(fn_ptab);
txt = strrep(txt0,'NaN','n/a');
if not(isequal(txt0,txt))
  fid = fopen(fn_ptab,'w');
  fprintf(fid,txt);
  fclose(fid);
end