function ptab = ter_readptab(fn_ptab)

  ptab = readtable(fn_ptab,'filetype','text','delim','tab','treat','n/a',...
    'datetime','text');
  if ismember('sex',ptab.Properties.VariableNames)
    if isnumeric(ptab.sex)
      if prod(isnan(ptab.sex))==1
        ptab.sex = repmat({'n/a'},size(ptab,1),1);
      end
    end
  end
  
end
