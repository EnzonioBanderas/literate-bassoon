t_demo = jsondecode(fileread('demographics.json'));
t_dass = jsondecode(fileread('dass21g.json'));
t_ehi  = jsondecode(fileread('ehi.json'));

clearvars t0
for h=1:3
  switch h
    case 1
      t00 = t_demo;
      ind0 = 2;
    case 2
      t00 = t_dass;
      ind0 = 8;
    case 3
      t00 = t_ehi;
      ind0 = 5;
  end
  label = fieldnames(t00);
  label = label(ind0:end);
  description = cell(size(label));
  levels      = repmat({''},size(label));
  value       = repmat({''},size(label));
  for i=1:numel(label)
    description{i} = t00.(label{i}).Description;
    if isfield( t00.(label{i}),'Levels')
      fn = fieldnames(t00.(label{i}).Levels);
      fn2 = fn;
      ind = cellfun(@(x) ~isempty(regexp(x,'^x\d','ONCE')),fn);
      fn2(ind) = cellfun(@(x) x(2:end),fn2(ind),'uni',0);
      str_lvl = sprintf('%s : %s',fn2{1},t00.(label{i}).Levels.(fn{1}));
      for j=2:numel(fn)
        str_lvl = sprintf('%s, %s : %s',str_lvl,fn2{j},...
          t00.(label{i}).Levels.(fn{j}));
      end
      levels{i} = str_lvl;
    end
  end
  
  t_temp = table(label,value,levels,description);
  if h==1
    t0 = t_temp;
    t0(end+1,1) = {''};
    t0(end,2) = {''};
    t0(end,3) = {''};
    t0(end,4) = {''};
    t0.label{end+1} = 'partcipant_id';
    t0(end,2) = {''};
    t0(end,3) = {''};
    t0(end,4) = {'participant ID'};
    t0 = t0([end,1:end-1],:);
  elseif h==2
    t0 = [t0;t_temp];
    t0(end+1,1) = {''};
    t0(end,2) = {''};
    t0(end,3) = {''};
    t0(end,4) = {''};
  else
    t0 = [t0;t_temp];
  end
end


writetable(t0,'template_sfbquestionnaire.xlsx','Sheet','sfbquestionnaire');


fl = ter_listFiles(pwd,'*.pdf')
