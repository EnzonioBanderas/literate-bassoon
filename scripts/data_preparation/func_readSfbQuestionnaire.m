function func_readSfbQuestionnaire(fp_s,fp_d)
  fn_ptab = fullfile(fp_d,'participants.tsv');
  %ptab = readtable(fn_ptab,'filetype','text','delim','tab','treat','n/a');
  ptab = ter_readptab(fn_ptab)
  fp_pt = fullfile(fp_d,'phenotype');
  
  fp0 = fileparts(mfilename('fullpath'));
  
  fn_jdass = fullfile(fp_pt,'dass21g.json');
  if not(exist(fn_jdass,'file')==2)
    copyfile(fullfile(fp0,'dass21g.json'),fn_jdass);
  end
  ji_dass = jsondecode(fileread(fn_jdass));
  f_dass = fieldnames(ji_dass);
  
  fn_jehi = fullfile(fp_pt,'ehi.json');
  if not(exist(fn_jehi,'file')==2)
    copyfile(fullfile(fp0,'ehi.json'),fn_jehi);
  end
  ji_ehi = jsondecode(fileread(fn_jehi));
  f_ehi = fieldnames(ji_ehi);
  
  fn_jdemo = fullfile(fp_pt,'demographics.json');
  if not(exist(fn_jdemo,'file')==2)
    copyfile(fullfile(fp0,'demographics.json'),fn_jdemo);
  end
  ji_demo = jsondecode(fileread(fn_jdemo));
  f_demo = fieldnames(ji_demo);
  
  
  fn_t_ehi = fullfile(fp_pt,'ehi.tsv');
  if exist(fn_t_ehi,'file')==2
    t_ehi = readtable(fn_t_ehi,'filetype','text','delim','tab','treat','n/a');
  else
    t_ehi = table('size',[0,numel(f_ehi)],...
      'variablename',['participant_id',f_ehi(2:end)'],...
      'variabletype',['cell','double','double','double',...
      repmat({'cell'},1,numel(f_ehi)-4)]);
  end
  t_ehi0 = t_ehi;
  
  fn_t_dass = fullfile(fp_pt,'dass21g.tsv');
  if exist(fn_t_dass,'file')==2
    t_dass = readtable(fn_t_dass,'filetype','text','delim','tab','treat','n/a');
  else
    t_dass = table('size',[0,numel(f_dass)],...
      'variablename',['participant_id',f_dass(2:end)'],...
      'variabletype',['cell','cell','cell','cell',...
      repmat({'double'},1,numel(f_dass)-4)]);
  end
  t_dass0 = t_dass;
  
  
  fn_t_demo = fullfile(fp_pt,'demographics.tsv');
  f_string = {
    'participant_id'
    'sex'
    'fem_hormcontr_med'       
    'fem_hormcontr_last'      
    'fem_hormcontr_next'      
    'fem_menstr_last'         
    'fem_menstr_next'         
    'sense_eyesight_deficit'  
    'health_meds_descr'       
    'health_psych_descr'      
    'health_neuro_descr'      
    'handedness'              
    'anc_mat_grandmother_nat' 
    'anc_mat_grandmother_reg' 
    'anc_mat_grandfather_nat' 
    'anc_mat_grandfather_reg' 
    'anc_pat_grandmother_nat' 
    'anc_pat_grandmother_reg' 
    'anc_pat_grandfather_nat' 
    'anc_pat_grandfather_reg' 
    'anc_mother_nat'          
    'anc_mother_reg'          
    'anc_father_nat'          
    'anc_father_reg'};
  if exist(fn_t_demo,'file')==2
    t_demo = readtable(fn_t_demo,'filetype','text','delim','tab','treat','n/a');
    for i=1:numel(f_string)
      if isnumeric(t_demo.(f_string{i}))
        t_demo.(f_string{i}) = num2cell(t_demo.(f_string{i}));
      end
    end
  else
    vn = ['participant_id',f_demo(2:end)'];
    vt = repmat({'cell'},1,numel(f_demo));
    vt(not(ismember(vn,f_string))) = {'double'};
    t_demo = table('size',[0,numel(f_demo)],...
      'variablename',vn,'variabletype',vt);
  end
  t_demo0 = t_demo;
  
  for i=1:size(ptab,1)
    pID = ptab.participant_id{i};
    fn_xls = fullfile(fp_s,pID,'documents',[pID '_sfbquestionnaire.xlsx']);
    if not(exist(fn_xls,'file') == 2)
      fprintf('not found: \n %s \n',fn_xls);
      continue
    end
    t1 = readtable(fn_xls);
    if not(ismember('label',fieldnames(t1)))
      fprintf('not proper format: \n %s \n',fn_xls);
      continue
    end
    
    ind = find(ismember(t_ehi.participant_id,pID));
    if isempty(ind)
      t_ehi.participant_id(end+1) = {pID};
      %t_ehi(end,2:end) = num2cell(nan(1,size(t_ehi,2)-1));
      ind = size(t_ehi,1);
    end
    for j=1:numel(f_ehi)
      ind1 = ismember(t1.label,f_ehi{j});
      if sum(ind1)==1
        t_ehi.(f_ehi{j})(ind) = t1.value(ind1);
      end
    end
    
    ind = find(ismember(t_dass.participant_id,pID));
    if isempty(ind)
      t_dass.participant_id(end+1) = {pID};
      ind = size(t_dass,1);
    end
    for j=1:numel(f_dass)
      ind1 = ismember(t1.label,f_dass{j});
      if sum(ind1)==1
        try
          t_dass.(f_dass{j})(ind) = str2double(t1.value{ind1});
        catch
          t_dass.(f_dass{j})(ind) = nan;
        end
      end
    end
    
    ind = find(ismember(t_demo.participant_id,pID));
    if isempty(ind)
      t_demo.participant_id(end+1) = {pID};
      ind = size(t_demo,1);
    end
    for j=1:numel(f_demo)
      ind1 = ismember(t1.label,f_demo{j});
      if sum(ind1)==1
        if iscell(t_demo.(f_demo{j})(ind)) 
          t_demo.(f_demo{j})(ind) = t1.value(ind1);
        else
          try 
            t_demo.(f_demo{j})(ind) = str2double(t1.value{ind1});
          catch
            t_demo.(f_dass{j})(ind) = nan;
          end
        end
      end
    end
    if not(isequal(t_demo.participant_id{ind},pID))
      error('participant_id mismatch : %s',fn_xls);
    end
    
    
    
  end

  
  %% calculate ehi scores and laterality coefficient
  for i=1:size(t_ehi,1)
    ind = find(contains(f_ehi,'ehi_q'));
    qval = '';
    for j=1:numel(ind)
      qval = [qval,t_ehi.(f_ehi{ind(j)}){i}]; %#ok<AGROW>
    end
    scoreR = sum(ismember(qval,'r')) + sum(ismember(qval,'n'));
    scoreL = sum(ismember(qval,'l')) + sum(ismember(qval,'n'));
    ehiscore = round(100*(scoreR-scoreL)/(scoreR+scoreL)); 
    t_ehi.ehi_lc(i)     = ehiscore;
    t_ehi.ehi_scoreR(i) = scoreR;
    t_ehi.ehi_scoreL(i) = scoreL;
  end
  t_ehi = sortrows(t_ehi);
  if not(isequaln(t_ehi,t_ehi0))
    fprintf('writing updated ehi table\n')
    writetable(t_ehi,fn_t_ehi,'filetype','text','delim','tab');
  end
  
  
  %% calculate dass scores 
  ind_d = not(cellfun(@isempty,regexp(f_dass,'_d$')));
  ind_a = not(cellfun(@isempty,regexp(f_dass,'_a$')));
  ind_s = not(cellfun(@isempty,regexp(f_dass,'_s$')));
  t_dass.dass21score_depression = sum(table2array(t_dass(:,ind_d)),2)*2;
  t_dass.dass21score_anxiety    = sum(table2array(t_dass(:,ind_a)),2)*2;
  t_dass.dass21score_stress     = sum(table2array(t_dass(:,ind_s)),2)*2;
  l_d = struct2cell(ji_dass.dass21score_depression.Levels);
  l_a = struct2cell(ji_dass.dass21score_anxiety.Levels);
  l_s = struct2cell(ji_dass.dass21score_stress.Levels);
  for i=1:size(t_dass,1)
    n_d = t_dass.dass21score_depression(i)+1;
    n_a = t_dass.dass21score_anxiety(i)+1;
    n_s = t_dass.dass21score_stress(i)+1;
    t_dass.dass21label_depression(i) = l_d(n_d);
    t_dass.dass21label_anxiety(i) = l_a(n_a);
    t_dass.dass21label_stress(i) = l_s(n_s);
  end
  t_dass = sortrows(t_dass);
  if not(isequaln(t_dass,t_dass0))
    fprintf('writing updated dass table\n')
    writetable(t_dass,fn_t_dass,'filetype','text','delim','tab');
  end
  
  
  %% demographics
  % force lower case letters
  t_demo.sex = lower(t_demo.sex);
  t_demo.handedness = lower(t_demo.handedness);
  
  % total education years
  ind_edu = find(contains(t_demo.Properties.VariableNames,'edu_y_') & ...
    not(contains(t_demo.Properties.VariableNames,'edu_y_total')));
  for i=1:numel(ind_edu)
    if iscell(t_demo.(t_demo.Properties.VariableNames{ind_edu(i)}))
      t_demo.(t_demo.Properties.VariableNames{ind_edu(i)}) = ...
        str2double(t_demo.(t_demo.Properties.VariableNames{ind_edu(i)}));
    end
  end
  if iscell(t_demo.edu_y_total)
    t_demo.edu_y_total = str2double(t_demo.edu_y_total);
  end
  
  edu_total = nansum(table2array(t_demo(:,ind_edu)),2);
  for i=1:size(t_demo,1)
    if isnan(t_demo.edu_y_total(i))
      t_demo.edu_y_total(i) = edu_total(i);
    elseif not(t_demo.edu_y_total(i) == edu_total(i))
      warning('mismatch in total education time: %s',...
        t_demo.participant_id{i});
    else
      t_demo.edu_y_total(i) = edu_total(i);
    end
  end
  
  f_anc = f_demo(contains(f_demo,'anc_'));
  str2rep = {
    'Albanien'          'Albanian'  
    'Albanisch'         'Albanian' 
    'Australien'        'Australian' 
    'Bosnisch'          'Bosnian'
    'Chinesisch'        'Chinese'
    'Deustch'           'German'         
    'Deutsch'           'German'
    'Deutschland'       'German'
    'Französisch'       'French'
    'Georgisch'         'Georgian'
    'Indien'            'Indian'
    'Indonesien'        'Indonesian'
    'Indonesisch'       'Indonesian'
    'Iran'              'Iranian'
    'Iranisch'          'Iranian'
    'Italienisch'       'Italian'
    'Kasachstan'        'Kazakhstanian'
    'Marokanisch'       'Moroccan'
    'Mazedonien'        'Macedonian'
    'Niederländisch'    'Dutch'
    'Pakistan'          'Pakistani'
    'Polen'             'Polish'
    'Polnisch'          'Polish'
    'Polnisch/Deutsch'  'Polish'
    'Rumänisch'         'Romanian'
    'Russisch'          'Russian'
    'Russland'          'Russian'
    'Spanien'           'Spanish'
    'Syrien'            'Syrian'
    'Syrisch'           'Syrian'
    'Südkorea'          'South Korean'
    'Tschetchenisch'    'Chechen'
    'Türkisch'          'Turkish'
    'Ukrainisch'        'Ukrainian'
    'marokanisch'       'Moroccan'
    'deutsch'           'German'
    'Ägyptisch'         'Egyptian'
    '?'                 'n/a'
    'Kurdisch'          'Kurdish'
    'kurdisch'          'Kurdish'
    };
  for i=1:numel(f_anc)
    for j=1:size(str2rep,1)
      t_demo.(f_anc{i}) = strrep(t_demo.(f_anc{i}),str2rep{j,:});
    end
  end
  
  for i=1:numel(f_string)
    for j=1:size(t_demo,1)
      if isempty(t_demo.(f_string{i}){j})
        t_demo.(f_string{i}){j} = 'n/a';
      end
    end
  end
  
  
  t_demo = sortrows(t_demo);
  if not(isequaln(t_demo,t_demo0))
    % some postprocessing is need to comply with BIDS, therefore first work
    % on a temporary file before schecking wether to overwrite the old one
    fn_temp = fullfile(fp0,'dempgraphics.temp.tsv');
    writetable(t_demo,fn_temp,'filetype','text','delim','tab');
    txt = fileread(fn_temp);
    txt2 = strrep(txt,'NaN','n/a');
    if not(isequal(txt,txt2))
      fid = fopen(fn_temp,'w');
      fprintf(fid,txt2);
      fclose(fid);
    end
    try
      txt0 = fileread(fn_t_demo);
    catch
      txt0 = '';
    end
    txt  = fileread(fn_temp);
    if not(isequaln(txt,txt0))
      fprintf('writing updated demographics table\n')
      copyfile(fn_temp,fn_t_demo);
    end
    delete(fn_temp);

  end
  
  
  
end
