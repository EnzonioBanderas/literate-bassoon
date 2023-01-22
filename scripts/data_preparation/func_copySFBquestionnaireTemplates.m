

fp0 = '/media/diskEvaluation/Evaluation/sfb1280a05study3/scripts/data_preparation';
fp_s = fullfile(fileparts(fileparts(fp0)),'sourcedata');
fn_template = fullfile(fp0,'template_sfbquestionnaire.xlsx');
fl = ter_listFiles(fp_s,'sub-*.pdf');
for i=1:numel(fl)
  [fp,pd] = fileparts(fileparts(fl{i}));
  if not(isequal(pd,'documents'))
    warning('folder structure unidentified for %s',fl{i});
    continue
  end
  pID = regexp(fl{i},'sub-[a-zA-Z0-9]*','match');
  pID = pID{end};
  fn_new1 = fullfile(fp,pd,[pID '_sfbquestionnaire.xlsx']);
  if exist(fn_new1,'file')==2
    fprintf('skipping existing finished xls sheet: %s\n',fn_new1)
    continue
  end
  fn_new2 = fullfile(fp,pd,[pID '_sfbquestionnaire_empty.xlsx']);
  if exist(fn_new2,'file')==2
    fprintf('skipping existing empty xls sheet: %s\n',fn_new2)
    continue
  end
  copyfile(fn_template,fn_new2);
end



fp0 = '/media/diskEvaluation/Evaluation/sfb1280a05study3/scripts/data_preparation';
fp_s = fullfile(fileparts(fileparts(fp0)),'sourcedata');
fn_template = fullfile(fp0,'template_study3questionnaire.xlsx');
fl = ter_listFiles(fp_s,'sub-*_study3questionnaire.pdf');
for i=1:numel(fl)
  [fp,pd] = fileparts(fileparts(fl{i}));
  if not(isequal(pd,'documents'))
    warning('folder structure unidentified for %s',fl{i});
    continue
  end
  pID = regexp(fl{i},'sub-[a-zA-Z0-9]*','match');
  pID = pID{end};
  fn_new1 = fullfile(fp,pd,[pID '_study3questionnaire.xlsx']);
  if exist(fn_new1,'file')==2
    fprintf('skipping existing finished xls sheet: %s\n',fn_new1)
    continue
  end
  fn_new2 = fullfile(fp,pd,[pID '_study3questionnaire_empty.xlsx']);
  if exist(fn_new2,'file')==2
    fprintf('skipping existing empty xls sheet: %s\n',fn_new2)
    continue
  end
  copyfile(fn_template,fn_new2);
end
