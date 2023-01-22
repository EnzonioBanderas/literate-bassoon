function func_convertDicom2Bids(fp_s,fp_d)


fn_ptab = fullfile(fp_d,'participants.tsv');
try
  ptab = readtable(fn_ptab,'FileType','text','delim','tab','treat','n/a');
  ptab0 = ptab;
catch
  warning('cloud not read participants.tsv');
end

dl = struct2table(dir(fullfile(fp_s,'sub-*')));
dirs2work = {}; 
for i=1:size(dl,1)
  fp1 = fullfile(dl.folder{i},dl.name{i},'ses-1','DICOM');
  fp12 = fullfile(fp_d,dl.name{i},'ses-1','func');
  fp2 = fullfile(dl.folder{i},dl.name{i},'ses-2','DICOM');
  fp22 = fullfile(fp_d,dl.name{i},'ses-2','func');
  if isfolder(fp1) 
    if ~isfolder(fp12)
      dirs2work = [dirs2work;fp1]; %#ok<AGROW>
    end
  end
  if isfolder(fp2) 
    if ~isfolder(fp22)
      dirs2work = [dirs2work;fp2]; %#ok<AGROW>
    end
  end
end
      
fprintf('These %d folders will be converted:\n',numel(dirs2work))
disp(dirs2work)

for i=1:numel(dirs2work)
  fprintf('Importing %s\n',dirs2work{i});
  sesID = regexp(dirs2work{i},'ses-[0-9]*','match');
  try
    ter_dcm2bids(dirs2work{i},fp_d,sesID{end}); 
  catch
    warning('error processing %s',dirs2work{i})
  end
end

% update ptab
try
  for i=1:numel(dirs2work)
    sesID = regexp(dirs2work{i},'ses-[0-9]*','match');
    pID   = regexp(dirs2work{i},'sub-[0-9a-zA-Z]*','match');
    pID = pID{end};
    label = ['date_' strrep(sesID{end},'-','_')];
    dl = struct2table(dir(fullfile(dirs2work{i})));
    dl = dl(dl.isdir & ~ismember(dl.name,{'.','..'}),:);
    fl = ter_listFiles(fullfile(dirs2work{i},dl.name{1}),'*');
    if ~ismember(pID,ptab.participant_id)
      ptab.participant_id{end+1} = pID;
      varnames = ptab.Properties.VariableNames;
      for j=1:numel(varnames)
        if iscell(ptab.(varnames{j})(end)) && ...
            isempty(ptab.(varnames{j}){end})
          ptab.(varnames{j})(end)={'n/a'};
        elseif isnumeric(ptab.(varnames{j})(end)) && ...
            isequal(ptab.(varnames{j})(end), 0)
          ptab.(varnames{j})(end)=nan;
        end
      end
    end
    ind = ismember(ptab.participant_id,pID);

    if numel(fl)>2
      dinfo = dicominfo(fl{end});
      tmp = datetime(dinfo.StudyDate,'format','yyyyMMdd');
      tmp.Format = 'yyyy-MM-dd';
      ptab.(label)(ind) = tmp;
      ptab.sex{ind} = lower(dinfo.PatientSex);
    end
    
  end
  try
    ptab = sortrows(ptab);
  catch
    fprintf('Cannot sort participants.tsv')
  end
  if ~isequaln(ptab,ptab0)
    writetable(ptab,fn_ptab,'filetype','text','delim','tab');
    txt = fileread(fn_ptab);
    txt2 = strrep(strrep(txt,'NaN','n/a'),'NaT','n/a');
    if ~isequal(txt,txt2)
      fid = fopen(fn_ptab,'w');
      fprintf(fid,txt2);
      fclose(fid);
    end
  end
catch
  fprintf('trouble updating participants.tsv')
end
  
  


