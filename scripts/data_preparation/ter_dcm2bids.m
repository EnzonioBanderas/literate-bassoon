function ter_dcm2bids(fp_dicomfolder,fp_bids,additionalInfo)
  
%   fn_dcm2niix = fullfile(fileparts(fileparts(mfilename('fullpath'))),...
%     'dcm2niix');
  fn_dcm2niix = evalc('!which dcm2niix');
  fn_dcm2niix = fn_dcm2niix(1:end-1);
  if exist(fn_dcm2niix,'file')~=2
    fn_dcm2niix = fullfile(fn_dcm2niix,'dcm2niix');
    if exist(fn_dcm2niix,'file')~=2
      error('dcm2niix not available in folder: %s',...
        fileparts(mfilename('fullpath'))); 
    end
  end
  
  fprintf('using dcm2niix: %s\n',fn_dcm2niix);
  
  pID = regexp(fp_dicomfolder,'sub-[a-zA-Z0-9]*','match');
  pID = pID{end};
  dir_temp = fullfile(fp_bids,pID,'temp');
  
  dtypes = {'anat' 'dwi' 'fmap' 'func'};
  stypes = {'gre','bold','T1w','dwi','sbref','epi'};
  addInfoLabels = {'task','acq','dir','inv','part','run'}; % run as last is important
  
  myCmd = sprintf(['%s -9 '...
    '-a n '...
    '-b y '...
    '-d 5 ' ...
    '-e n '...
    '-f %s_%%d ' ...
    '-g n '...
    '-i n '...
    '-l n '...
    '-m 2 '...
    '-o %s '...
    '-p y '...
    '-r n '...
    '-s n '...
    '-t n '...
    '-v 0 '...
    '-w 2 '...
    '-x n '...
    '-z y '...
    '%s'],...
    strrep(fn_dcm2niix,' ' ,'\ '),...
    pID, strrep(dir_temp,' ' ,'\ '),strrep(fp_dicomfolder,' ','\ '));
  if exist(dir_temp,'dir')~=7
    mkdir(dir_temp)
    system(myCmd);
  else
   dl =  struct2table(dir(dir_temp));
   dl = dl(~ismember(dl.name,{'.','..'}),:);
   if isempty(dl)
     system(myCmd);
   else
     warning(['temp folder not empty, working with content instead '...
       'of converting\n  %s'],dir_temp);
   end
  end
  fl =struct2table(dir(dir_temp));
  fl = fl(~fl.isdir,:);
  fl = cellfun(@fullfile,fl.folder,fl.name,'uni',0);
  cellfun(@delete,fl(contains(fl,{'AAhead','localizer','_MPR_'})));
  fl = fl(~contains(fl,{'AAhead','localizer','_MPR_'}));
  
  t1 = table(fl);
  t1.Properties.VariableNames{1} = 'fn_old';
  for i=1:numel(fl)
    fn_old = fl{i};
    [~,fn,fe] = fileparts(fn_old);
    if strcmpi(fe,'.gz')
      [~,fn,fe_tmp] = fileparts(fn);
      fe = [fe_tmp,fe]; %#ok<AGROW>
    end
    if contains(fn,'mp2rage') && contains(fn,'UNI-DEN')
      strrep(fn,'mp2rage','mp2rageUNIDEN');
    elseif contains(fn,'mp2rage') && contains(fn,'UNI')
      strrep(fn,'mp2rage','mp2rageUNI');
    end
      
    ses = regexp(fn_old,'ses-[a-zA-Z0-9]*','match');
    if isempty(ses) && nargin>2
      ses = regexp(additionalInfo,'ses-[a-zA-Z0-9]*','match');
    end
    if ~isempty(ses)
      ses = ses{end};
    else
      ses = '';
    end
    
    inv = regexp(fn_old,'INV[0-9]*','match');
    if isempty(inv) && nargin>2
      inv = regexp(additionalInfo,'inv-[a-zA-Z0-9]*','match');
    end
    if ~isempty(inv)
      inv = inv{end};
      inv = strrep(inv,'INV','inv-');
    else
      inv = '';
    end
    
    %2do for mp2rage with multiple phase images adopt part properly depending on
    %phase, also T1maps
    if ~isempty(inv)
      part = 'part-mag';
    else
      part = '';
    end
   
    %labels0 = strsplit(fn,'_');
    ind_dtype = cellfun(@(x) contains(lower(fn),x),dtypes);
    if sum(ind_dtype)==2 && ismember('fmap',dtypes(ind_dtype))
      ind_dtype = ismember(dtypes,'fmap'); 
    elseif sum(ind_dtype)~=1 
      warning(['type of volume data (anat, func, ...) '...
        'could not be identified:\n %s'],fn);
      continue
    end
    dtype = dtypes{ind_dtype};
    
    ind_stype = cellfun(@(x) contains(lower(fn),x),lower(stypes));
    if sum(ind_stype)~=1
      warning(['sequence type of volume data (bold, T1w, ...) '...
        'could not be identified:\n %s'],fn);
      continue
    end
    stype = stypes{ind_stype};
    useFor = '';
    if ~isempty(inv) && contains(fn_old, 'mp2rage')
      stype = 'MPRAGE';
    elseif strcmpi(dtype,'fmap') && ismember(stype,{'bold','dwi'})
      useFor = stype;
      stype = 'epi';
    end

    t1.dtype{i}  = dtype;
    t1.stype{i}  = stype;
    t1.ses{i}    = ses;
    t1.useFor{i} = useFor;
    %t1.inv{i}   = inv;
    %t1.part{i}  = part;
    
    addInfo = ['_' inv '_' part];
    for j=1:numel(addInfoLabels) 
      t1.(addInfoLabels{j}){i} = '';
      myExpr = [addInfoLabels{j} '-[a-zA-Z0-9]*'];
      tmp = regexp(fn,myExpr,'match');
      if ~isempty(tmp)
        t1.(addInfoLabels{j}){i} = strrep(tmp{end},[addInfoLabels{j} '-'],'');
        if ~(strcmpi(dtype,'fmap') && ...
            ismember(addInfoLabels(j),{'task'}))
          addInfo= [addInfo '_' tmp{end}]; %#ok<AGROW>
        end
      end
    end
    while contains(addInfo,'__')
      addInfo = strrep(addInfo,'__','_');
    end
    
    t1.addInfo{i} = addInfo;
    t1.fn{i}    = fn;
    t1.fe{i}    = fe;
        
  end
  %fp_out = fullfile(fp_bids,pID,ses,dtypes{ind_dtype});
  %disp(t1);

  ind_json = ismember(t1.fe,'.json');
  ind_nii  = ismember(t1.fe,{'.nii','.nii.gz'});
  ses = unique(t1.ses);
  
  % identify classic fieldmap files and award proper labels 
  % (magnitude 1 2 and phasediff), and update phasediff json information;
  for i=1:numel(ses)
    ind_ses = ismember(t1.ses,ses{i});
    
    ind_fmap = ismember(t1.dtype,'fmap');
    ind_gre = ismember(t1.stype,'gre');
    infos = unique(t1.addInfo(ind_ses&ind_fmap&ind_gre,:));
    for j=1:numel(infos)
      ind = ismember(t1.addInfo,infos(j))&ind_ses&ind_fmap&ind_gre;
      if sum(ind)~=6
        error('mismatch of file names')
      end
      for k=find(ind)'
        if strcmpi(t1.fn{k}(end-1:end),'ph')
          t1.stype{k} = 'phasediff';
        else
          t1.stype{k} = 'magnitude';
        end
      end
      ind_j_ph = find(ismember(t1.stype,'phasediff')&ind&ind_json);
      ind_j_m1 = find(ismember(t1.stype,'magnitude')&ind&ind_json);
      ind_j_m2 = ind_j_m1(2);
      ind_j_m1 = ind_j_m1(1);
      jinfo_ph = jsondecode(fileread(t1.fn_old{ind_j_ph}));
      jinfo_m1 = jsondecode(fileread(t1.fn_old{ind_j_m1}));
      jinfo_m2 = jsondecode(fileread(t1.fn_old{ind_j_m2}));
      if jinfo_m1.EchoTime<jinfo_m2.EchoTime
        t1.stype(ismember(t1.fn,t1.fn(ind_j_m1)))={'magnitude1'};
        t1.stype(ismember(t1.fn,t1.fn(ind_j_m2)))={'magnitude2'};
        jinfo_ph.EchoTime1 = jinfo_m1.EchoTime;
        jinfo_ph.EchoTime2 = jinfo_m2.EchoTime;
      else
        t1.stype(ismember(t1.fn,t1.fn(ind_j_m1)))={'magnitude2'};
        t1.stype(ismember(t1.fn,t1.fn(ind_j_m2)))={'magnitude1'};
        jinfo_ph.EchoTime1 = jinfo_m2.EchoTime;
        jinfo_ph.EchoTime2 = jinfo_m1.EchoTime;
      end
      ter_savePrettyJson(t1.fn_old{ind_j_ph},jinfo_ph);      
    end
  end
  
  
  t1.fn_new = cellfun(@(x,y,z,A,B) ...
    strrep(fullfile(y,x,[pID '_' y z '_' A B]),'__','_'),...
    t1.dtype,t1.ses,t1.addInfo,t1.stype,t1.fe,'uni',0);
  
  % dwi and bold fieldmaps might carry the same names, try to solve that
  t_fmap = t1(ismember(t1.dtype,'fmap'),:);
  for i=1:size(t_fmap,1)
    ind = find(ismember(t_fmap.fn_new ,t_fmap.fn_new(i)));
    if numel(ind)>1
      if numel(unique(t_fmap.useFor(ind))) == numel(ind)
        differences = t_fmap.useFor(ind);
      elseif numel(unique(t_fmap.task(ind))) == numel(ind)
        differences = t_fmap.task(ind);
      else
        differences = cellfun(@(x,y) [x y],...
          t_fmap.useFor(ind),t_fmap.task(ind));
      end
      
      for j=1:numel(ind)
        ind1 = ind(j);
        
        if ~isempty(differences{j})
          if isempty(t_fmap.acq{ind1})
            t_fmap.addInfo{ind1} = strrep(['_acq-' differences{j}...
              t_fmap.addInfo{ind1} ],'__','_');
          else
            t_fmap.addInfo{ind1} = strrep(t_fmap.addInfo{ind1}, ...
              ['_acq-' t_fmap.acq{ind1} ],...
              [ '_acq-' t_fmap.acq{ind1} differences{j}]);
          end
        end
      end
      t1(ismember(t1.dtype,'fmap'),:) = t_fmap;
      t1.fn_new = cellfun(@(x,y,z,A,B) ...
        strrep(fullfile(y,x,[pID '_' y z '_' A B]),'__','_'),...
        t1.dtype,t1.ses,t1.addInfo,t1.stype,t1.fe,'uni',0);
      t_fmap = t1(ismember(t1.dtype,'fmap'),:);
    end
  end
  
  
  
  % update fieldmap json info with IntendedFor fields
  for i=1:numel(ses)
    ind_ses = ismember(t1.ses,ses{i});
    ind_fm = ind_ses & ismember(t1.dtype,'fmap') & ind_json;
    for k=find(ind_fm)'
      if ismember(t1.useFor(k),'dwi') 
        ind_2go = ismember(t1.dtype,'dwi') & ind_ses & ind_nii;
      else
        ind_2go = ismember(t1.dtype,'func') & ind_ses & ind_nii & ...
          ismember(t1.task,t1.task(k));
      end
      if sum(ind_2go) > 0
        jinfo0 = jsondecode(fileread(t1.fn_old{k}));
        jinfo = jinfo0;
        jinfo.IntendedFor = t1.fn_new(ind_2go);
        fprintf('%s intended for \n',t1.fn_new{k})
        disp(jinfo.IntendedFor );
        if ~isequal(jinfo,jinfo0)
          ter_savePrettyJson(t1.fn_old{k},jinfo)
        end
      end
    end
  end
  
  for i=1:size(t1,1)
    t1.fn_new{i} = fullfile(fp_bids,pID,t1.fn_new{i});
  end
  for i=1:size(t1,1)
    if sum(ismember(t1.fn_new,t1.fn_new(i)))==1
      if ~isfolder(fileparts(t1.fn_new{i}))
        mkdir(fileparts(t1.fn_new{i}));
      end
      movefile(t1.fn_old{i},t1.fn_new{i});
    else
      warning('Non-unique new filename for : %s\n',t1.fn_old{i})
    end
  end
  
  fl = struct2table(dir(dir_temp));
  if size(fl,1)==2
    rmdir(dir_temp)
  elseif size(fl,1)==3 && strcmpi(fl.name{3},'niitable.csv')
    rmdir(dir_temp,'s')
  else
    writetable(t1,fullfile(dir_temp,'niitable.csv'));
  end
end

