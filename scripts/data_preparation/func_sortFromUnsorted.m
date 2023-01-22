function func_sortFromUnsorted(fp_us,fp_s,fp_dis,varargin)

if ~ismember('silent',lower(varargin))
  fprintf('\n%s\n%s -Sorting unsorted sourcefile folder\n',...
    repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'));
end
if ~isfolder(fp_us)
  fprintf('no unsorted sourcedata folder found\n')
  return
end
tic

[~,tmp] = fileparts(fp_s);
fp_diss = fullfile(fp_dis,tmp);
copyFiles = true;
if ismember('move',lower(varargin))
  copyFiles = false;
end

if ~ismember({'skipphysio'},lower(varargin))
  if ~contains('silent',lower(varargin))
    fprintf('Sorting physio files\n');
  end
  if nargin < 1
    fp_us='/media/diskEvaluation/Evaluation/sfb1280a05study3/sourcefiles_unsorted';
    fp_s = fullfile(fileparts(fp_us),'sourcedata');
  end
  fl = ter_listFiles(fp_us,'*.acq');
  fl = fl(~contains(fl,'.mat.acq'));
  fl = fl(~contains(lower(fl),'test'));

  str2replace = {
    'Sub'               'sub'
    'sub_'              'sub-'
    'Ses'               'ses'
    'ses_'              'ses-'
    '-ses'              '_ses'
    'ses1'              'ses-1'
    'ses2'              'ses-2'
    'Quest'             'quest'
    'questionaire'      'questionnaire'
    'questionnaires'    'questionnaire'
    'Checklist'         'checklist'
    'CheckList'         'checklist'
    'Sfbquest'          'sfbquest'
    'SFbquest'          'sfbquest'
    'SFBquest'          'sfbquest'
    'SFBQuest'          'sfbquest'
    };

  for i=1:numel(fl)
    [~,fn0,fe0] = fileparts(fl{i});
    fn0 = [fn0 fe0]; %#ok<AGROW>
    for j=1:size(str2replace,1)
      fn0 = strrep(fn0,str2replace{j,:});
    end
    pID = regexp(fn0,'sub-[a-zA-Z0-9]*','match');
    if isempty(pID)
      warning('Unidentified participant ID for file: %s',fl{i});
      continue
    end
    pID = pID{end};
    sID = regexp(fn0,'ses-[a-zA-Z0-9]*','match');
    if isempty(sID)
      warning('Unidentified session ID for file: %s',fl{i});
      continue
    end
    sID = sID{end};
    fp_new = fullfile(fp_s,pID,sID,'physio');
    fn_new = fullfile(fp_new,fn0);
    fn_new_dis = fullfile(fp_diss,pID,sID,'physio',fn0);
    if exist(fn_new,'file')~=2 && exist(fn_new_dis,'file')~=2
      % sort if not yet present
      if ~isfolder(fp_new)
        mkdir(fp_new);
      end
      if copyFiles
        copyfile(fl{i},fn_new);
      else
        movefile(fl{i},fn_new);
      end
      disp(fn_new);
    elseif ~copyFiles && ~isequal(fl{i},fn_new)
      delete(fl{i})
    end
    
    fn_mat = strrep(fl{i},'.acq','.mat');
    if exist(fn_mat,'file')~=2
      warning('Missing mat file for %s',fl{i})
    else
      fn_newmat     = strrep(fn_new,    '.acq','.mat');
      fn_newmat_dis = strrep(fn_new_dis,'.acq','.mat');
      if exist(fn_newmat,'file')~=2 && exist(fn_newmat_dis,'file')~=2
        %movefile(fn_mat, fn_newmat);
        if copyFiles
          copyfile(fn_mat, fn_newmat);
        else
          movefile(fn_mat, fn_newmat);
        end
        disp(fn_newmat);
      elseif ~copyFiles && ~isequal(fn_mat,fn_newmat)
        delete(fn_mat)
      end
    end
  end
end

%% sort log files 
if ~ismember('skiplogs',lower(varargin))
  if ~ismember('silent',lower(varargin))
    fprintf('Sorting log files\n');
  end
  fl = ter_listFiles(fp_us,'sub-*.txt');
  for i=1:numel(fl)
    [~,fn,fe] = fileparts(fl{i});
    fID = ter_parseFname(fn);
    if ~isfield(fID,'pID')
      continue
    end
    pID = fID.pID;
    if ~isfield(fID,'ses')
      sID = '';
    else
      sID = fID.ses;
    end
    fn_old = fl{i};
    fn_new = fullfile(fp_s,pID,sID,'stimlogs',[fn,fe]);
    fn_new_dis = fullfile(fp_diss,pID,sID,'stimlogs',[fn,fe]);
    if ~isequal(fl{i},fn_new) && ~isequal(fl{i},fn_new_dis)
      if exist(fn_new,'file')~=2 && exist(fn_new_dis,'file')~=2
        if ~isfolder(fileparts(fn_new))
          mkdir(fileparts(fn_new))
        end
        if copyFiles
          copyfile(fn_old,fn_new)
        else
          movefile(fn_old,fn_new)
        end
        disp(fn_new);
      elseif exist(fn_new,'file')==2
        warning('file already exists : %s',fn_new);
      elseif exist(fn_new_dis,'file')==2
        warning('file already exists : %s',fn_new_dis);
      end
    end
  end
end

%% sort log files 
if ~ismember('skiplogs',lower(varargin))
  if ~ismember('silent',lower(varargin))
    fprintf('Sorting eyetracker files\n');
  end
  fl = ter_listFiles(fp_us,'sub-*ses-*run-*_eyetrack*');
  for i=1:numel(fl)
    [~,fn,fe] = fileparts(fl{i});
    fID = ter_parseFname(fn);
    if ~isfield(fID,'pID')
      continue
    end
    pID = fID.pID;
    if ~isfield(fID,'ses')
      sID = '';
    else
      sID = fID.ses;
    end
    fn_old = fl{i};
    fn_new = fullfile(fp_s,pID,sID,'eyetracking',[fn,fe]);
    fn_new_dis = fullfile(fp_diss,pID,sID,'eyetracking',[fn,fe]);
    if ~isequal(fl{i},fn_new) && ~isequal(fl{i},fn_new_dis)
      if exist(fn_new,'file')~=2 && exist(fn_new_dis,'file')~=2
        if ~isfolder(fileparts(fn_new))
          mkdir(fileparts(fn_new))
        end
        if copyFiles
          copyfile(fn_old,fn_new)
        else
          movefile(fn_old,fn_new)
        end
        disp(fn_new);
      elseif exist(fn_new,'file')==2
        warning('file already exists : %s',fn_new);
      elseif exist(fn_new_dis,'file')==2
        warning('file already exists : %s',fn_new_dis);
      end
    end
  end
end



%% Sort Dicoms
if ~ismember({'skipdicoms'},lower(varargin))
  fprintf('Sorting dicom files\n');
  %warning('off','images:dicominfo:fileVRDoesNotMatchDictionary')
  %ter_sortDicom('presorted',false,'structure','subjectgivenname','inputdir',fp0)

  filemap = table();
  filemap.pID={};
  filemap.studyID={};
  filemap.sessionID={};
  filemap.fn_old={};
  filemap.fn_new={};
  filemap.ind_dcmdir = zeros(0,1);
  filemap.ind_pat    = zeros(0,1);
  filemap.ind_study  = zeros(0,1);
  %fl = ter_listAllExtension(fp_us,'',6,0);
  fl_dcmdir = ter_listFiles(fp_us,'DICOMDIR');
  dcmdir = cell(size(fl_dcmdir));
  for i=1:numel(fl_dcmdir)
    dcmdir{i} = images.dicom.parseDICOMDIR(fl_dcmdir{i}); 
    dcmdir{i}.path = fileparts(fl_dcmdir{i});
  end
  warning('off','MATLAB:table:RowsAddedExistingVars');

  for i=1:numel(dcmdir)
    for j=1:numel(dcmdir{i}.Patients)
      tmp = strsplit(dcmdir{i}.Patients(j).Payload.PatientName,'^');
      pID = tmp{end};
      if ~contains(pID,'sub-')
        pID = ['sub-' upper(pID)];
      end
      for k=1:numel(dcmdir{i}.Patients(j).Studies)
        studyID = dcmdir{i}.Patients(j).Studies(k).Payload.StudyID;
        for l=1:numel(dcmdir{i}.Patients(j).Studies(k).Series)
          serinfo = dcmdir{i}.Patients(j).Studies(k).Series(l).Payload;
          fp_new = fullfile(fp_s,pID,['study-' studyID],'DICOM',...
            sprintf('%s_ser-%03.0f_%s',pID,serinfo.SeriesNumber,...
            serinfo.SeriesDescription));
          for m=1:numel(dcmdir{i}.Patients(j).Studies(k).Series(l).Images)
            fn_old = fullfile(dcmdir{i}.path,strrep(...
              dcmdir{i}.Patients(j).Studies(k).Series(l).Images(m).Payload.ReferencedFileID,...
              '\',filesep));
            [~,fn,fe] = fileparts(fn_old);
            fn_new    = fullfile(fp_new,[fn fe]);
            filemap.fn_old{end+1} = fn_old;
            filemap.fn_new{end}   = fn_new;
            filemap.studyID{end} = studyID;
            filemap.pID{end} = pID;
            filemap.ind_dcmdir(end) = i;
            filemap.ind_pat(end)    = j;
            filemap.ind_study(end)  = k;
          end
          for m=1:numel(dcmdir{i}.Patients(j).Studies(k).Series(l).SRDocuments)
            fn_old = fullfile(dcmdir{i}.path,strrep(...
              dcmdir{i}.Patients(j).Studies(k).Series(l).SRDocuments(m).Payload.ReferencedFileID,...
              '\',filesep));
            [~,fn,fe] = fileparts(fn_old);
            fn_new    = fullfile(fp_new,[fn fe]);
            filemap.fn_old{end+1} = fn_old;
            filemap.fn_new{end}   = fn_new;
            filemap.studyID{end} = studyID;
            filemap.pID{end} = pID;
            filemap.ind_dcmdir(end) = i;
            filemap.ind_pat(end)    = j;
            filemap.ind_study(end)  = k;
          end
        end
      end
    end
  end
  toc
  participants = unique(filemap.pID);
  for i=1:numel(participants)
    pID = participants{i};
    studies = unique(filemap.studyID(ismember(filemap.pID,pID)));
    if numel(studies)==2
      for j=1:numel(studies)
        ind = ismember(filemap.pID,pID)&ismember(filemap.studyID,studies(j));
        sesID = sprintf('ses-%d',j);
        filemap.sessionID(ind) = {sesID};
        filemap.fn_new(ind) = strrep(filemap.fn_new(ind),...
          ['study-' studies{j}],sesID);
        filemap.fn_new_dis(ind) = strrep(filemap.fn_new(ind),fp_s,fp_diss);
      end
    else
      warning('incomplete number of studies for subject %s', pID);
      filemap = filemap(not(ismember(filemap.pID,pID)),:);
    end
  end
  for i=1:size(filemap,1)
    if exist(filemap.fn_old{i},'file')==2
      if exist(filemap.fn_new{i},'file')~=2 && ...
          exist(filemap.fn_new_dis{i},'file')~=2
        if ~isfolder(fileparts(filemap.fn_new{i}))
          mkdir(fileparts(filemap.fn_new{i}));
        end
        if copyFiles
          copyfile(filemap.fn_old{i},filemap.fn_new{i}); 
        else
          movefile(filemap.fn_old{i},filemap.fn_new{i});
        end
        %disp(filemap.fn_new{i});
      end
    else
      %if exist(filemap.fn_new{i},'file')~=2 && ...
      %    exist(filemap.fn_new_dis{i},'file')~=2
      %  warning('file does not exist: %s',filemap.fn_old{i});
      %end
    end
  end
  inds_dd = unique(filemap.ind_dcmdir);
  for i=1:numel(inds_dd)
    ind1 = ismember(filemap.ind_dcmdir,inds_dd(i));
    inds_pat = unique(filemap.ind_pat(ind1));
    for j=1:numel(inds_pat)
      ind2 = ismember(filemap.ind_pat,inds_pat(j)) & ind1;
      inds_study = unique(filemap.ind_study(ind2));
      for k=1:numel(inds_study)
        ind3 = ismember(filemap.ind_study,inds_study(k)) & ind2;
        fp = fileparts(fileparts(filemap.fn_new{find(ind3,1,'first')}));
        try
          fp_dis = fileparts(fileparts(...
            filemap.fn_new_dis{find(ind3,1,'first')}));
        catch
          disp('häh?')  
        end
        clearvars dicomdir
        dicomdir.Patients.Payload = ...
          dcmdir{inds_dd(i)}.Patients(inds_pat(j)).Payload;
        dicomdir.Patients.Studies = ...
          dcmdir{inds_dd(i)}.Patients(inds_pat(j)).Studies(inds_study(k));
        for l=1:numel(dicomdir.Patients.Studies.Series)
          for m=1:numel(dicomdir.Patients.Studies.Series(l).Images)
            refID = dicomdir.Patients.Studies.Series(l).Images(m).Payload.ReferencedFileID;
            ind = contains(filemap.fn_old,strrep(refID,'\',filesep));
            refID_new = strrep(filemap.fn_new{ind},fp,'');
            if strcmp(refID_new(1),filesep)
              refID_new = refID_new(2:end);
            end
            dicomdir.Patients.Studies.Series(l).Images(m).Payload.ReferencedFileID = ...
              strrep(refID_new,filesep,'\');
          end
        end
        fn_new     = fullfile(fp,    'dicomdir.mat');
        fn_new_dis = fullfile(fp_dis,'dicomdir.mat');
        if exist(fn_new,'file')~=2 && exist(fn_new_dis,'file')~=2
          if ~isfolder(fp)
            mkdir(fp)
          end
          save(fn_new,'dicomdir');
        end
      end
    end
  end
end



  %% sort documents
  
  if ~ismember({'skipdocs'},lower(varargin))
    fl = ter_listFiles(fp_us,{'*.pdf','*.xlsx','.docx'});
    fl2 = fl;
    for j=1:size(str2replace,1)
      fl2 = strrep(fl2,str2replace{j,:});
    end
    for i=1:numel(fl)
      [~,fn,fe] = fileparts(fl2{i});
     
      % process randomlists
      if strcmpi(fn,'Randomisierungsliste')
        finfo = dir(fl2{i});
        fn_new = fullfile(fp_s,'documents',...
          ['randomlist_' datestr(finfo.datenum,'yyyy-mm-dd') fe]);
        if exist(fn_new,'file')
          delete(fl2{i})
        else
          movefile(fl2{i},fn_new);
        end
        continue
      end
     
      % process all other files
      pID = regexp(fl2{i},'sub-[0-9a-zA-Z]*','match');
      if isempty(pID)
        warning('Cannot identify participant ID for file:\n   %s\n',fl{i});
        continue
      else
        pID = pID{end};
      end
      
      fn_new = fullfile(fp_s,pID,'documents',[fn fe]);
      if exist(fn_new,'file')~=2
        if ~isfolder(fileparts(fn_new))
          mkdir(fileparts(fn_new));
        end
        if copyFiles
          copyfile(fl{i},fn_new); 
        else
          movefile(fl{i},fn_new);
        end
      elseif ~copyFiles && ~isequal(fl{i},fn_new)
        delete(fl{i});
      end
      
    end
    
  end


  %% delete empty folders
  ter_deleteEmptyFolders(fp_us);
  
  %% clear away unneeded DIOCMDIR files
  for j=1:2
    fl = ter_listFiles(fp_us,'DICOMDIR');
    for i=1:numel(fl)
      dl = struct2table(dir(fullfile(fileparts(fl{i}),'*')));
      dl = dl(~ismember(dl.name,{'..','.'}),:);
      if size(dl,1)==1
        fn_new = strrep(fl{i},fileparts(fl{i}),...
          fullfile(fp_dis,'DICOMDIRS'));
        if exist(fileparts(fn_new),'dir')~=7
          mkdir(fileparts(fn_new));
        end
        movefile(fl{i},fn_new);
      end
    end
    if numel(fl)>1
      ter_deleteEmptyFolders(fp_us);
    end
  end

end


