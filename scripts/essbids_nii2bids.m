function [t_converted,t_excluded] = essbids_nii2bids(fp_dicomfolder,fp_bids,additionalInfo,varargin)
%% ESSBIDS_DICOM2BIDS
%
% Converts all DICOM files in a folder and it's subfolders into .nii.gz 
% files using Chris Rordan's dcm2niix, and sorts them into a BIDS 
% compatible data structure. Sorting is based on BIDS labels in protocol
% names. Typical modifications to json fields, like IntendedFor for fmaps 
% are applied, and, if missing, task descriptor jsons are created.
% "participants.tsv" is appended with participant_id, sex and age gathered
% from DICOM fields. Optional input allows modification of BIDS fields.
% Intended to convert data of a single subject per run.
%
%  Inputs: 
%    fp_dicomfolder : folder of a single subject dicom data
%    fp_bids        : bids group level folder for output
%    additionalInfo : text string or structure with additional BIDS field
%                     information, will take precedence over file name
%                     level information
%
%  Optional inputs:
%    filenameMods : n x 2 cell array of regular expressions to
%                    find/replace in output names
%    ignoreSerNum : numerical array of series numbers to remove immediately
%
%  Outputs: One overview tables covering the converted files and a second 
%           covering the files where conversion failed 
%
%  Examples:
%    fp_dcm = './sourcedata/sub-001/DICOM'
%    dp_bids = './rawdata'
%    additionalinfo = 'ses-mri';
%    filenameMods = { ...
%      '_dti_','_dwi_';
%      '_qsm_','_anat_'}
%    %%
%   
%
%
%


  
  %% presets
  % todo :check on config file here
  types_suffix = {'gre','bold','T1w','T2w','dwi','sbref','epi','FLAIR',...
    'MP2RAGE','UNIT1','MEGRE'};
  files2ignore = {'AAhead','localizer','_MPR_','_TRACE','_TENSOR_B0',...
    '_FA.','_ColFA.','_ADC.','AAHead','tfl_b1map','_Aspire_','aspire','_T2star_', ...
    'dummy', 'REMOVE'};


  %% verify dcm2niix availability
  [a,m,v] = essbids_checkToolAvailability('dcm2niix'); %#ok<ASGLU>
  if not(a.dcm2niix)
    error('Please make dcm2niix available before running dicom2bids');
  end
  fprintf('%s\n',m.dcm2niix)
  try
    if ispc()
      fn_dcm2niix = strtrim(evalc('!where.exe dcm2niix'));
    else
      fn_dcm2niix = strtrim(evalc('!which dcm2niix'));
    end
    fprintf('Using %s\n',fn_dcm2niix)
  catch
    fn_dcm2niix = 'dcm2niix';
  end
  
  
  %% verify input
  if not(isfolder(fp_dicomfolder))
    error('essbids:dicom2bids:input','Dicom folder does not exist')
  end
  default_ignoreSerNum = [];
  default_filenameMods = {};
  p = inputParser;
  addParameter(p,'ignoresernum',default_ignoreSerNum,@isnumeric);
  addParameter(p,'filenamemods',default_filenameMods,@iscell);
  lowVarargin = varargin;
  ind_str = cellfun(@ischar,lowVarargin);
  lowVarargin(ind_str) = lower(lowVarargin(ind_str));
  parse(p,lowVarargin{:});
  ignoreSerNum = p.Results.ignoresernum(:);  
  filenameMods = p.Results.filenamemods;
  if mod(numel(filenameMods),2)~=0
    filenameMods = default_filenameMods;
    warning('filenameMods number uneven, ignoring mods');
  
  end
    
  
   %% gather bids flags from input
  bf_dir = essbids_parseLabel(fp_dicomfolder,'dropFinfo',true);
  if nargin < 3
    additionalInfo = '';
  end
  if nargin < 4
    filenameMods = {};
  end
  if isempty(additionalInfo)
    bf_add = struct();
  elseif ischar(additionalInfo)
    bf_add = essbids_parseLabel(additionalInfo,'dropFinfo',true);
  elseif isstruct(additionalInfo)
    bf_add = additionalInfo;
  else
    bf_add = struct();
  end
  fn_add = fieldnames(bf_add);
  bf0 = bf_dir;
  for i=1:numel(fn_add)
    bf0.(fn_add{i}) = bf_add.(fn_add{i});
  end
  bf0 = essbids_parseLabel(essbids_buildFileName(bf0,...
    'relativePath',true),'dropFinfo',true);
  % in minimum a participant id has to be identified
  if not(isfield(bf0,'sub'))
    error('Participant information is missing for this folder:\n %s\n',...
      fp_dicomfolder)
  end
  pID = sprintf('sub-%s',bf0.sub);
  
  
  %% gather some information
  tab_ov = dicomOverview(fp_dicomfolder);
  if numel(unique(tab_ov.name))>1
    error('essbids:dicom2bids:MultipleSubjects',...
      ['There are multiple names detected in this dataset.\n'...
      'Please provide only data from one subject per function call.\n'...
      '%s'],evalc('disp(tab_ov)'));
  end
    
  
  %% update participants.tsv
  fn_ptab = fullfile(fp_bids,'participants.tsv');
  if exist(fn_ptab,'file')==2
    ptab = essbids_readTsv(fn_ptab);
  else 
    ptab = table('Size',[0,3],...
      'VariableTypes',{'cell','categorical','double'},...
      'VariableNames',{'participant_id','sex','age'});
    jinfo_ptab.sex.Description = 'participant sex';
    jinfo_ptab.sex.Levels = struct('m','male','f','female','o','other');
    jinfo_ptab.age.Description = 'participant age at first MRI session';
    jinfo_ptab.age.Units = 'years';
    ptab=addprop(ptab,{'JsonSidecar'},{'table'});
    ptab.Properties.CustomProperties.JsonSidecar = jinfo_ptab;
    essbids_writeTsv(fn_ptab,ptab)
  end
  ind_pid = ismember(ptab.participant_id,pID);
  if sum(ind_pid)==0
    wsts=warning('query','MATLAB:table:RowsAddedExistingVars');
    warning('off','MATLAB:table:RowsAddedExistingVars');
    ptab.participant_id{end+1}= pID;
    warning(wsts.state,'MATLAB:table:RowsAddedExistingVars');
    ptab.sex(end) = tab_ov.sex(1);
    ptab.age(end) = tab_ov.age(1);
    essbids_writeTsv(fn_ptab,ptab)
  else
    
  end
  
  %% convert dicoms or use temp folder
  dir_temp = fullfile(fp_bids,pID,'temp');
%   % build dcm2niix command 
%   myCmd = sprintf(['%s -9 '...
%     '-a n '...
%     '-b y '...
%     '-d 5 ' ...
%     '-e n '...
%     '-f %s_%%d ' ...
%     '-g n '...
%     '-i n '...
%     '-l n '...
%     '-m 2 '...
%     '-o %s '...
%     '-p y '...
%     '-r n '...
%     '-s n '...
%     '-t n '...
%     '-v 0 '...
%     '-w 2 '...
%     '-x n '...
%     '-z y '...
%     '%s'],...
%     strrep(fn_dcm2niix,' ' ,'\ '),...
%     pID, strrep(dir_temp,' ' ,'\ '),strrep(fp_dicomfolder,' ','\ '));
%   
%   if exist(dir_temp,'dir')~=7
%     fprintf('creating temporary directory for nii conversion:\n %s\n',...
%       dir_temp);
%     mkdir(dir_temp)
%     system(myCmd);
%   else
%     dl =  struct2table(dir(dir_temp));
%     dl = dl(~ismember(dl.name,{'.','..'}),:);
%     if isempty(dl)
%      system(myCmd);
%     else
%      warning(['temp folder not empty, working with content instead '...
%        'of converting\n  %s'],dir_temp);
%     end
%   end
  
  %% process the files in the temporary folder
  if exist(fullfile(dir_temp, 'fl.tsv'), 'file')
    fl = readtable(fullfile(dir_temp, 'fl.tsv'), ...
        'Delimiter', '\t', 'FileType', 'text');
  else
%     fl = [struct2table(dir(fullfile(dir_temp,'*.nii*'))); ...
%           struct2table(dir(fullfile(dir_temp,'*.bval*'))); ...
%           struct2table(dir(fullfile(dir_temp,'*.bvec*')))];
    fl = struct2table(dir(fullfile(dir_temp,'*.nii*')));
  end
  if length(unique(fl.name))~=length(fl.name)
    error('renaming not done properly, duplicates probably in name column!')
  end
  
  fl_renamed = fl(cellfun(@(x,y) ~strcmp(x, y), fl.name, fl.name_old), :);
  if ~isempty(fl_renamed)
    fl_renamed.name_json = cellfun(@(x) strrep(x, '.nii.gz', '.json'), fl_renamed.name, 'uni', 0);
    fl_renamed.name_old_json = cellfun(@(x,y) strrep(x, '.nii.gz', '.json'), fl_renamed.name_old, 'uni', 0);
    cellfun(@(x, y, z) movefile(fullfile(x, y), fullfile(x, [z, '_TEMP'])), ...
        fl_renamed.folder, fl_renamed.name_old, fl_renamed.name, 'uni', 0)
    cellfun(@(x, y, z) movefile(fullfile(x, y), fullfile(x, [z, '_TEMP'])), ...
        fl_renamed.folder, fl_renamed.name_old_json, fl_renamed.name_json, 'uni', 0)
    % Use '_TEMP' suffix in order to make sure that you are not overwriting
    % any files in the renaming process
    cellfun(@(x, z) movefile(fullfile(x, [z, '_TEMP']), fullfile(x, z)), ...
        fl_renamed.folder, fl_renamed.name, 'uni', 0)
    cellfun(@(x, z) movefile(fullfile(x, [z, '_TEMP']), fullfile(x, z)), ...
        fl_renamed.folder, fl_renamed.name_json, 'uni', 0)
  end
  
  if isempty(fl)
    error('essbids:dicom2bids:TempFolder',[...
      'Temp folder is not empty, but contains no nii files.\n'...
      'Please check content and delete the folder before restarting\n'...
      'the program.\n  %s'],dir_temp)
  end
  fl = fl(~fl.isdir,:);
  if isempty(fl)
    fprintf('No nii files found in %s\n',dir_tmp)
    return
  end
  if iscell(fl.name)
    fl = cellfun(@fullfile,fl.folder,fl.name,'uni',0);
  else
    % char array in case of a single file
    fl = {fullfile(fl.folder,fl.name)};
  end
  cellfun(@delete,fl(contains(fl,files2ignore)));
  fl = fl(~contains(fl,[files2ignore,'niitable.csv']));
  % modify filennames according to replacement rules provided in input
  for i=1:size(filenameMods,1)
    for j=1:numel(fl)
      fn_old = fl{j};
      [fp,fn,fe] = essbids_fileparts(fn_old);
      fn_tmp=fn;
      while contains(fn_tmp,filenameMods{i,1})
        fn_tmp = regexprep(fn_tmp,filenameMods{i,:});
      end
      fn_new = fullfile(fp,[fn_tmp fe]);
      if not(isequal(fn_old,fn_new))
        movefile(fn_old,fn_new);
        fl{j} = fn_new;
        [~,fnn] = essbids_fileparts(fn_new);
        fl2 = dir(fullfile(fp,[ fn '.*']));
        for k=1:numel(fl2)
          fn_old2 = fullfile(fl2(k).folder,fl2(k).name);
          fn_new2 = fullfile(fl2(k).folder,strrep(fl2(k).name,fn,fnn));
          movefile(fn_old2,fn_new2)
        end
      end
    end
  end
  
  try
    fl_json = struct2table(dir(fullfile(dir_temp,'*')));
    fl_json = fl_json(~fl_json.isdir,:);
    fl_json = cellfun(@fullfile,fl_json.folder,fl_json.name,'uni',0);
    cellfun(@delete,fl_json(contains(fl_json,files2ignore)));
  catch
  end
  
  t1 = table(fl);
  t1.Properties.VariableNames{1} = 'fn_old';
%   t1.AcquisitionNumber = nan(numel(fl),1);
  t1.SeriesNumber = nan(numel(fl),1);
  t1.fn_new = repmat({''},numel(fl),1);
  t1.bidsflags = repmat({''},numel(fl),1);
  t1.jinfo     = repmat({''},numel(fl),1);
  
  %% first read through
  for i=1:numel(fl)
    fn_old = fl{i};
    %fn_old = strrep(fn_old,'_INV1','_inv-1');
    %fn_old = strrep(fn_old,'_INV2','_inv-2');
    %fn_old = strrep(fn_old,'_UNI','_UNIT1');
        
    bf1 = essbids_parseLabel(fn_old);
    if not(isfield(bf1,'type'))
      bf1.type = 'unknownType';
    end
    
    for j=1:numel(fn_add)
      bf1.(fn_add{j}) = bf_add.(fn_add{j});
    end
    if isfield(bf1,'suffix')
      suffix = bf1.suffix;
    else
      suffix = '';
    end
    switch suffix
      case 'INV1'
        bf1.inv = '1';
        bf1.suffix = 'MP2RAGE';
      case 'INV2'
        bf1.inv = '2';
        bf1.suffix = 'MP2RAGE';
      case 'UNI'
        bf1.suffix = 'UNIT1';
        other
      otherwise
        tmp = strsplit(bf1.fname,'_');
        tmp = tmp(ismember(lower(tmp),lower(types_suffix)));
        if not(isempty(tmp))
          bf1.suffix = types_suffix{ismember(lower(types_suffix),...
            lower(tmp(end)))};
        end
    end
   
    try
      if strcmpi(bf1.type,'dwi')&& contains(bf1.fname,'fmap')
        bf1.type='fmap';
        bf1.suffix = 'epi';
      end
    catch
    end
    
    t1.jinfo{i} = jsondecode(fileread(fullfile(bf1.fpath,...
      [bf1.fname '.json'])));
    if isfield(t1.jinfo{i},'ImageType')
      if ismember('phase',lower(t1.jinfo{i}.ImageType))
        %disp(fn_old)
        %disp(t1.jinfo{i}.ImageType)
        if strcmp(bf1.type,'fmap')
          bf1.suffix = 'phasediff';
           % will be corrected in next step in case of two phase images
        else
          bf1.part = 'phase';
        end   
      end
    else
      warning('No ImageType field for : %s',bf1.fname)
    end
    t1.fn_new{i} = essbids_buildFileName(bf1,'relativePath',true);
    t1.bidsflags{i}=bf1;
    if isfield(t1.jinfo{i},'SeriesNumber')
      t1.SeriesNumber(i) = t1.jinfo{i}.SeriesNumber;
    end
%     if isfield(t1.jinfo{i},'AcquisitionNumber')
%       t1.AcquisitionNumber(i) = t1.jinfo{i}.AcquisitionNumber;
%     end
  end

  
  %% remove unwanted series numbers
  ind_rmSer = ismember(t1.SeriesNumber,ignoreSerNum)
  if sum(ind_rmSer)>0
    fprintf('Removing %d series by manual demand: %s\n',...
      sum(ind_rmSer),strtrim(evalc('disp(t1.SeriesNumber(ind_rmSer)'')')))
    t1 = t1(not(ind_rmSer),:);
    fl2rm = fl(ind_rmSer);
    for i=1:numel(fl2rm)
      [fp,fn] = essbids_fileparts(fl2rm{i});
      delete(fullfile(fp,[fn '.*']));
    end
    fl = fl(not(ind_rmSer));
  end
  
      
  %% second read through combining different echos, phase and magnitude etc.
  for i=1:size(t1,1)
    if isfield(t1.bidsflags{i},'part')
      if not(strcmp(t1.bidsflags{i}.part,'phase'))
        continue
      end
      % check wether a corresponing mag exists
      for j=1:size(t1,1)
        if not(isequal(t1.jinfo{i}.ProtocolName,t1.jinfo{j}.ProtocolName))
          continue
        end
        if not(isequal(t1.jinfo{i}.ImageOrientationPatientDICOM,...
            t1.jinfo{j}.ImageOrientationPatientDICOM))
          continue
        end
        try
          ismag = ismember('m',lower(t1.jinfo{j}.ImageType));
        catch
          ismag = false;
        end
        if not(ismag)
          continue
        end
        timediff = seconds(...
          datetime(t1.jinfo{j}.AcquisitionTime,'InputFormat','HH:mm:ss.S') - ...
          datetime(t1.jinfo{i}.AcquisitionTime,'InputFormat','HH:mm:ss.S'));
        if timediff<5 
          t1.bidsflags{j}.part = 'mag';
        end
      end
    end
  end
  
    
  %% identify multi-echo sequences, fieldmaps, set proper echo number
  serNum = unique(t1.SeriesNumber);
  for i=1:numel(serNum)
    ind = find(ismember(t1.SeriesNumber,serNum(i)));
    if numel(ind)==1
      continue
    end
    if isfield(t1.jinfo{ind(1)},'EchoNumber')
      te = nan(size(ind));
      for j=1:numel(ind)
        te(j)=t1.jinfo{ind(j)}.EchoTime;
      end
      te_s = sort(unique(te),'ascend');
      for j=1:numel(te_s)
        ind1 = ind(ismember(te,te_s(j)));
        ind1 = ind1(1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(t1.bidsflags{ind1}.type,'fmap')
          if ismember('m',lower(t1.jinfo{ind1}.ImageType))
            t1.bidsflags{ind1}.suffix = sprintf('magnitude%1.0f',j);
            t1.jinfo{ind1}.EchoNumber = j;
            fprintf('%0.5f : magnitude%d\n',te_s(j),j)
          elseif ismember('phase',lower(t1.jinfo{ind1}.ImageType))
            t1.bidsflags{ind1}.suffix = sprintf('phase%1.0f',j);
            t1.jinfo{ind1}.EchoNumber = j;
            fprintf('%0.5f : phase%d\n',te_s(j),j)
          end
        else
          t1.bidsflags{ind1}.echo = sprintf('%.0f',j);
          t1.jinfo{ind1}.EchoNumber = j;
          if strcmpi(t1.bidsflags{ind1}.suffix,'gre')
            t1.bidsflags{ind1}.suffix = 'MEGRE';
          end
          fprintf('%0.5f : echo-%d\n',te_s(j),j)
        end        
      end
    end
  end
  
  % add required json field to phasediff in classic fieldmap
  for i=1:size(t1,1)
    if not(isfield(t1.bidsflags{i},'suffix'))
      continue
    elseif not(strcmp(t1.bidsflags{i}.suffix,'phasediff'))
      continue
    end
    % check wether a corresponing magnitude1 and magnitude2
    tmp=struct;
    for j=1:size(t1,1)
      if not(isfield(t1.bidsflags{j},'suffix'))
        continue
      elseif not(ismember(t1.bidsflags{j}.suffix,...
          {'magnitude1','magnitude2'}))
        continue
      end
      if not(isequal(t1.jinfo{i}.ProtocolName,t1.jinfo{j}.ProtocolName))
        continue
      end
      if not(isequal(t1.jinfo{i}.ImageOrientationPatientDICOM,...
          t1.jinfo{j}.ImageOrientationPatientDICOM))
        continue
      end
      try
        ismag = ismember('m',lower(t1.jinfo{j}.ImageType));
      catch
        ismag = false;
      end
      if not(ismag)
        continue
      end
      timediff = seconds(...
        datetime(t1.jinfo{j}.AcquisitionTime,'InputFormat','HH:mm:ss.S') - ...
        datetime(t1.jinfo{i}.AcquisitionTime,'InputFormat','HH:mm:ss.S'));
      if timediff<5 
        tmp.(t1.bidsflags{j}.suffix) = t1.jinfo{j}.EchoTime;
      end
    end
    if not(numel(fieldnames(tmp))==2)
      warning('essbids:dicom2bids:VolumeMismatch',...
        'Phasediff requires exactly two magnitude volumes, found: %d',...
        numel(fieldnames(tmp)));
      continue
    end
    t1.jinfo{i}.EchoTime1 = tmp.magnitude1;
    t1.jinfo{i}.EchoTime2 = tmp.magnitude2;
  end
  
%   %% update new names
%   for i=1:numel(fl)
%     t1.fn_new{i} = essbids_buildFileName(t1.bidsflags{i},...
%       'relativePath',true);
%   end
  
  %% update fieldmap json info with IntendedFor fields
  for i=1:size(t1,1)
    if not(strcmp(t1.bidsflags{i}.type,'fmap'))
      continue
    end
    bf = t1.bidsflags{i};
    if strcmp(bf.suffix,'epi')
      if isfield(bf,'task')
        types = {'func'};
      elseif contains(bf.fname,'dwi')
        types = {'dwi'};
      else
        types = {'func'};
      end
      if isfield(bf,'dir')
        bf.dir = bf.dir(end:-1:1);
      end
    else
      types = {'func'};
    end
    fn = fieldnames(bf);
    fn = fn(not(ismember(fn,...
      {'suffix','type','fname','fpath','extension'})));
    fl_IntendedFor = {};
    for j=1:size(t1,1)
      if not(ismember(t1.bidsflags{j}.type,types))
        continue
      end
      bf2 = t1.bidsflags{j};
      isMatch = true;
      for k=1:numel(fn)
        if isfield(bf2,fn{k})
          isMatch = isMatch && isequal(bf.(fn{k}),bf2.(fn{k}));
        else
          isMatch = false;
        end
      end
      if isMatch
        fl_IntendedFor = [fl_IntendedFor;t1.fn_new{j}]; %#ok<AGROW>
      end
    end
    for j=1:numel(fl_IntendedFor)
      % paths should be relative to subject folder
      tmp = strsplit(fl_IntendedFor{j},filesep);
      fl_IntendedFor{j} = fullfile(tmp{2:end}); %#ok<AGROW>
    end
    if isempty(fl_IntendedFor)
      warning('essbids:dicom2bids:fieldmap',...
        '%s : No files found to apply fieldmap to',t1.fn_old{i});
    else
      fprintf('Fieldmap %s intended for:\n%s',t1.fn_new{i},...
        evalc('disp(fl_IntendedFor)'))
    end
    t1.jinfo{i}.IntendedFor = sort(fl_IntendedFor);
    if isfield(t1.bidsflags{i},'task')
      if isfield(t1.bidsflags{i},'acq')
        t1.bidsflags{i}.acq = [t1.bidsflags{i}.acq t1.bidsflags{i}.task];
      else
        t1.bidsflags{i}.acq= t1.bidsflags{i}.task;
      end
      t1.bidsflags{i}=rmfield(t1.bidsflags{i},'task');
    end
  end
  
  %% update new names
  for i=1:numel(fl)
    t1.fn_new{i} = essbids_buildFileName(t1.bidsflags{i},...
      'relativePath',true);
  end
  
  %% pass along only uniquely labelled files of known type
  % leave other files in temp folder and save mat file with parsed metadata
  ind_unique = true(size(t1,1),1);
  for i=1:size(t1,1)
    ind = ismember(t1.fn_new,t1.fn_new(i));
    ind_unique(i) = sum(ind)==1;
  end
  t_excluded = t1(not(ind_unique),:);
  t1 = t1(ind_unique,:);
  if size(t_excluded,1)>0
    warning('No unique file names could be derived for:\n%s',...
      evalc('disp(t_excluded.fn_old)'));
  end
  ind_unknownType = cellfun(@(x) contains(x,...
    [ filesep 'unknownType' filesep]),t1.fn_new);
  if sum(ind_unknownType)>0
    t3 = t1(ind_unknownType,:);
    t1 = t1(not(ind_unknownType),:);
    warning('No known dataype could be identified for:\n%s',...
      evalc('disp(t3.fn_old)'));
    t_excluded = [t_excluded;t3];
  end
  if size(t_excluded,1)>0
    %writetable(t_excluded,fullfile(dir_temp,'niitable.csv'));
    save(fullfile(dir_temp,'niitable.mat'),'t_excluded');
  end
  
  
  
  
  
  %% moving the files
  fn_bi = fullfile(fp_bids,'.bidsignore');
  if exist(fn_bi,'file')==2
    txt_bi0 = fileread(fn_bi);
  else
    txt_bi0 = '';
  end
  txt_bi = txt_bi0;
  for i=1:size(t1,1)
    fn_old = t1.fn_old{i};
    [fp,fn] = essbids_fileparts(fn_old);
    fl2move = dir(fullfile(fp,[fn '.*']));
    
    [fp_new,fn_new] = essbids_fileparts(t1.fn_new{i});
    fp_new = fullfile(fp_bids,fp_new);
    if not(isfolder(fp_new))
      mkdir(fp_new)
    end
    for j=1:numel(fl2move)
      fn1 = fullfile(fl2move(j).folder,fl2move(j).name);
      [~,~,fe] = essbids_fileparts(fl2move(j).name);
      fn2 = fullfile(fp_new,[fn_new,fe]);
      movefile(fn1,fn2)
      if strcmpi(fe,'.json')
        jinfo_tmp = jsondecode(fileread(fn2));
        if not(isequal(jinfo_tmp,t1.jinfo{i}))
          fprintf('Updating json file : %s\n',fn2)
          essbids_writeJson(fn2,t1.jinfo{i})
        end
      elseif strcmp(t1.bidsflags{i}.type,'fmap') && ...
          (strcmpi(fe,'.bval') || strcmpi(fe,'.bvec'))
        % ignore bval and bvec files in fmap folders
        if not(contains(txt_bi,'*/fmap/*.bval'))
          txt_bi = sprintf('%s\n**/fmap/*.bval',txt_bi);
        end
        if not(contains(txt_bi,'**/fmap/*.bvec'))
          txt_bi = sprintf('%s\n**/fmap/*.bvec',txt_bi);
        end
      end
    end
  end
  if not(isequal(txt_bi,txt_bi0))
    fid=fopen(fn_bi,'w');
    fprintf(fid,'%s\n',strtrim(txt_bi));
    fclose(fid);
  end
  
  
  %% add task descriptors
  tasks = cell(size(t1,1),1);
  for i=1:size(t1,1)
    if isfield(t1.bidsflags{i},'task') && ...
        not(strcmpi(t1.bidsflags{i}.suffix,'sbref'))
      tasks{i} = [t1.bidsflags{i}.task '_' t1.bidsflags{i}.suffix];
    end
  end
  tasks = tasks(not(cellfun(@isempty,tasks)));
  tasks = unique(tasks);
  for i=1:numel(tasks)
    fn_task = fullfile(fp_bids,sprintf('task-%s.json',tasks{i}));
    if exist(fn_task,'file')==2
      continue
    end
    jinfo_task.TaskName = sprintf('long name for the task labeled %s',...
      tasks{i});
    jinfo_task.Instructions    = ...
      'Please fill the instructions given in here';
    jinfo_task.TaskDescription = ...
      'Please fill in a brief description of the task here';
    fprintf('New task descriptor file added: %s\n',fn_task)
    essbids_writeJson(fn_task,jinfo_task);
  end
  
  
  %% add dataset_description.json
  fn_descr = fullfile(fp_bids,'dataset_description.json');
  if not(exist(fn_descr,'file')==2)
    fprintf(['Generating new dataset_description.json file with dummy '...
      'content.\nPlease fill in needed information here.\n  %s\n'],...
      fn_descr);
    descr= struct('Name','Name of the experiment',...
      'BIDSVersion','1.6.0');
    descr.Authors = {'Author A. One','Author B. Two'};
    descr.Funding = {'Grant1 with ID number','Grant2 with ID number'};
    essbids_writeJson(fn_descr,descr);
  end
  
  
  %% write to readme
  fn_rm = fullfile(fp_bids,'README');
  fid = fopen(fn_rm,'a');
  fprintf(fid,'%s : data from dicoms added to participant %s\n',...
    datestr(now,'yyyy-mm-dd HH:MM:ss'),pID);
  fclose(fid);
  t_converted = t1;
  
  %% rename and move fl.tsv to record manual filename changes
  movefile(fullfile(dir_temp, 'fl.tsv'), ...
      fullfile(dir_temp, 'fl_used.tsv'))
  
  %% for each of the possible 4 QSM files, rename them and move them to anat
%   QSM_search_strings = {'Brainmask', 'Susceptibility', ...
%       'TissuePhase', 'UnwrappedPhase'};
%   QSM_replace_strings = {'acq-qsm_Brainmask', 'acq-qsm_Chimap', ...
%       'acq-qsm_TissuePhase', 'acq-qsm_UnwrappedPhase'};
%   for iQSM = 1:length(QSM_search_strings)
%       fp_QSM = dir(fullfile(dir_temp, ['*_', QSM_search_strings{iQSM}, '*nii*']));
%       if exist(fullfile(fp_QSM.folder, fp_QSM.name), 'file')
%           fp_QSM_new = strrep(fp_QSM.name, QSM_search_strings{iQSM}, QSM_replace_strings{iQSM});
%           fp_QSM_new = fullfile(fp_QSM.folder, '..', 'anat', fp_QSM_new);
%           fp_QSM_old = fullfile(fp_QSM.folder, fp_QSM.name);
%           movefile(fp_QSM_old, fp_QSM_new)
%           
%           % update t_excluded
%           t_exluded(ismember(t_excluded.fn_old, fp_QSM_old), :) = [];
%           
%           % same for json
%           fp_QSM_old = func_changeFileName_PrePostFileext(fp_QSM_old, '', '', 'json');
%           fp_QSM_new = func_changeFileName_PrePostFileext(fp_QSM_new, '', '', 'json');
%           movefile(fp_QSM_old, fp_QSM_new)
%       end
%   end
%   % Note: Chimap according to bids should not be in ppb but in ppm!

  %% cleanup
  if isempty(t_excluded)
      try
          delete(fullfile(dir_temp, 'niitable.mat'))
      catch
      end
  else
      save(fullfile(dir_temp, 'niitable.mat'), 't_excluded')
  end
  try
    rmdir(dir_temp);
  catch
  end
  
  
  
end


%% helper functions
function tab_ov = dicomOverview(fp_dicomfolder)
  
  dl = struct2table(dir(fullfile(fp_dicomfolder,'**')));
  dl = dl(dl.isdir==0,:);
  [~,~,dl.ext] = cellfun(@fileparts,dl.name,'uni',0);
  dl = dl(not(ismember(dl.name,'DICOMDIR')),:);
  dl = dl(ismember(lower(dl.ext),{'','.ima','.dcm'}),:);
  
  uni_dirs = unique(dl.folder);
  tic
  t0=table('Size',[numel(uni_dirs),4],...
    'VariableNames',{'name','sex','age','studydate'},...
    'VariableTypes',{'cell','categorical','cell','cell'});
  for i=1:numel(uni_dirs)
    ind = find(ismember(dl.folder,uni_dirs(i)),1,'first');
    fn_dcm = fullfile(dl.folder{ind},dl.name{ind});
    obj = images.internal.dicom.DICOMFile(fn_dcm);
    t0.name{i} = obj.getAttributeByName('PatientName');
    t0.sex(i)  = lower(obj.getAttributeByName('PatientSex'));
    t0.age{i}  = obj.getAttributeByName('PatientAge');
    t0.studydate{i} = [obj.getAttributeByName('StudyDate') 'T' ...
      obj.getAttributeByName('StudyTime')];
  end
  toc

  t0 = unique(t0);
  t0.age = cellfun(@age2double,t0.age);
  t0.studydate = cellfun(@datestr2datetime,t0.studydate);
  tab_ov = sortrows(t0,'studydate','ascend');
  
  
end

function ageDouble = age2double(ageString)
  % converts DICOM age default pattern ('025Y') to double, nan on error
  ageDouble = nan;
  try
    ageDouble = str2double(strrep(upper(ageString),'Y',''));
  catch
  end
end

function studyDateTime = datestr2datetime(studydateString)
  % converts studydate string to datetime
  studyDateTime = NaT;
  try
    studyDateTime = datetime(studydateString,'InputFormat',...
      'yyyyMMdd''T''HHmmss.S','Format','yyyy-MM-dd''T''HH:mm:SS');
  catch
  end
end