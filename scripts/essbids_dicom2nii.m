function essbids_dicom2nii(fp_dicomfolder,fp_bids,additionalInfo,varargin)
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
  % build dcm2niix command 
  myCmd = sprintf(['%s -9 '... % %s and shared
    '-f %s_%%d ' ... % %s
    '-o %s '... % %s
    '-z y '... % default=n (n=no), shared
    '-p y '... % default
    '-r n '... % default
    '-s n '... % default
    '-v 0 '... % default
    '-w 2 '... % default
    '-x n '...% default
    '-g n '... % default
    '-i n '... % default
    '-a n '... % default
    '-b y '... % default
    '-d 5 ' ... % default
    '-e n '... % default
    '-t n '... % ?
    '-l o '... % default=o (o=original), actually using n
    '%s'],...
    strrep(fn_dcm2niix,' ' ,'\ '),...
    pID, strrep(dir_temp,' ' ,'\ '),strrep(fp_dicomfolder,' ','\ '));
%     '-m 2 '... % default
  
  if exist(dir_temp,'dir')~=7
    fprintf('creating temporary directory for nii conversion:\n %s\n',...
      dir_temp);
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
  
  
  %% save list of nii files as table with possible rename column
  fl = struct2table(dir(fullfile(dir_temp,'*.nii*')));
  fl.name_old = fl.name;
  fl = [fl(:, 1), fl(:, end), fl(:, 2:end-1)];
  writetable(fl, fullfile(dir_temp, 'fl.tsv'), ...
      'Delimiter', '\t', 'FileType', 'text')

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