function func_antsWrapper4study7(pID,fp_d)
  
  % work in progress as of 2022-07-27, please wait till i finish this enzo
  if nargin < 1
    pID = 'sub-Z7T7173';
  end
  if nargin < 2
    fp_d = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata';
  end
  fp0     = '/media/diskEvaluation/Evaluation/sfb1280a05study7/';
  fp_temp = '/media/dataVol0/Evaluation/temporary_files/sfb1280a05study7/';
  
  fp_home = pwd;
  fp_de = fullfile(fileparts(fp_d),'derivatives');
  fp_t1map = fullfile(fp_de,'MP2RAGE_utils',pID);
  fp_out = fullfile(fp_de,'ants',pID);

  % search terms
  strID_func  = [pID '*_task-fear*_bold.nii.gz']; % for func files
  strID_fmap  = [pID '*fear*_epi.nii.gz'];        % for fmap files
  strID_t1map = [pID '*calcT1mapT1_clean.nii.gz'];% for t1map files
  
  % directories for ants processing, header for ants scripts
  fp_ants = fullfile(fp0,'scripts','packages','ANTs','install','bin');
  fp_ascr = fullfile(fp0,'scripts','packages','sk_ants_scripts');
  antsScriptHdr = sprintf([
    'export ANTSPATH=%s/\n'...
    'export SKSCRIPTS=%s/\n'...
    'export PATH=${SKSCRIPTS}:${ANTSPATH}:$PATH\n\n'],fp_ants,fp_ascr);

  
  %% collect functional files
  fl_func = essbids_listFiles(fullfile(fp_d,pID,'**',strID_func));
  if numel(fl_func)==0
    error('no valid functionals found');
  end
  % check whether output already exists
  [~,fn,fe] = essbids_fileparts(fl_func);
  fl_final0 = cellfun(@(x,y) fullfile(fp_out,...
    [x '_bet_MoCorr_DistCorr' y]),fn,fe,'uni',0);
  if prod(cellfun(@(x) exist(x,'file')==2,fl_final0))==1
    fprintf('Subject already processed: %s\n',pID);
    return;
  end
  

 
  
  %% collect all the proper opposed phase images
  % create a cell with two columns, left bold images, right opposed phase 
  % fmaps. Keep this pattern up for the rest of this script
  fl_f1 = [fl_func,cell(size(fl_func))];
  fl_fmap = essbids_listFiles(fullfile(fp_d,pID,'**',strID_fmap));
  fl_func_relPath = cellfun(@(x) strsplit(x,[filesep pID filesep]),...
    fl_func,'uni',0);
  fl_func_relPath = cellfun(@(x) x{2},fl_func_relPath,'uni',0);
  for i=1:numel(fl_fmap)
    try
      jinfo = jsondecode(fileread(essbids_findJsonSidecar(fl_fmap{i})));
      fl_f1(ismember(fl_func_relPath,jinfo.IntendedFor),2)=fl_fmap(i);
    catch
      warning('trouble attributing fmap %s',fl_fmap{i});
    end
  end
  ind_empty = cellfun(@isempty,fl_f1(:,2));
  if prod(~ind_empty)==0
    fprintf('Missing fieldmaps for:\n');
    disp(fl_func(ind_empty))
    error('Incomplete set of opposed phases, doublecheck IntendedFor fields')
  end
 


  %% find t1 map
  fl_t1map = essbids_listFiles(fullfile(fp_t1map,'**',strID_t1map));
  if numel(fl_t1map)==0
    error('no t1 map found here: %s',fp_t1map);
  elseif numel(fl_t1map)>1
    disp(fl_t1map);
    warning('multiple t1 maps available, using first one');
  end
  fn_t1map0 = fl_t1map{1};



  %% build all the filnames needed while processing
  % following pattern describes above, left column func, right column fmap
  fp_t = fullfile(fp_temp,pID);
  [~,fn,fe] = essbids_fileparts(fl_f1);
  bf = essbids_parseLabel(fl_f1);
  fl_f2 = cellfun(@(x,y) fullfile(fp_t,[x y]),fn,fe,'uni',0);
  fl_av  = cellfun(@(x,y) fullfile(fp_t,[x '_avg' y]),fn,fe,'uni',0);
  fl_abr = cellfun(@(x,y) fullfile(fp_t,[x '_avgbrain' y]),fn,fe,'uni',0);
  fl_msk = cellfun(@(x,y) fullfile(fp_t,[x '_avgbrain_mask' y]),fn,fe,...
    'uni',0);
  fl_bet = cellfun(@(x,y) fullfile(fp_t,[x '_bet' y]),fn,fe,'uni',0);
  fl_dct = cellfun(@(x,y) fullfile(fp_t,[x '_bet_DistCorr_template0' ...
    y]),fn,fe,'uni',0);
  fl_final = cellfun(@(x,y) fullfile(fp_t,[x '_bet_MoCorr_DistCorr' y]),...
    fn(:,1),fe(:,1),'uni',0);
  [~,fn,fe] = essbids_fileparts(fn_t1map0);
  fn_t1map = fullfile(fp_t,[fn,fe]);


  %% create temporary directory and copy files to it
  if not(isfolder(fp_t))
    fprintf('Creating temporary folder : %s\n',fp_t);
    mkdir(fp_t);
  end
  cd(fp_t);
  % copy files into temporary folder
  for i=1:numel(fl_f1)
    if exist(fl_f2{i},'file')==2
      continue
    end
    copyfile(fl_f1{i},fp_t);
  end
  if exist(fn_t1map,'file')~=2
    copyfile(fn_t1map0,fp_t);
  end


  %% generate brain extracted functionals
  for i=1:numel(fl_f2)
    if exist(fl_bet{i},'file')==2
      continue
    end
    % average functional images (same for bold and fmaps)
    syscmd('fslmaths %s -Tmean %s',fl_f2{i},fl_av{i});
    % generate brain masks
    syscmd('bet %s %s -f 0.1 -m -o',fl_av{i},fl_abr{i});
    % dilate the brain masks
    syscmd('fslmaths %s -dilF %s',fl_msk{i},fl_msk{i});
    % generated brain extracted functionals
    syscmd('fslmaths %s -mul %s %s',fl_f2{i},fl_msk{i},fl_bet{i});
  end



  %% apply distortion correction
  for i=1:size(fl_f2,1)
    if exist(fl_dct{i},'file')==2
      continue
    end
    fn_scr = fullfile(fp_t,sprintf('script_distortioncorr_%d.sh',i));
    jinfo = jsondecode(fileread(essbids_findJsonSidecar(fl_f1{i})));
    TR = jinfo.RepetitionTime;
    % from script documentation: 
    % -f: Fixed Image (created if not specified. Specify as /path/to/data/fixed.nii.gz)
    % -x: Fixed Mask (created if not specified. Specify as /path/to/data/fixedMask.nii.gz)
    % -t: TR in seconds (e.g. 2.0)
    % -a: Path to 4D fMRI data (required. Specify as /path/to/data/fMRI.nii.gz)
    % -b: Path to opposite PE data (required. Specify as /path/to/data/oPE.nii.gz)
    % -n: Number of threads (default = '16'. Specify more threads if available.)
    % -s: staged processing
    fid = fopen(fn_scr,'w');
    fprintf(fid,['%s'...
      'sk_ants_Realign_Estimate.sh'...
      ' -n 20 -t %s -f %s -x %s -a %s -b %s'],...]
      antsScriptHdr,sprintf('%f',TR), ...
      fl_abr{i,1}, ...  % -f "fixed" image, i.e., average functional masked
      fl_msk{i,1}, ...  % -x brain mask
      fl_bet{i,1}, ...  % -a brain extracted functional (bold) volumes 
      fl_bet{i,2});     % -b brain extracted opposed phase (fmap) volumes)
    fclose(fid);
    syscmd('chmod +x %s',fn_scr);
    system(fn_scr);
  end

  
  %% calculate corregistration to anatomy 
  %
  list_prefixes = cellfun(@(x) fullfile(fp_t,[pID ...
    '_reg-t1map2func_ses-' x.ses '_run-' x.run '_']),bf(:,1),'uni',0);
  fl_affine = cellfun(@(x) [x '0GenericAffine.mat'],list_prefixes,'uni',0);
  for i=1:size(fl_dct,1)
    if exist(fl_affine{i},'file')==2
      continue
    end
    fn_scr = fullfile(fp_t,sprintf('script_calcCoreg_%d.sh',i));
    fid = fopen(fn_scr,'w');
    fprintf(fid,['%s'...
      'antsRegistrationSyN.sh \\\n'...
      '-e 42            `# Fix random seed to an int value` \\\n'...
      '-d 3             `# ImageDimension: 2 or 3, here 3 obviously` \\\n'...
      '-f %s `# Fixed image(s) or source image(s) or reference image(s)` \\\n'...
      '-m %s `# Moving image(s) or target image(s)`\\\n'...
      '-o %s `# OutputPrefix: A prefix that is prepended to all output files.` \\\n'...
      '-t r             `# transform type (default = ''s''), r=rigid(1 stage)`\\\n'...
      '-p f             `# precision type (default = ''d''), f: float, d: double` \\\n'...
      '-j 1             `# use histogram matching, 0=false or 1=true (default = 0)` \\\n'...
      '-n 12            `# Number of threads`\n\n'],...
      antsScriptHdr,fl_dct{i},fn_t1map,list_prefixes{i});
    fclose(fid);
    syscmd('chmod +x %s',fn_scr);
    system(fn_scr);
  end


  %% apply all deformations
  for i=1:size(fl_bet,1)
    if exist(fl_final{i},'file')==2
      continue
    end
    jinfo = jsondecode(fileread(essbids_findJsonSidecar(fl_f1{i})));
    TR = jinfo.RepetitionTime;
    fn_scr = fullfile(fp_t,sprintf('script_applyDeformations_%d.sh',i));
    fid = fopen(fn_scr,'w');
    % -u: Ref ITK Registration .mat/.txt file
    % -v: Ref ITK Registration .nii.gz warp file
    % -x: Final ITK Registration .mat/.txt file
    % -y: Final ITK Registration .nii.gz warp file -> used diffently in
    %       editTE version
    % -i: Use inverse of Final ITK
    % -r: Reference Image
    % -f: Anatomical Image
    % -t: TR in seconds (e.g. 2.0)
    % -a: Path to 4D fMRI data (e.g. /path/to/data/fMRI.nii.gz)
    % -n: Number of threads (default = '16')
    fprintf(fid,['%s'...
      'sk_ants_Realign_Reslice_editTE.sh \\\n'...
      '-t %s   `# TR in seconds` \\\n'...
      '-a %s   `# Path to 4D fMRI data` \\\n'...
      '-y %s   `# affine file` \\\n'...
      '-n 12   `# Number of threads` \n\n'],...
      antsScriptHdr,TR,fl_bet{i},fl_affine{i});
    fclose(fid);
    syscmd('chmod +x %s',fn_scr);
    system(fn_scr);
  end

  if not(isfolder(fp_out))
    mkdir(fp_out)
  end
  for i=1:numel(fl_final)
    if exist(fl_final0{i},'file')~=2
      fprintf('copy from temporary directory folder : %s\n',fl_final{i});
      copyfile(fl_final{i},fp_out);
    end
  end
  


  %% wrap up 
  cd(fp_home);
  %delete temporary older
  %rmdir('s',fp_t)

end



function syscmd(cmdString,varargin)
  par = varargin;
  for i=1:numel(par)
    par{i} = strrep(par{i},' ','\ ');
  end
  cmdstr = sprintf(cmdString,par{:});
  fprintf([cmdstr,'\n']);
  system(cmdstr);
end