function func_correctDataset(fp_s,fp_dis,fp_d)
  
  % 2021-08-09 Escape signes caught in fprintf (strrep(...,'\','\\'))

  if nargin < 1
    fp_s = evalin('base','fp_s');
  end
  if nargin < 2
    fp_dis = evalin('base','fp_dis');
  end
  if nargin < 3 
    fp_d = evalin('base','fp_d');
  end
    
  
  fp_s   = strip(fp_s,  'right',filesep);
  fp_dis = strip(fp_dis,'right',filesep);
  fp_d   = strip(fp_d,  'right',filesep);


  subj2rename = {
    % old name         new name
    'sub-MM7EL15'      'sub-MM07EL15' % typo subject name
    'sub-EL06En19'     'sub-EL06EN19'
    'sub-NN06Nz06'     'sub-NN06NZ06'
    'sub-BA04DA18'     'sub-BI04DA18'
    'sub-CE07R018'     'sub-CE07RO18'
    'sub-CE06EL20'     'sub-CI06EL20'
    'sub-TT09CH06'     'sub-TT08CH06'
    'sub-yKE06EF23'    'sub-KE06EF23'
    };

  dicom2discard = {
    % pID            session    series numbers
    %'sub-HA06IM01'  'ses-1'     3:11   % started up with some resting state
    };
  logs2discard = {
    %'HA06IM01-Study3_phase-20201205-120457-VTQNA.log'
    'sub-RG11LD22_ses-1_run-2_datetime-20210527T170150_stimlog.txt' % very briefs
    'sub-MM07EL15_ses-3_run-1_datetime-20210513T174059_stimlog.txt'
    'sub-KO04PE18_ses-3_run-1_datetime-20210611T102225_stimlog.txt' % wrong session selected 
    'sub-KO04PE18_ses-2_run-1_datetime-20210611T102119_questlog.txt'
    'sub-KO04PE18_ses-2_run-1_datetime-20210611T102119_questanswers.txt'
    'sub-KO04PE18_ses-2_run-1_datetime-20210611T101858_stimlog.txt'
    'sub-AS06LD14_ses-2_run-1_datetime-20210615T120933_stimlog.txt'
    'sub-LI06LE12_ses-3_run-2_datetime-20210623T142123_stimlog.txt'
    'sub-LI06UZ04_ses-3_run-2_datetime-20210710T121825_stimlog.txt' % restarted before any event
    };
  
  %files to rename
  files2rename = {...
    %fullfile(fp_s,'sub-NA04RI07/ses-1/stimlogs/NA04RI07-Study3_phase-20210126-212638-ZPUWM.log') ...
    };
  
  %folders2rename
  folders2rename = {
    %fullfile(fp_s,'sub-TA08OS26','ses-1','DICOM','sub-TA08OS26_ser-012_func_bold_dir-AP_run-3_ses-1_task-fear') ...
    %fullfile(fp_s,'sub-TA08OS26','ses-1','DICOM','sub-TA08OS26_ser-012_func_bold_dir-AP_run-2_ses-1_task-fear') ;
    };
  
  %% subjects to remove
  subjects2remove = {
    'sub-blablub'   % dummy run
    'sub-test'   % dummy run
    'sub-Test1'   % dummy run
    'sub-timmann'   % dummy run
    'sub-ab05cd10'
    'sub-ab11cd12'
    'sub-ab12cd12'
    'sub-ab30cd01'
    'sub-ab31cd01'
    'sub-AB02CD01'
    'sub-AB02CD02'
    'sub-DE05HU70'
    'sub-EN07IT14'
    'sub-MM07EL08'
    %'sub-EN08NS16'  % ses-2 run-1 acquired in wrong direction 
    %'sub-ES06AR01'  % ses-1 tasks in wrong PE direction
    'sub-TE05TH12' % pilot / test dataset
    'sub-RG11LD22' % pilot / test dataset
    'sub-WA04RD30' % pilot / test dataset
    'sub-MM07EL15' % pilot / test dataset
    'sub-EN07IT24' % pilot / test dataset
    'sub-IN06AN02' % only session 1 acquired
    'sub-HA05ED14' % only session 1 acquired
    'sub-ER06OH20' % high DASS-21 score
    'sub-KE06EF23' % high DASS-21 score
    'sub-LI06LE12' % high DASS-21 score
    'sub-YA05US17' % high DASS-21 score
    'sub-NA05EY02' % no active substance in blood serum found
    'sub-IA05LL02' % no active substance in blood serum found
    'sub-NN06NZ06' % high dass-21 score
    };
  
  
  %% remove subjects
  for i=1:numel(subjects2remove)
    fp1 = fullfile(fp_s,subjects2remove{i});
    fp2 = fullfile(fp_dis,'sourcedata',subjects2remove{i});
    if exist(fp1,'dir')==7
      if exist(fp2,'dir')~=7
        if not(isfolder(fileparts(fp2)))
          mkdir(fileparts(fp2))
        end
        movefile(fp1,fp2)
      else
        fprintf(['Conflict removing sourcedata, please copy manually:\n'...
          '   %s\n'],fp1)
      end
    end
    fp1 = fullfile(fp_d,subjects2remove{i});
    fp2 = fullfile(fp_dis,'rawdata',subjects2remove{i});
    if exist(fp1,'dir')==7
      if exist(fp2,'dir')~=7
        if not(isfolder(fileparts(fp2)))
          mkdir(fileparts(fp2))
        end
        movefile(fp1,fp2)
      else
        fprintf(['Conflict removing rawdata, please copy manually:\n'...
          '   %s\ņ'],fp1)
      end
    end
  end
  fl_tsv = ter_listFiles(fp_d,'*.tsv',2);
  fl_tsv = fl_tsv(~contains(fl_tsv,'sub-'));
  for i=1:numel(fl_tsv)
    wtmp = warning('query','MATLAB:textscan:AllNatSuggestFormat');
    warning('off','MATLAB:textscan:AllNatSuggestFormat');
    try
      t1 = ter_readBidsTsv(fl_tsv{i});
      %readtable(fl_tsv{i},'filetype','text','delim','tab');
    catch
      fprintf('Cannot read file: %s\n',fl_tsv{i})
      continue;
    end
    warning(wtmp.state,'MATLAB:textscan:AllNatSuggestFormat');
    ind2rm = ismember(t1.participant_id,subjects2remove);
    if sum(ind2rm)>0
      t2 = t1(~ind2rm,:);
      t3 = t1(ind2rm,:);
      fn_rm = strrep(fl_tsv{i},fileparts(fp_d),fp_dis);
      if exist(fn_rm,'file')~=2
        if ~isfolder(fileparts(fn_rm))
          mkdir(fileparts(fn_rm));
        end
        writetable(t3,fn_rm,'filetype','text','delim','tab');
        writetable(t2,fl_tsv{i},'filetype','text','delim','tab');
      else
        t4 = ter_readBidsTsv(fn_rm);
        try
          try
            t3 = [t4;t3];
          catch
          end
          ter_writeBidsTsv(t3,fn_rm);
          ter_writeBidsTsv(t2,fl_tsv{i});
          %writetable(t3,fn_rm,'filetype','text','delim','tab');
          %writetable(t2,fl_tsv{i},'filetype','text','delim','tab');
        catch
          fprintf('trouble removing tsv info, please check manually: %s',...
            fl_tsv{i})
        end
      end
      
    end
  end
  
  
  
  
  
  %% correct ptab
  % multiple stimsettingratings 
  fn_ptab = fullfile(fp_d,'participants.tsv');
  ptab0 = ter_readBidsTsv(fn_ptab);
  ptab = ptab0;
  ind = ismember(ptab.participant_id,'sub-NN07RT25');
  ptab.stim_intensity(ind) = 3.24;
  ptab.stim_sensory_threshold(ind) =  0.7; 
  if not(isequaln(ptab0,ptab))
    ter_writeBidsTsv(ptab,fn_ptab)
  end
  
  
  %% 1st rename subjects
  % renames everything but dicoms
  for j=1:size(subj2rename,1)
    pID_old = subj2rename{j,1};
    pID_new = subj2rename{j,2};
    nold = strrep(pID_old,'sub-','');
    nnew = strrep(pID_new,'sub-','');
    isMod = false;
    try
      fp_discard2 = fullfile(fp_dis,'sourcedata_renamed');
      if ~isfolder(fp_discard2)
        mkdir(fp_discard2);
      end
      try
        copyfile(fullfile(fp_s,pID_old),fullfile(fp_discard2,pID_old));
      catch
      end
      if isfolder(fullfile(fp_s,pID_old))
        movefile(fullfile(fp_s,pID_old),fullfile(fp_s,pID_new));
      end
    catch
    end
    fl2rename = ter_listFiles(fp_s,['*' nold '*']);
    fl2rename = fl2rename(cellfun(@(x) exist(x,'file')==2,fl2rename));
    fl_rename2 = strrep(fl2rename,nold,nnew);
    for i=1:numel(fl2rename)
      if exist(fl_rename2{i},'file')==2
        delete(fl2rename{i})
      else
        if ~isfolder(fileparts(fl_rename2{i}))
          mkdir(fileparts(fl_rename2{i}));
        end
        movefile(fl2rename{i},fl_rename2{i});
        isMod=true;
      end
    end
    fl2mod = ter_listFiles(fp_s,{['*' nnew '*.txt'],['*' nnew '*.log']});
    for i=1:numel(fl2mod)
      cont = fileread(fl2mod{i});
      contmod = strrep(cont,nold,nnew);
      if ~isequal(cont,contmod)
        fprintf('modifing file :%s\n',fl2mod{i});
        fid = fopen(fl2mod{i},'w');
        fprintf(fid,strrep(contmod,'\','\\'));
        fclose(fid);
        isMod = true;
      end
    end
    if isMod
      func_sortFromUnsorted(fullfile(fp_s,pID_new),fp_s,fp_dis,...
        'skipdicoms','silent','move');
    end
  end
  
  %movefile('/media/diskEvaluation/Evaluation/sfb1280a05study3/sourcedata/sub-DI05CH11/ses-1/DICOM','/media/diskEvaluation/Evaluation/sfb1280a05study3/sourcedata/sub-DT05CH11/ses-1/DICOM');
  %movefile('/media/diskEvaluation/Evaluation/sfb1280a05study3/sourcedata/sub-DI05CH11/ses-2/DICOM','/media/diskEvaluation/Evaluation/sfb1280a05study3/sourcedata/sub-DT05CH11/ses-12/DICOM');
  %dl = ter_listFiles('/media/diskEvaluation/Evaluation/sfb1280a05study3/sourcedata/sub-DT05CH11','*sub-DI05CH11*');
  %for i=1:numel(dl)
  %  movefile(dl{i},strrep(dl{i},'sub-DI05CH11','sub-DT05CH11'));
  %end
  
  
  %% remove files 
  for i=1:size(dicom2discard,1)
    fp = fullfile(fp_s,dicom2discard{i,1},dicom2discard{i,2},'DICOM');
    for j=dicom2discard{i,3}
      fl = dir(fullfile(fp,sprintf('*ser-%03.0f*',j)));
      if ~isempty(fl)
        fp_old = fullfile(fl(1).folder,fl(1).name);
        fp_new = fullfile(fp_dis,'sourcedata',dicom2discard{i,1},...
          dicom2discard{i,2},'DICOM',fl(1).name);
        mkdir(fileparts(fp_new))
        movefile(fp_old,fp_new);
      end
    end
  end
  for i=1:numel(logs2discard)
    finfo = ter_parseFname(logs2discard{i});
    pID = finfo.pID;
    sID = finfo.ses;
    fp1 = fullfile(fp_s,pID,sID,'stimlogs');
    fp2 = fullfile(fp_dis,'sourcedata',pID,sID,'stimlogs');
    fn_old = fullfile(fp1,logs2discard{i});
    fn_new = fullfile(fp2,logs2discard{i});
    if exist(fn_old,'file')==2
      if ~isfolder(fileparts(fn_new))
        mkdir(fileparts(fn_new));
      end
      fprintf('moving to "discarded/sourcedata" folder:\n  %s\n',fn_old);
      movefile(fn_old,fn_new)
    end
  end
  
  
  %% remove superflous event in a event table
  fn = fullfile(fp_d,'sub-RG06AS24','ses-3','beh',...
    'sub-RG06AS24_ses-3_task-fear_run-2_events.tsv');
  if exist(fn,'file')==2
    t0 = readtable(fn,'filetype','text','delim','tab');
    ind = t0.onset==800.8245 & ismember(t0.trial_type,'ShockOut');
    if sum(ind)==1
      t1 = t0(ind~=1,:);
      writetable(t1,fn,'FileTYpe','text','delim','tab');
    end
  end
  % similar thing happened here in in between a CS- and CS+ presentation 
  % in late acquisition phase 
  fn = fullfile(fp_d,'sub-MP05KE06','ses-1','beh',...
    'sub-MP05KE06_ses-1_task-fear_run-2_events.tsv');
  if exist(fn,'file')==2
    t0 = readtable(fn,'filetype','text','delim','tab');
    ind = t0.onset>433.90 & t0.onset<434.00 & ...
      ismember(t0.trial_type,'ShockOut');
    if sum(ind)==1
      t1 = t0(ind~=1,:);
      writetable(t1,fn,'FileType','text','delim','tab');
    end
  end
  
    %     ind = ismember(t0.trial_type,'US_Digitimer_Out');
%     t1 = t0(~ind,:);
%     if ~isequal(t0,t1)
%       writetable(t1,fn,'filetype','text','delim','tab');
%     end
%   end
  
  
  %   fn = fullfile(fp_d,'sub-EN05DK02','ses-1','func',...
%     'sub-EN05DK02_ses-1_task-fear_dir-AP_run-3_events.tsv');
%   if exist(fn,'file')==2
%     t0 = readtable(fn,'filetype','text','delim','tab');
%     ind = ismember(t0.trial_type,'US_Digitimer_Out');
%     t1 = t0(~ind,:);
%     if ~isequal(t0,t1)
%       writetable(t1,fn,'filetype','text','delim','tab');
%     end
%   end
  

  %% rename files
  for i=1:size(files2rename,1)
    fn_old = files2rename{i,1};
    fn_new = files2rename{i,2};
    if exist(fn_old,'file')==2
      if exist(fn_new,'file')==2
        warning('file cannot be renamed, already exists : %s',fn_new);
      else
        movefile(fn_old,fn_new);
      end
    end
  end

  
  %% rename folders
  for i=1:size(folders2rename,1)
    fn_old = folders2rename{i,1};
    fn_new = folders2rename{i,2};
    if exist(fn_old,'dir')==7
      if exist(fn_new,'dir')==7
        warning('folder cannot be renamed, already exists : %s',fn_new);
      else
        movefile(fn_old,fn_new);
      end
    end
  end
  
  

  %% remove single US events before actual events
  fn_tsv = '/media/diskEvaluation/Evaluation/sfb1280a05study6/rawdata/sub-ES06RT24/ses-3/beh/sub-ES06RT24_ses-3_task-fear_run-2_events.tsv';
  t1 = ter_readBidsTsv(fn_tsv);
  if strcmp(t1.trial_type{1},'ShockOut')
    ter_writeBidsTsv(t1(2:end,:),fn_tsv);
  end
  
  
  
  
  %% correct EDA data scaling
  % acquisition software was configured to 10 microsiemens/V, 
  t_scale = readtable(fullfile(fp_s,'documents',...
    'scaleing_skinconductance_SiemensPerVolt.xlsx'),'sheet','scaling');
  fp_discard2 = fullfile(fp_dis,'before_skinconductance_rescaling');
  if not(isfolder(fp_discard2))
    mkdir(fp_discard2);
  end
  for i=1:size(t_scale,1)
    pID = t_scale.participant_id{i};
    if not(isfolder(fullfile(fp_d,pID)))
      continue
    end
    fl = ter_listFiles(fullfile(fp_d,pID),[pID '*_physio.tsv.gz']);
    for j=1:numel(fl)
      tmp = ter_parseFname(fl{j});
      if exist(fullfile(fp_discard2,[tmp.fname,tmp.extension]),'file')==2
        continue
      end
      vn = ['SiemensPerVolt_' tmp.ses([1:3,5]) '_' tmp.run([1:3,5])];
      sc = t_scale.(vn)(i)/10; % 
      if sc == 1
        continue
      end
      fprintf(['Subject %2.0f of %2.0f: %s Rescaling skinconductance '...
        'in %s by factor of %4.2f'],i,size(t_scale,1),pID,tmp.fname,sc);
      data = ter_readBidsTsv(fl{j});
      %plot(data.skinconductance);hold on;
      data.skinconductance = sc * data.skinconductance;
      %plot(data.skinconductance);hold on;
      movefile(fl{j},fp_discard2);
      ter_writeBidsTsv(data,fl{j});
      fprintf(' ... done.\n');
    end
  end
  % modifying one datset where skincondctance amplification was swichted 
  % within session
  fl = ter_listFiles(fullfile(fp_d,'sub-HT05UL11','ses-1','beh'),...
    'sub-HT05UL11*ses-1*run-2*_physio.tsv.gz');
  tmp = ter_parseFname(fl{1});
  if not(exist(fullfile(fp_discard2,[tmp.fname,tmp.extension]),'file')==2)
    fprintf('Subject %s : Correcting switch of amplification in %s',...
        pID,tmp.fname);
    data = ter_readBidsTsv(fl{1});
    ind1 = find(data.skinconductance==max(data.skinconductance),1,'first');
    ind2 = find(data.skinconductance==max(data.skinconductance),1,'last');
    ind2 = ind2+3;
    data.skinconductance(1:ind2) = 0.5 * data.skinconductance(1:ind2);
    data.skinconductance(ind1:ind2) = nan;
    movefile(fl{1},fp_discard2);
    ter_writeBidsTsv(data,fl{1});
    fprintf(' ... done.\n');
  end
  
  %% correct events tables
  % there were two typos in the original protocol files
  fl_e_hab = ter_listFiles(fp_d,'sub-*ses-1*run-1*_events.tsv');
  fp_discard2 = fullfile(fp_dis,'events_before_correction');
  if not(isfolder(fp_discard2))
    mkdir(fp_discard2);
  end
  for i=1:numel(fl_e_hab)
    et0 = ter_readBidsTsv(fl_e_hab{i});
    et=et0;
    et.trial_type(9)=strrep(et.trial_type(9),'trial-7','trial-8');
    if not(isequal(et,et0))
      fprintf('Correct events table: %s\n',fl_e_hab{i});
      movefile(fl_e_hab{i},fp_discard2);
      ter_writeBidsTsv(et,fl_e_hab{i});
    end
  end
  fl_e_acq = ter_listFiles(fp_d,'sub-*ses-1*run-2*_events.tsv');
  for i=1:numel(fl_e_hab)
    et0 = ter_readBidsTsv(fl_e_acq{i});
    et = et0;
    ind1 = find(contains(et0.trial_type,'trial-3_'));
    ind2 = find(contains(et0.trial_type,'trial-2_'));
    if ind1<ind2
      et.trial_type(ind1)=strrep(et.trial_type(ind1),'trial-3_','trial-2_');
      et.trial_type(ind2)=strrep(et.trial_type(ind2),'trial-2_','trial-3_');
    end
    if not(isequal(et,et0))
      fprintf('Correct events table: %s\n',fl_e_acq{i});
      movefile(fl_e_acq{i},fp_discard2);
      ter_writeBidsTsv(et,fl_e_acq{i});
    end
  end
  
  
end