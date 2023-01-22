function func_readBiopacMat(fp_s,fp_d)
% wrapper function to read biopac physio files to rawdata and create 
% first event list from biopac digital channels 
% 
% inputs : fp_s  -  path to sourcedata folder
%          fp_d  -  path to rawdata folder

  fprintf('\n%s\n%s - Collecting BIOPAC file data, save to "rawdata" \n',...
    repmat('=',72,1),datestr(now,'yyyy-mm-dd hh:MM:ss'))

  % settings for ter_BiopacPhysio2BIDSphysio
  addInfo_fear = 'task-fear_run-1';
  %addInfo_rest = 'task-rest_run-1';
  %minTrig = 40;     % minimum expected number of triggers within a session
  %TR_fear = 1.790;  % repetition time in fear conditioning 
  %TR_rest = 1.430;  % repetition time in resting state
  TR_fear = 0;
  minTrials = 5;
  minRunSpacing = 60; % seconds in between runs
  
  fl_mat = ter_listFiles(fp_s,'sub-*_ses-*_physio.mat');
  fl_mat = fl_mat(~contains(fl_mat,'logs'));
  fl_mat = fl_mat(~contains(fl_mat,'_tables'));
  fl_mat = fl_mat(contains(fl_mat,[filesep 'physio' filesep]));

  
  try
    ftab = struct2table(cellfun(@ter_parseFname,fl_mat));
  catch
    warning('Not all files follow the same pattern:\n');
    disp(fl_mat);
    ftab = table('Size',[0,2],'VariableType',{'cell','cell'},...
      'VariableName',{'pID','ses'});
    for j=1:numel(fl_mat)
      tmp = ter_parseFname(fl_mat{j});
      if isfield(tmp,'pID') && isfield(tmp,'ses')
        ftab(end+1,1) = {tmp.pID}; %#ok<AGROW>
        ftab(end,2)   = {tmp.ses};
      else
        warning('Could not properly identify ID and session for: \n%s\n',...
          fl_mat{j});
      end
    end
  end
  if not(isequal(size(unique(ftab)),size(ftab)))
    warning('Multiple identical file names for, skipping these')
    %2do implemenmt skipping, now just exit
    ftab_uni = unique(ftab);
    for j=1:size(ftab,1)
      ind_ft = ismember(ftab,ftab_uni(j,:));
      if sum(ind_ft)>1
        disp(ftab_uni(j,:))
      end
    end
    
    error('for now just error')
  end
  
  pIDs = ftab.pID;
  sIDs = ftab.ses;
  
  % skip sub-EN05NK27_day1.mat, in favor of sub-EN05NK27_day1_resampled.mat
  %fl_mat = fl_mat(~contains(fl_mat,'sub-EN05NK27_day1.mat'));
  %[~,fn_mat] = cellfun(@fileparts,fl_mat,'uni',0);
  %pIDs = cellfun(@(x) regexp(x,'sub-[0-9a-zA-Z]*','match'),fn_mat);
  %sIDs = cellfun(@(x) regexp(x,'ses-[0-9a-zA-Z]*','match'),fn_mat);
  %fl_out = cellfun(@(x,y) fullfile(fp_da,x,y,[x '_' y]),pIDs,sess,...
  fl_out = cellfun(@(x,y) fullfile(fp_d,x,y,'beh',[x '_' y '_task-fear']),...
    pIDs,sIDs,'uni',0);
  while max(contains(fl_out,'__'))==1
    fl_out = strrep(fl_out,'__','_');
  end

% for i=1:numel(fl_mat)
%   if exist([fl_out{i} '_physio.tsv.gz'],'file')~=2
%     fprintf('Converting mat to tsv: %s\n',fn_mat{i});
%     ter_BiopacPhysio2BIDSphysio(fl_mat{i},fp_d,TR,addInfo);
%   end
% end
  [~,fn_mat] = cellfun(@fileparts,fl_mat,'uni',0);
  %pIDs = cellfun(@(x) regexp(x,'sub-[0-9a-zA-Z]*','match'),fn_mat);
  %sIDs = cellfun(@(x) regexp(x,'ses-[0-9a-zA-Z]','match'),fn_mat);

  fp_out = cellfun(@(x,y) fullfile(fp_d,x,y,'beh'),pIDs,sIDs,...
    'uni',0);

  for i=1:numel(fl_mat)
    if ~isfolder(fp_out{i})
      mkdir(fp_out{i});
    end
    if numel(ter_listFiles(fp_out{i}, '*_physio.json')) == 0
      fprintf('Importing biopac mat to tsv: %s\n',fn_mat{i});
%    addInfo = 'task-fear_dir-AP';
%    if strcmpi(sIDs{i}, 'ses-2')
%      addInfo = 'task-fear_dir-AP_run-1';
%    end
      %[fn_tsv,physiotabs] = 
      ter_BiopacPhysio2BIDSphysio(fl_mat{i},fp_d,...
        TR_fear,addInfo_fear,minTrials,minRunSpacing);
%       nTrig = cellfun(@(x) sum(x.trigger),physiotabs);
%       fn_bold = strrep(fn_tsv,'_physio','_bold.nii.gz');
%       for j=1:numel(fn_tsv)
%         if exist(fn_bold{j},'file')~=2
%           warning('bold file not found: %s',fn_bold{j})
%         end
%         try
%           nVol = str2double(evalc(sprintf('!fslval %s dim4',...
%             fn_bold{j})));
%           if nVol ~= nTrig(j)
%             warning('mismatch %s: %d triggers, %d volumes',fn_tsv{j},...
%               nTrig(j),nVol);
%           end
%         catch
%         end
%       end
      
      %addInfo = 'task-rest_dir-AP';
%       [fn_tsv2,physiotabs2] = ter_BiopacPhysio2BIDSphysio(fl_mat{i},fp_d,...
%         TR_rest,addInfo_rest,minTrials);
%       
%       nTrig2 = cellfun(@(x) sum(x.trigger),physiotabs2);
%       fn_bold2 = strrep(fn_tsv2,'_physio','_bold.nii.gz');
%       for j=1:numel(fn_tsv2)
%         if exist(fn_bold2{j},'file')~=2
%           warning('bold file not found: %s',fn_bold2{j})
%         end
%         try
%           nVol = str2double(evalc(sprintf('!fslval %s dim4',...
%             fn_bold2{j})));
%           if nVol ~= nTrig2(j)
%             warning('mismatch %s: %d triggers, %d volumes',fn_tsv2{j},...
%               nTrig2(j),nVol);
%           end
%         catch
%         end
%       end
      
      %fn_tsv = [fn_tsv;fn_tsv2]; %#ok<AGROW>
      %save(strrep(fl_mat{i},'.mat','_tables.mat'),'fn_tsv','physiotabs');
    end
  end