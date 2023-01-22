%% preface
fp_s = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/sourcedata2';
ws1 = warning('query','essbids:InvalidFieldName');
ws2 = warning('query','essbids:MultiplesDataType');
warning('off','essbids:InvalidFieldName');
warning('off','essbids:MultiplesDataType');



%% sort presentation logs
dl = struct2table(dir('/media/ELHuserdata/TemporaryFileExchange/sfb1280a05study7/preslog/sub-*'));
dl = dl(not(dl.isdir),:);          % no directories
dl = dl(dl.datenum > now - 21,:);  % only recent date, i.e. past 3 weeks
dl = dl(not(contains(dl.name,{'test','dummy','sub-_'})),:);
fl = cellfun(@fullfile,dl.folder,dl.name,'uni',0);
for i=1:numel(fl)
  fn0 = fl{i};
  [~,fn,fe] = fileparts(fn0);
  
  bf = essbids_parseLabel(fn0);
  fn_new = fullfile(fp_s,['sub-' bf.sub],'preslog',[fn,fe]);
  if not(contains(bf.sub,'Z7T'))
    continue
  end
  if not(exist(fn_new,'file')==2)
    fprintf('copy %s \n to  %s\n',fn0,fn_new);
    if not(isfolder(fileparts(fn_new)))
      mkdir(fileparts(fn_new))
    end
    copyfile(fn0,fn_new)
  end
end



%% sort physio files
dl = struct2table(dir('/media/ELHuserdata/TemporaryFileExchange/sfb1280a05study7/physio_data/sub-*'));
dl = dl(not(dl.isdir),:);          % no directories
dl = dl(dl.datenum > now - 21,:);  % only recent date, i.e. past 3 weeks
dl = dl(not(contains(dl.name,{'test','dummy','sub-_','sub-Z7T7213_sub-Z7T7214'})),:);
fl = cellfun(@fullfile,dl.folder,dl.name,'uni',0);

for i=1:numel(fl)
  fn0 = fl{i};
  [~,fn,fe] = fileparts(fn0);
  bf = essbids_parseLabel(fn0);
  mySuffix = '';
  if not(isfield(bf,'suffix'))
    if strcmpi(bf.extension,'.acq')
      mySuffix = '_physio';
    else
      mySuffix = [ '_' bf.extension(2:end)];
    end
  end
  fn_new = fullfile(fp_s,['sub-' bf.sub],'physio_data',[fn,mySuffix,fe]);
  if not(contains(bf.sub,'Z7T'))
    continue
  end
  if not(exist(fn_new,'file')==2)
    fprintf('copy %s \n to  %s\n',fn0,fn_new);
    if not(isfolder(fileparts(fn_new)))
      mkdir(fileparts(fn_new))
    end
    copyfile(fn0,fn_new)
  end
end



%% sort eyetracker
%/media/ELHuserdata/TemporaryFileExchange/sfb1280a05study7/eyetracker_data/
dl = struct2table(dir('/media/ELHuserdata/TemporaryFileExchange/sfb1280a05study7/eyetracker_data/sub-*'));
dl = dl(not(dl.isdir),:);          % no directories
dl = dl(dl.datenum > now - 21,:);  % only recent date, i.e. past 3 weeks
dl = dl(not(contains(dl.name,{'test','dummy','sub-_'})),:);
fl = cellfun(@fullfile,dl.folder,dl.name,'uni',0);
for i=1:numel(fl)
  fn0 = fl{i};
  [~,fn,fe] = fileparts(fn0);
  
  bf = essbids_parseLabel(fn0);
  fn_new = fullfile(fp_s,['sub-' bf.sub],'eyetracker_data',[fn,fe]);
  if not(contains(bf.sub,'Z7T'))
    continue
  end
  if not(exist(fn_new,'file')==2)
    fprintf('copy %s \n to  %s\n',fn0,fn_new);
    if not(isfolder(fileparts(fn_new)))
      mkdir(fileparts(fn_new))
    end
    copyfile(fn0,fn_new)
  end
end



%% rawdata
dl = struct2table(dir('/media/ELHuserdata/TemporaryFileExchange/sfb1280a05study7/mri_rawdata/*.dat'));
dl = dl(not(dl.isdir),:);          % no directories
dl = dl(dl.datenum > now - 21,:);  % only recent date, i.e. past 3 weeks
dl = dl(not(contains(dl.name,{'test','dummy','sub-_'})),:);
fl = cellfun(@fullfile,dl.folder,dl.name,'uni',0);
% to do, previous header reader does currently not work properly
%rdinfo = ter_rawDataInfo(fl);

%% DICOM
% copy .zip
% extract
% preprocess func_preprocessDICOM.m

%% end
warning(ws1.state,ws1.identifier);
warning(ws2.state,ws2.identifier);
fprintf('reseting access rights\n');
evalc('!/media/diskEvaluation/Evaluation/sfb1280a05study7/script_setProjectRights.sh');


