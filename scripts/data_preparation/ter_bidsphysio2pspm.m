function [mlb,pspmout] = ter_bidsphysio2pspm(fn_physio,fp_pspm_eval,varargin)

%% preface 
mlb     = {};
pspmout = {};

if ~isfolder(fp_pspm_eval)
  error('please create the target folder')
end

% delimter-seperated text file import was introduced with v5.0.0
if isempty(which('pspm'))
  error('Please add PsPM 5.0.0 or above to MATLAB paths');
end
[~,v]=pspm_version;
if str2double(v(1:min(strfind(v,'.'))-1)) < 5
  error('Current PsPM version is %s, please update to 5.0.0 or above',v);
end

if iscell(fn_physio)
  cellfun(@(x) ter_bids2pspm(x,fp_pspm_eval),fn_physio)
  return
end

trim2options = {'trim2starttime','trim2fmritrig','notrim','trim2events'};
trim2selection = trim2options(ismember(trim2options,lower(varargin)));
if isempty(trim2selection)
  trim2selection = trim2options{1};
else
  trim2selection = trim2selection{1};
end


%% find and verify the three needed files physio/beh, json and events
if exist(fn_physio,'file')==2 && ~isfolder(fileparts(fn_physio))
  fn_physio = which(fn_physio); % ensures that filename contains the path
end
if isempty(fn_physio)
  error('physio file input invalid')
end

[~,fn,fe] = fileparts(fn_physio);
if isequal(lower(fe),'.gz')
  fn_json = strrep(fn_physio,'.tsv.gz','.json');
  fn_tsv = fullfile(fp_pspm_eval,fn);
elseif isequal(lower(fe),'.tsv')
  fn_json = strrep(fn_physio,'.tsv','.json');
  fn_tsv = fullfile(fp_pspm_eval,[fn fe]);
else
  error('physio filename should end with extension .tsv or .tsv.gz')
end
if exist(fn_json,'file')~=2
  error('json sidecar missing: %s',fn_json);
end

fn_ev = [fn_physio(1:max(strfind(fn_physio,'_'))) 'events.tsv'];
recID = regexp(fn_ev,'recording-[a-zA-Z0-9]*','match');
if ~isempty(recID)
  fn_ev = strrep(fn_ev,['_' recID{end}],'');
end
if exist(fn_ev,'file')~=2
  error('events table missing, cannot process pspm model without it')
end


%% prepare decompressed tsv file if needed
if exist(fn_tsv,'file')~=2
  if isequal(lower(fe),'.gz')
    gunzip(fn_physio,fp_pspm_eval);
  else
    copyfile(fn_physio,fn_tsv);
  end
end


%% read info from json sidecar and events table
jinfo = jsondecode(fileread(fn_json));
sf = jinfo.SamplingFrequency;
ch = jinfo.Columns;
t0 = jinfo.StartTime;

et = readtable(fn_ev,'filetype','text','delim','tab','treat','n/a');
t0_events = min(et.onset) -2 -t0;
tmax_events = max(et.onset + et.duration)+15 -t0;

switch trim2selection
  case 'trim2fmritrig'
    pt = readtable(fn_tsv,'filetype','text','delim','tab','treat','n/a');
    pt = table2array(pt(:,ismember(lower(ch),'trigger')));
    t0_fmri   = (find(pt,1,'first')-1)/sf;
    %tmax_fmri = (find(pt,1,'last')-1)/sf;
end


%% setup and run data import
dsv.datafile = cellstr(fn_tsv);
for i=1:numel(ch)
  channel_name = lower(ch{i});
  switch channel_name
    case 'trigger'
      dsv.importtype{i}.marker.chan_nr.chan_nr_spec = i;
      dsv.importtype{i}.marker.flank_option         = 'ascending';
      dsv.importtype{i}.marker.sample_rate          = sf;
    case 'cardiac'
      dsv.importtype{i}.ecg.chan_nr.chan_nr_spec    = i;
      dsv.importtype{i}.ecg.sample_rate             = sf;
    case 'skinconductance'
      dsv.importtype{i}.scr.chan_nr.chan_nr_spec    = i;
      dsv.importtype{i}.scr.sample_rate             = sf;
      dsv.importtype{i}.scr.scr_transfer.none       = true;
    case 'pulseoximeter'
      dsv.importtype{i}.ppu.chan_nr.chan_nr_spec    = i;
      dsv.importtype{i}.ppu.sample_rate             = sf;
    case 'respiratory'
      dsv.importtype{i}.resp.chan_nr.chan_nr_spec   = i;
      dsv.importtype{i}.resp.sample_rate            = sf;
    otherwise
      warning('channel mapping not yet defined for channel %s',...
        channel_name);
  end
end
dsv.delimiter = '\t';
dsv.header_lines = 0;
dsv.channel_names_line = 0;
dsv.exclude_columns = 0;
mlb{1}.pspm{1}.prep{1}.import.datatype.dsv = dsv;
mlb{1}.pspm{1}.prep{1}.import.overwrite = true;
pspmout(1) = pspm_jobman('run',mlb(1));


%% setup and run trimming if needed
switch trim2selection
  case 'trim2starttime'
    t_offset = 0;
  case 'notrim'
    t_offset = -t0;
  case 'trim2fmritrig'
    t_offset = t0_fmri+t0;
  case 'trim2events'
    t_offset = t0_events+t0;
end
switch trim2selection
  case 'trim2starttime'
    % trimming output file to timepoint 0
    mlb{2}.pspm{1}.prep{1}.trim.datafile          = pspmout{1};
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_file.from = -t0;
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_file.to   = Inf;
    mlb{2}.pspm{1}.prep{1}.trim.overwrite         = true;
  case 'trim2fmritrig'
    % alternatively trim to fmri trigger pulses, should be the same at
    % front, but also trim the end 
    mlb{2}.pspm{1}.prep{1}.trim.datafile          = pspmout{1};
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_mrk.from = 0;
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_mrk.to   = 2;
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_mrk.mrk_chan.chan_def = 0;
    mlb{2}.pspm{1}.prep{1}.trim.overwrite        = true;
  case 'trim2events'
    mlb{2}.pspm{1}.prep{1}.trim.datafile          = pspmout{1};
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_file.from = t0_events;
    mlb{2}.pspm{1}.prep{1}.trim.ref.ref_file.to   = tmax_events;
    mlb{2}.pspm{1}.prep{1}.trim.overwrite         = true;
end
if numel(mlb)==2
  warning('off','ID:marker_out_of_range');
  pspmout(2) = pspm_jobman('run',mlb(2));
  warning('on','ID:marker_out_of_range');
end


%% delete tsv file in pspm folder if not identical with source file
if ~isequal(fn_tsv,fn_physio)
  delete(fn_tsv);
end


% %% save onsets file for GLM analysis
% fn_data = pspmout{end}{1};
% fn_onsets = [fn_data(1:max(strfind(fn_data,'_'))) 'onsets.mat'];
% names  = unique(et.trial_type);
% onsets = cell(size(names));
% %durations = cell(size(names));
% for i=1:numel(names)
%   onsets{i} = et.onset(ismember(et.trial_type,names(i)))-t_offset;
%   %durations{i} = et.duration(ismember(et.trial_type,names(i)));
% end
% save(fn_onsets,'names','onsets')%,'durations'


%% save temporal trimming offset to text file
fn_data = pspmout{end}{1};
fn_t0   = [fn_data(1:max(strfind(fn_data,'_'))) 't0.txt'];
fid = fopen(fn_t0,'w');
fprintf(fid,'%d',t_offset);
fclose(fid);



