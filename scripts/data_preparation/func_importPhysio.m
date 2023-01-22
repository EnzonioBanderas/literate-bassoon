function func_importPhysio(fp_d,fp_de)

fn_ptab = fullfile(fp_d,'participants.tsv');
ptab = readtable(fn_ptab,'filetype','text','delim','tab','treat','n/a');


suffixes_physio = {...
  '_physio.tsv',...
  '_physio.tsv.gz',...
  '_beh.tsv',...
  '_beh.tsv.gz'};

for i=1:size(ptab,1)
  pID = ptab.participant_id{i};
  dl = struct2table(dir(fullfile(fp_d,pID,'ses-*')));
  if isempty(dl)
    sess = {''};
  else
    sess = cellstr(dl.name);
  end
  for j=1:numel(sess)
    
    %% find event tables in func for beh subfolders
    sID = sess{j};
    fp2search = {
      fullfile(fp_d,pID,sID,'func');
      fullfile(fp_d,pID,sID,'beh')};
    fp2search = fp2search(cellfun(@isfolder,fp2search));
    if isempty(fp2search)
      continue
    end
    fn_ev = ter_listFiles(fp2search,'sub-*_events.tsv',0);
    if isempty(fn_ev)
      continue
    end
    fn_ph = cell(size(fn_ev));
    for k=1:numel(fn_ev)
      for l=1:numel(suffixes_physio)
        fn_test = strrep(fn_ev{k},'_events.tsv',suffixes_physio{l});
        if exist(fn_test,'file')==2
          fn_ph{k} = fn_test;
          continue;
        end
      end
    end
    ind_empty = cellfun(@isempty,fn_ph);
    if sum(ind_empty)>0
      warning('no physio files could be found for these event tables:')
      disp(fn_ev(ind_empty))
    end
    fn_ph = fn_ph(~ind_empty);
    %fn_ev = fn_ev(~ind_empty);
    if isempty(fn_ph)
      continue
    end
    
    fp_out = fullfile(fp_de,'pspm',pID,'data');
    if ~isfolder(fp_out)
      mkdir(fp_out)
    end
    for k=1:numel(fn_ph)
      [~,fn,fe] = fileparts(fn_ph{k});
      if strcmpi(fe,'.gz')
        [~,fn] = fileparts(fn);
      end
      fn_pspm = fullfile(fp_out, ['pspm_' fn '.mat']);
      if exist(fn_pspm,'file')~=2
        ter_bidsphysio2pspm(fn_ph{k},fp_out,'trim2events');
      end
    end
    
  end
  
end
  


