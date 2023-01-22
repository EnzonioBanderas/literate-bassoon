function func_defaceAnatomicals(fp_d,fp_dis)
  %% deface anatomicals using pydeface
  fprintf('\n%s\n%s : Defacing anatomical MRI volumes\n',...
    repmat('=',1,72),datestr(now,'yyyy-mm-dd hh:MM:ss'));
  fl = ter_listFiles(fp_d,'sub-*acq-mprage_*T1w.nii.gz');
  fp_disc = fullfile(fp_dis,'non-defaced_T1w');
  if ~isfolder(fp_disc)
    mkdir(fp_disc)
  end
  
  fprintf(['Processing %d files, this will take approximately '...
    '5-10 minutes per file'],numel(fl));
  
  for i=1:numel(fl)
    fp = fileparts(fl{i});
    fprintf('%s : Started %s\n',datestr(now,'hh:MM:ss'),...
      strrep(fl{i},fp,''));
    fl2 = ter_listFiles(fp,'*.nii.gz');
    fl2 = fl2(~ismember(fl2,fl(i)));

    myCmd = sprintf('pydeface %s --nocleanup',strrep(fl{i},' ','\ '));
    if numel(fl2)>0
      myCmd = [myCmd ' --applyto']; %#ok<AGROW>
      for j=1:numel(fl2)
        myCmd= sprintf('%s %s',myCmd,strrep(fl2{j},' ','\ '));
      end
    end
  
    disp(myCmd);
    system(myCmd);

    fl2 = [fl(i);fl2]; %#ok<AGROW>
    cellfun(@(x) movefile(x,fp_disc),fl2)
    cellfun(@(x) copyfile(strrep(x,'nii.gz','json'),fp_disc),fl2);
    cellfun(@(x) movefile(strrep(x,'nii.gz','json'),...
      strrep(x,'.nii.gz','_defaced.json')),fl2);
    fl2 = ter_listFiles(fp,'*_defaced*');
    for j=1:numel(fl2)
      [~,fn,fe] = fileparts(fl2{j});
      fn_new = strrep(fn,'_defaced','');
      binfo_acq = regexp(fn,'acq-[a-zA-Z0-9]*','match');
      if isempty(binfo_acq)
        tmp = strsplit(fn_new,'_');
        tmp = [tmp(1:end-1),'acq-defaced',tmp(end)];
        fn_new = tmp{1};
        for k=2:numel(tmp)
          fn_new = [fn_new '_' tmp{k}]; %#ok<AGROW>
        end
      else
        fn_new = strrep(fn_new,binfo_acq{end},[binfo_acq{end} 'defaced']);
      end
      movefile(fl2{j},fullfile(fp,[fn_new fe]));
    end
    fl3 = ter_listFiles(fp,'*_pydeface*');
    cellfun(@(x) movefile(x,fp_disc),fl3);
  end
  
end