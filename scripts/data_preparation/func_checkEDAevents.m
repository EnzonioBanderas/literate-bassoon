function allOk = func_checkEDAevents(fp_de)

  fprintf('Checking EDA event tables in derivatives folder\n');
  fp_eda = fullfile(fp_de,'EDAevaluation');
  fl = ter_listFiles(fp_eda,'sub-*_eventTable.mat',3);
  t0 = table('size',[numel(fl),4],...
    'VariableType',{'cell','cell','cell','double'},...
    'VariableName',{'participant_id','ses','run','nEvents'});
  t0.nEvents = nan(size(t0.nEvents));

  tmp = cell(size(fl));
  for i=1:numel(fl)
    tmp{i} = load(fl{i});
    t0.nEvents(i) = size(tmp{i}.eventtable,1);
    tmp2 = ter_parseFname(fl{i});
    t0.participant_id{i} = tmp2.pID;
    t0.ses{i} = tmp2.ses;
    t0.run{i} = tmp2.run;
  end
  t0 = sortrows(t0,{'ses','run'});

  allOk = true;
  sID = unique(t0.ses);
  for i=1:numel(sID)
    ind_s = ismember(t0.ses,sID{i});
    rID = unique(t0.run(ind_s));
    for j=1:numel(rID)
      ind_r = find(ind_s & ismember(t0.run,rID{j}));
      ind_err = ind_r(not(t0.nEvents(ind_r)==median(t0.nEvents(ind_r))));
      if numel(ind_err)>0
        allOk = false;
        fprintf('Outlier(s) in %s %s (should be %d events)\n',...
          sID{i},rID{i},median(t0.nEvents(ind_r))),
        disp(t0(ind_err,:))
      end
    end
  end

  if allOk
    fprintf('   Ok.\n')
  end
  
end

