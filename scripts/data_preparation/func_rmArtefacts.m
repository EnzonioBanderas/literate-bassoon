function sc = func_rmArtefacts(fn_physio)

%'/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata/sub-Z7T7179/ses-1/func/sub-Z7T7179_ses-1_task-fear_run-2_physio.tsv.gz'
  fn_events = strrep(fn_physio,'_physio.tsv.gz','_events.tsv')
  t1= essbids_readTsv(fn_physio);
  et= essbids_readTsv(fn_events);

  sf = t1.Properties.CustomProperties.JsonSidecar.SamplingFrequency;
  sc = t1.skinconductance;
  t = t1.Properties.CustomProperties.Time;
  et = et(ismember(et.trial_type,'US'),:);

  TR = 1.620;
  tol = 0.02;
  
  
  clf
  plot(t,sc)
  sc1 = sc;
  hold on
  fprintf('removing %d US artefacts\n',size(et,1));
  for i=1:size(et,1)
    ind = t>=et.onset(i) & t <= et.onset(i) + et.duration(i);
    plot(t(ind),sc(ind),'r')
    mean_prev = mean(sc(t>=et.onset(i)-0.050 & t<=et.onset(i)));
    mean_post = mean(sc(t>=et.onset(i)+0.250 & t<=et.onset(i)+.300)) ;
    ind2rep = t>=et.onset(i) & t<=et.onset(i)+.250;
    sc1(ind2rep) = linspace(mean_prev,mean_post,sum(ind2rep));
  end
  plot(t,sc1,'g')

  [~,loc]=findpeaks(sc1,'minPeakProminence',.05,'MinPeakDistance',300);
  i = 2;
  while i< numel(loc)
    if not(t(loc(i))-t(loc(i-1))>=TR-tol || t(loc(i+1))-t(loc(i))<=TR+tol)
      loc(i) = [];
    else
      i=i+1;
    end  
  end
  fprintf('removing %d TR artefacts\n',numel(loc))
  

  for i=1:numel(loc)
    mean_pre  = mean(sc1(loc(i)+(-180:-150)));
    mean_post = mean(sc1(loc(i)+(150:180)));
    sc1(loc(i)+(-150:150)) = linspace(mean_pre,mean_post,301);
  end
  plot(t,sc1,'y')

end




