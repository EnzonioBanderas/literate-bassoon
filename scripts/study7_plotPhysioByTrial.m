function study7_plotPhysioByTrial()

  fn_physio = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata/sub-Z7T7179/ses-1/func/sub-Z7T7179_ses-1_task-fear_run-1_physio.tsv.gz';
  fn_events = strrep(fn_physio,'_physio.tsv.gz','_events.tsv');
  fn_eyetrack = '/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata2/sub-Z7T7179/ses-1/func/sub-Z7T7179_ses-1_task-fear_acq-stxtr1620_run-1_eyetrack.tsv.gz';
  
  fn_diamond = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/diamond.png';
  fn_square  = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/square.png';
  fn_fix     = '/media/diskEvaluation/Evaluation/sfb1280a05study7/misc/FixCross.png';
  img_d = imread(fn_diamond);
  img_s = imread(fn_square);
  img_m = (img_s+img_d)/2;
  %img_f = imread(fn_fix);

  ax_screen = subplot(5,2,1:4);
  imshow(img_m);
  hold on;
  %plot([0 1920],[0 1080])
  
  bf = essbids_parseLabel(fn_physio);
  myTitle = sprintf('sub-%s  ses-%s  sub-%s  trial-',bf.sub,bf.ses,bf.run);




  myXLims = [110,130];

  t_phys = essbids_readTsv(fn_physio);
  t_phys.time = t_phys.Properties.CustomProperties.Time;
  t_et = essbids_readTsv(fn_eyetrack);
  t_et.time = t_et.Properties.CustomProperties.Time;
  t_et.ps_A = t_et.EyeA_PupilHeight.*t_et.EyeA_PupilWidth;
  t_et.ps_B = t_et.EyeB_PupilHeight.*t_et.EyeB_PupilWidth;
  
  t_ev = essbids_readTsv(fn_events);
  uni_ev = cellstr( unique(t_ev.trial_type));
  t_ev2 = table('Size',[size(t_phys,1),numel(uni_ev)+1],...
    'VariableTypes',{'double','double','double'},...
    'VariableNames',['time',uni_ev']);
  t_ev2.time = t_phys.time;
  
  for i=1:numel(uni_ev)
    lo = 1.1 * (i-1);
    hi = lo + 1;
    t_ev2.(uni_ev{i})(:)=lo;
    ind = find(ismember(t_ev.trial_type,uni_ev(i)));
    for j=1:numel(ind)
      ind1 = t_ev2.time >= t_ev.onset(ind(j)) & ...
        t_ev2.time <= t_ev.onset(ind(j))+t_ev.duration(ind(j));
      t_ev2.(uni_ev{i})(ind1) = hi;
    end
  end


  

  
  ax_ps = subplot(5,2,5:6);
  plot(t_et.time,t_et.ps_A);
  hold on;
  plot(t_et.time,t_et.ps_B)
  ylabel('pupilsize [a.u.]');
  legend({'eyeA','eyeB'});
  

  ax_eda = subplot(5,2,7:8);
  plot(t_phys.time,t_phys.skinconductance);
  ylabel('EDA [µS]');
  
  ax_ev = subplot(5,2,9:10);
  for i=1:numel(uni_ev)
    plot(t_ev2.time,t_ev2.(uni_ev{i}));
    hold('on');
  end
  ylim([-0.2,numel(uni_ev)*1.1+0.1]);
  ylabel('events');
  legend(uni_ev);

 
  tIDs = unique(t_ev.trial_index);
  
  
  for i=1:numel(tIDs)
    ind0 = find(ismember(t_ev.trial_index,tIDs(i)),1,'first');
    myXLims = t_ev.onset(i) + [-5,15]; 
    axes(ax_eda);
    xlim(myXLims);
    axes(ax_ev);
    xlim(myXLims);
    axes(ax_ps);
    xlim(myXLims);
  
    axes(ax_screen);
    ind = t_et.time>= myXLims(1) & t_et.time<= myXLims(2);
    plot(t_et.EyeA_X_Gaze,t_et.EyeA_Y_Gaze)
    title(sprintf('%s%d',myTitle,tIDs(i)));
    pause
  end