function sc = func_rmEdaArtefacts(fn_physio)

  %'/media/diskEvaluation/Evaluation/sfb1280a05study7/dumpHereForSorting/rawdata/sub-Z7T7179/ses-1/func/sub-Z7T7179_ses-1_task-fear_run-2_physio.tsv.gz'
  
  %% read in data and important information
  fn_events = strrep(fn_physio,'_physio.tsv.gz','_events.tsv');
  t1= essbids_readTsv(fn_physio);
  et= essbids_readTsv(fn_events);
  sf = t1.Properties.CustomProperties.JsonSidecar.SamplingFrequency;
  sc = t1.skinconductance;
  sc0 = sc;
  t = t1.Properties.CustomProperties.Time;
  et = et(ismember(et.trial_type,'US'),:);
  bf = essbids_parseLabel(fn_physio);

  %plot first plot
  clf
  plot(t,sc0);
  title(bf.fname,'Interpreter','none');
  hold on;
  
  %% remove US artefacts
  fprintf('removing %d US artefacts\n',size(et,1));
  for i=1:size(et,1)
    ind = t>=et.onset(i) & t <= et.onset(i) + et.duration(i);
    plot(t(ind),sc(ind),'r')
    mean_prev = mean(sc(t>=et.onset(i)-0.050 & t<=et.onset(i)));
    mean_post = mean(sc(t>=et.onset(i)+0.250 & t<=et.onset(i)+.300)) ;
    ind2rep = t>=et.onset(i) & t<=et.onset(i)+.250;
    sc(ind2rep) = linspace(mean_prev,mean_post,sum(ind2rep));
  end
  %plot(t,sc1,'g')

  %% remove TR artefact
  TR = 1.620;
  t_smooth_pre  = 0.100;
  t_smooth_post = 0.100;
  t_pre  = -0.020;
  t_post =  0.150;

  for i=1:floor(max(t)/TR)
    t_low = (i-1)*TR +t_pre;
    mean_low = mean(sc( t>=t_low-t_smooth_pre & t<=t_low));
    t_high = (i-1)*TR +t_post;
    mean_high = mean(sc( t>=t_high & t<=t_high+t_smooth_post));
    ind = t<=t_high & t>=t_low;
    sc(ind) = linspace(mean_low,mean_high,sum(ind));
    plot(t(ind),sc(ind),'r')
  end

  %% finaly apply FIT lowpass filter
  N    =           62;    % Order
  Fc   =        10.00;    % Cutoff Frequency in Hz
  flag =      'scale';    % Sampling Flag
  
  % Create the window vector for the design algorithm.
  win   = blackman(N+1);
  b     = fir1(N, Fc/(sf/2), 'low', win, flag);
  Hd    = dfilt.dffir(b);
  delay = mean(grpdelay(Hd,N,sf));
  
  % use appending the last value to compensate for group delay 
  % (like described with zero-padding in MATLAB bandpass function)
  %   data_filtered = filter(Hd,[data;zeros(delay,1)]);
  sc_filtered = [ones(delay,1)*sc(1);sc;ones(delay,1)*sc(end)];
  sc_filtered = filter(Hd,sc_filtered);
  sc_filtered = sc_filtered(2*delay+1:end);
  sc = sc_filtered;
  

  plot(t,sc,'g');

  
%   tol = 0.04;
%   [~,loc]=findpeaks(s1,'minPeakProminence',.03,'MinPeakDistance',300);
%   i = 2;
%   while i< numel(loc)
%     if not(t(loc(i))-t(loc(i-1))>=TR-tol || t(loc(i+1))-t(loc(i))<=TR+tol)
%       loc(i) = [];
%     else
%       i=i+1;
%     end  
%   end
%   
%   i=1;
%   
%   while i<numel(loc)
%     td =  t(loc(i+1)) - t(loc(i));
%     if t(loc(i+1)) - t(loc(i)) > TR+tol
%       tmp = round(linspace(loc(i),loc(i+1),round(td/TR)+1));
%       tmp = tmp(2:end-1)';
%       loc = sort([loc;tmp]);
%       fprintf('%d TR interpolated\n',numel(tmp));
%     end
%     i=i+1;
%     
%   end
%   fprintf('removing %d TR artefacts\n',numel(loc))
%   
%   for i=1:numel(loc)
%     mean_pre  = mean(s1(loc(i)+(-180:-150)));
%     mean_post = mean(s1(loc(i)+(150:180)));
%     s1(loc(i)+(-150:150)) = linspace(mean_pre,mean_post,301);
%   end
%   plot(t,s1,'y')
% 
% 
%   sc = s1;
end




