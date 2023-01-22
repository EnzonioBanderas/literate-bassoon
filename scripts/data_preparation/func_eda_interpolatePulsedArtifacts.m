function sc = func_eda_interpolatePulsedArtifacts(fn_physio)
%% FUNC_EDA_INTERPOLATEPULSEDARTIFACTS
% function is build to interpolate TR and US artefacts from EDA data
  

  %% preset 
  TR = 1.620;
  t_smooth  = 0.075;
  t_pre_TR  = 0.020;
  t_post_TR = 0.150;
  t_pre_US  = 0.020;
  t_post_US = 0.250;



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
  
  index_data2replace = false(size(sc));
  %% calculate areas to smooth out US artifacts
  fprintf('removing %d US artefacts\n',size(et,1));
  ind_US = false(size(sc));
  for i=1:size(et,1)
    ind = t>=et.onset(i) & t <= et.onset(i) + et.duration(i);
    ind_US(ind)=true;
    ind2rep = t>=et.onset(i) -t_pre_US & t<=et.onset(i)+t_post_US;
    index_data2replace(ind2rep)=true;
  end
  sc_tmp = sc;
  sc_tmp(~ind_US) = nan;
  plot(t,sc_tmp,'y')

  %% remove TR artefact

  for i=1:floor(max(t)/TR)
    t_low = (i-1)*TR -t_pre_TR;
    %mean_low = mean(sc( t>=t_low-t_smooth_pre & t<=t_low));
    t_high = (i-1)*TR +t_post_TR;
    %mean_high = mean(sc( t>=t_high & t<=t_high+t_smooth_post));
    ind = t<=t_high & t>=t_low;
    index_data2replace(ind)=true;
    %sc(ind) = linspace(mean_low,mean_high,sum(ind));
    %plot(t(ind),sc(ind),'r')
  end
  
  
  %% combine multiple smoothing segments when they are not far enoufg apart
  ind_ons = find([0;diff(index_data2replace)==1]);
  ind_end = find([0;diff(index_data2replace)==-1]);
  i = 1;
  while i<numel(ind_ons)
    td = t(ind_ons(i+1))-t(ind_end(i));
    if td <= t_smooth
      ind_ons(i+1) = [];
      ind_end(i)   = [];
    else 
      i=i+1;
    end
  end

  sc_rep = nan(size(sc));
  for i=1:numel(ind_ons)
    t1_low   = t(ind_ons(i))-t_smooth;
    t1_high  = t(ind_ons(i));
    mean_low = mean(sc(t>=t1_low & t<=t1_high));
    t2_low   = t(ind_end(i));
    t2_high  = t(ind_end(i))+t_smooth;
    mean_high = mean(sc(t>=t2_low & t<=t2_high));
    if isnan(mean_low)
      mean_low = mean_high;
    elseif isnan(mean_high)
      mean_high = mean_low;
    end
    t1 = median(t(t>=t1_low & t<=t1_high));
    t2 = median(t(t>=t2_low & t<=t2_high));
    if isempty(t1)
      t1=t(1);
    elseif isempty(t2)
      t2=t(end);
    end
    ind2rep = t>=t1 & t<=t2;
    sc_rep(ind2rep) = linspace(mean_low,mean_high,sum(ind2rep));
  end
  plot(t,sc_rep,'r')
   sc(~isnan(sc_rep)) = sc_rep(~isnan(sc_rep));
  


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
legend('raw','US','interpolation','final')
  
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




