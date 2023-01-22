function data_filtered = ter_bandpassAnalogBiopac(data,samplingfreq)

  Fs   = samplingfreq;    % Sampling Frequency in Hz
  N    =         2^16;    % Order
  Fc1  =         0.05;    % First Cutoff Frequency in Hz
  Fc2  =        10.00;    % Second Cutoff Frequency in Hz
  flag =      'scale';    % Sampling Flag
  
  % Create the window vector for the design algorithm.
  win   = blackman(N+1);
  b     = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
  Hd    = dfilt.dffir(b);
  delay = mean(grpdelay(Hd,N,Fs));
  
  % use appending the last value to compensate for group delay 
  % (like described with zero-padding in MATLAB bandpass function)
  %   data_filtered = filter(Hd,[data;zeros(delay,1)]);
  data_filtered = [ones(delay,1)*data(1);data;ones(delay,1)*data(end)];
  data_filtered = filter(Hd,data_filtered);
  data_filtered = data_filtered(2*delay+1:end);
  
end
  
