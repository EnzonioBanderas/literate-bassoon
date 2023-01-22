function [pulse_signal_time, pulse_signal] = ...
    func_pulseTable2pulseSignal(pulse_table, pulse_length)
%FUNC_PULSETABLE2PULSESIGNAL Summary of this function goes here
%   Detailed explanation goes here
if ~exist('pulse_length', 'var')
    pulse_length = 30;
end
pulse_startTime = min(pulse_table.time);
pulse_ind = round(pulse_table.time-pulse_startTime)+1;

pulse_signal = zeros(1, max(pulse_ind) + 4 * pulse_length);
pulse_signal_time = (1:max(pulse_ind))+pulse_startTime;
pulse_signal_time = [(pulse_signal_time(1) - (pulse_length * 2 - 1) * mean(diff(pulse_signal_time))):mean(diff(pulse_signal_time)):pulse_signal_time(1), ...
    pulse_signal_time, ...
    pulse_signal_time(end):mean(diff(pulse_signal_time)):(pulse_signal_time(end) + (pulse_length * 2 -1) * mean(diff(pulse_signal_time)))];
pulse_ind = pulse_ind + (pulse_length * 2);

pulse_signal(pulse_ind') = 1;
pulse_kernel = [zeros(1, pulse_length), ones(1, pulse_length)];
pulse_signal = conv(pulse_signal, pulse_kernel, 'same'); % full makes more sense
end

