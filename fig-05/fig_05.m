% Script to generate plots in Figure 5 from:
%
% INSERT ARTICLE BIBLIOGRAPHY WHEN PUBLISHED
%
% Caption: Plots of time domain (top) noisy ECG signal (60 bpm) corrupted 
% with an SNR = −5dB and (bottom) its enhanced counterpart obtained after 
% modulation spectrum domain based filtering.

%% Auxiliary scripts
addpath('./support/');

%% Load data
data = load('./data/data.mat');
x = data.signal;
fs=256; % Hz
time_vector = [1:length(x)]' ./ fs; 

%% Enhance signal
x_enhanced = Enhancement_ECG(x,fs);

%% Plot figure
% === Figure
figure('units','normalized','outerposition',[0 0 1 1])

% === Time series
ax1 = subplot(2,1,1);
plot(time_vector, x)
title('Noisy ECG signal')
ax2 = subplot(2,1,2);
plot(time_vector, x_enhanced)
title('Enhanced ECG signal')
linkaxes([ax1, ax2], 'x');
xlabel([ax1, ax2], 'Time (s)');

