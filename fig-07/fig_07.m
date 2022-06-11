% Script to generate plots in Figure 4 from:
%
% INSERT ARTICLE BIBLIOGRAPHY WHEN PUBLISHED
%
% Caption: Time domain representation (top plots) of a clean (left) and noisy 
% (right) y-axis accelerometer signal segment (sampled at 30 Hz) and their 
% corresponding spectra (middle plots) and modulation spectrograms (bottom plots).

%% Load data 
clean_data = load('./data/clean_y_acc.mat');
clean_x = clean_data.AccelY2;
noisy_data = load('./data/noisy_y_acc.mat');
noisy_x = noisy_data.AccelY2;
fs=30;
time_vector = [1:length(clean_x)]' ./ fs; 

%% Compute PSDs
clean_psd = rfft_psd(clean_x, fs);
noisy_psd = rfft_psd(noisy_x, fs);

%% Compute modulation spectrograms
clean_modspectr = strfft_modulation_spectrogram(clean_x(1: 136), fs, 8, 0.25*8, 16, [], 2, []);
noisy_modspectr = strfft_modulation_spectrogram(noisy_x(1: 136), fs, 8, 0.25*8, 16, [], 2, []);                               

%% Plot figure
figure('units','normalized','outerposition',[0 0 1 1])

% === Time series
ax1 = subplot(3,2,1);
plot(time_vector, clean_x, 'k');
title('Clean y-axis accelerometer signal')
ax2 = subplot(3,2,2);
plot(time_vector, noisy_x, 'k');
title('Noisy y-axis accelerometer signal')
linkaxes([ax1, ax2])
xlabel([ax1, ax2], 'Time (s)');

% === Power spectrum
ax1 = subplot(3,2,3);
plot(clean_psd.freq_axis, clean_psd.PSD, 'k');
ax2 = subplot(3,2,4);
plot(noisy_psd.freq_axis, noisy_psd.PSD, 'k');
linkaxes([ax1, ax2])
xlabel([ax1, ax2], 'Frequency (s)');
ylabel([ax1, ax2], 'Power');

% === Modulation spectrograms
color_lims = [-60, -30];
ax1 = subplot(3,2,5);
plot_modulation_spectrogram_data(clean_modspectr, [], [], [], color_lims, 'inferno');
ax2 = subplot(3,2,6);
plot_modulation_spectrogram_data(noisy_modspectr, [], [], [], color_lims, 'inferno');
xlabel([ax1, ax2], 'Modulation frequency (Hz)');
ylabel([ax1, ax2], 'Conventional frequency (Hz)');

