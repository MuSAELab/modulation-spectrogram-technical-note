% Script to generate plots in Figure 2 from:
%
% INSERT ARTICLE BIBLIOGRAPHY WHEN PUBLISHED
%
% Caption: Plots of (top-left) three snippets of synthetic ECG corrupted 
% by noise at SNR levels of 30 dB, 5dB and -5dB; (bottom-left) their 
% respective time-frequency plots; and (right) the modulation spectrograms 
% for the three signals.

%% Load data
ecg_data = load('./data.mat');
ecg_x    = ecg_data.sig';
ecg_fs   = ecg_data.fs;
ecg_time = [1:length(ecg_x)]' ./ ecg_fs; 

% Data consist of three 5-seconds segments of synthetic ECG
% Each segment was contaminated with pink noise to have a SNR of 30, 5 and -5 dB
ecg_30 = ecg_x(( 0 * ecg_fs)+1 : ( 5 * ecg_fs));
ecg_05 = ecg_x(( 5 * ecg_fs)+1 : (10 * ecg_fs));
ecg_n5 = ecg_x((10 * ecg_fs)+1 : (15 * ecg_fs));

%% Empirical values for the STFFT transformation
win_size_sec  = 0.125;   % seconds
win_over_sec  = 0.09375; % seconds
win_size_smp  = round(win_size_sec * ecg_fs); % samples
win_over_smp  = round(win_over_sec * ecg_fs); % samples
win_shft_smp  = win_size_smp - win_over_smp;
nfft_factor_1 = 64;
nfft_factor_2 = 4;

%% Compute STFT spectrogram of entire signal
ecg_spectrogram = strfft_spectrogram(ecg_x, ecg_fs, win_size_smp, win_shft_smp, nfft_factor_1);

%% Modulation spectrograms
ecg_modspect_30 = strfft_modulation_spectrogram(ecg_30, ecg_fs, win_size_smp, win_shft_smp, ...
                                                nfft_factor_1, [], nfft_factor_2, []);
ecg_modspect_05 = strfft_modulation_spectrogram(ecg_05, ecg_fs, win_size_smp, win_shft_smp, ...
                                                nfft_factor_1, [], nfft_factor_2, []);
ecg_modspect_n5 = strfft_modulation_spectrogram(ecg_n5, ecg_fs, win_size_smp, win_shft_smp, ...
                                                nfft_factor_1, [], nfft_factor_2, []);

%% Plot figure
% Plot parameters for spectrogram and modulation spectrograms
color_map = 'inferno';
freq_lims = [0, 64];
% Plot parameters for modulation spectrograms
color_lims = [-100, -20];

% === Figure
figure('units','normalized','outerposition',[0 0 1 1])

% === Time series
subplot(5, 4, ([1:3,5:7]))
plot(ecg_time, ecg_x, 'Color', 'k');
hold on;
title('ECG signal with increasing pink noise');
xlabel('Time (s)');
ylabel('ECG amplitude');
% Identigy different SNR segments
text(2.5,  -1.8, 'SNR = 30 dB');
plot([5, 5], get(gca,'ylim'), 'r--');
text(7.5,  -1.8, 'SNR = 5 dB');
plot([10, 10], get(gca,'ylim'), 'r--');
text(12.5, -1.8, 'SNR = -5 dB');
% Add and hide colorbar, to align x-axis with plot below
set(colorbar,'visible','off'); 

% === Spectrogram
subplot(5, 4, ([13:15,17:19]))
plot_spectrogram_data(ecg_spectrogram, [], [], freq_lims, [], color_map);
title('ECG signal spectrogram')
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% === Modulation spectrograms
ax1 = subplot(3, 4, 4);
plot_modulation_spectrogram_data(ecg_modspect_30, [], freq_lims, [], color_lims, color_map);
title('ECG modulation spectrogram (SNR = 30 dB)')
ax2 = subplot(3, 4, 8);
plot_modulation_spectrogram_data(ecg_modspect_05, [], freq_lims, [], color_lims, color_map);
title('ECG modulation spectrogram (SNR = 5 dB)')
ax3 = subplot(3, 4, 12);
plot_modulation_spectrogram_data(ecg_modspect_n5, [], freq_lims, [], color_lims, color_map);
title('ECG modulation spectrogram (SNR = -5 dB)')
xlabel([ax1, ax2, ax3], 'Modulation frequency (Hz)');
ylabel([ax1, ax2, ax3], 'Conventional frequency (Hz)');

