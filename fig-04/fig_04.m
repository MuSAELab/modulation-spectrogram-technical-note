% Script to generate plots in Figure 4 from:
%
% INSERT ARTICLE BIBLIOGRAPHY WHEN PUBLISHED
%
% Caption: Modulation spectrogram for speech in clean (left) 
% and reverberant (right, RT60 = 0.8 s) environments.

%% Load data
% Clean speech
clean_file= './data/speech_clean.wav';
[clean_x, clean_fs] = audioread(clean_file);
n_samples = length(clean_x);
% Reverberant speech
rebrv_file= './data/speech_reverberant.wav';
[revrb_x, revrb_fs] = audioread(rebrv_file);
revrb_x = revrb_x(1:n_samples);

%% Modulation spectrograms
clean_modspectr = strfft_modulation_spectrogram(clean_x, clean_fs, ...
                                                1024, 0.25*1024, 1, [], 16, []);

revrb_modspectr = strfft_modulation_spectrogram(revrb_x, revrb_fs, ...
                                                1024, 0.25*1024, 1, [], 16, []);
                                            
%% Plot figure
% === Figure
figure('units','normalized','outerposition',[0 0 1 1])

% === Modulation spectrograms
ax1 = subplot(1,2,1);
plot_modulation_spectrogram_data(clean_modspectr, [], [], [], [-105, -65], 'inferno')
title('Clean neutral speech')
ax2 = subplot(1,2,2);
plot_modulation_spectrogram_data(revrb_modspectr, [], [], [], [-105, -65], 'inferno')
title('RT-0.8sec')
xlabel([ax1, ax2], 'Modulation frequency (Hz)');
ylabel([ax1, ax2], 'Conventional frequency (Hz)');

