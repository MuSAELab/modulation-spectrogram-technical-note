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
% Reverberant speech
rebrv_file= './data/speech_reverberant.wav';
[revrb_x, revrb_fs] = audioread(rebrv_file);

%% Modulation spectrograms
clean_modspectr = strfft_modulation_spectrogram(clean_x(1: 17152), clean_fs, ...
                                                1024, 0.25*1024, 1, [], 16, []);

revrb_modspectr = strfft_modulation_spectrogram(revrb_x(1: 17152), revrb_fs, ...
                                                1024, 0.25*1024, 1, [], 16, []);
                                            
%% AMA
% Parameters

figure()
ax1 = subplot(2,2,1);
plot_modulation_spectrogram_data(clean_modspectr, [], [], [], [-105, -65], 'inferno')
title('Clean neutral speech')
ax2 = subplot(2,2,2);
plot_modulation_spectrogram_data(revrb_modspectr, [], [], [], [-105, -65], 'inferno')
title('RT-0.8sec')
xlabel([ax1, ax2], 'Modulation frequency (Hz)');
ylabel([ax1, ax2], 'Conventional frequency (Hz)');

%% Atlas
% [B,F,T] = MODSPECGRAM(    A,       Fs,  ASEWINDOW,BASEOVERLAP,BASENFFT,MODWINDOW,MODNFFT)
subplot(2,2,3)
[B,F,T]=modspecgram(clean_x,clean_fs, 1024,0.75*1024, 1024, 64,512);
imagesc(T,F,20*log10(abs(B)));
axis xy
caxis([-10 30])
colorbar
colormap('inferno')
xlabel('Modulation frequency (Hz)')
ylabel('Acoustic Frequeny (Hz)')
title('clean')

