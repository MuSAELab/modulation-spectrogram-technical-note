% Script to generate plots in Figure 6 from:
%
% INSERT ARTICLE BIBLIOGRAPHY WHEN PUBLISHED
%
% Caption: Measuring breathing rate from the breathing-related 
% modulations in an ECG signal.

%% Load data
brData = load('./data/data.mat');

%% Process signals
% ECG and Breating data (x and r respectively)
x  = brData.ch0(1: 60000);
r  = brData.ch1(1: 60000);
fs = brData.fs;

% Frequency limits for QRS complex (Hz)
f_ini_qrs = 1;
f_fin_qrs = 40;

% FIR filters
filter_qrs = filter_bandpass_fir(f_ini_qrs, f_fin_qrs, fs);
filter_lp  = filter_lowpass_fir(2, fs);
xf = filtfilt(filter_qrs.Numerator, 1, x);
rf = filtfilt(filter_lp.Numerator, 1, r);

% Use only 60 seconds
xf = xf((fs*30)+1 : (fs*90));
rf = rf((fs*30)+1 : (fs*90));
xf = xf - mean(xf);
rf = rf - mean(rf);
rf = -rf; 

% Compute modulation spectrogram
x_ms = wavelet_modulation_spectrogram(xf, fs, 6, [0.5: 0.5: 40], [], [], {'ECG'});

% Compute average of power across QRS BW ()
[~, ix_f_ini_qrs] = min(abs(x_ms.freq_axis - f_ini_qrs)); 
[~, ix_f_fin_qrs] = min(abs(x_ms.freq_axis - f_fin_qrs)); 
avg_power = mean(x_ms.power_modulation_spectrogram(ix_f_ini_qrs:ix_f_fin_qrs, :));

%% Plot figure
% Plot parameters for spectrogram and modulation spectrograms
color_map = 'inferno';
mod_freq_lims = [0.1, 0.666];   % in Hz, and equal to 6 - 36 breaths per min


freq_lims = [0, 64];


% Plot parameters for modulation spectrograms
color_lims = [-100, -20];

% === Figure
figure('units','normalized','outerposition',[0 0 1 1])

% === Time series
ax = subplot(5,5,[1:5]);
% Plot only first 30 seconds
xfp = xf(1 : 30 * fs);
rfp = rf(1 : 30 * fs);
time_vector = (0:size(xfp,1)-1)./fs;
plot(time_vector, [xfp, (rfp*3)+1.5], 'k')
xlabel('Time (s)')
yticks([0 1.5])
yticklabels({'ECG','Breathing'})
ax.YGrid = 'on';

% === Modulation spectrogram
subplot(5,5,[7:9,12:14])
plot_modulation_spectrogram_data(x_ms, [], [0, 40], mod_freq_lims, [-60, -25], 'inferno')
xlabel('Modulation frequency (Hz)')
ylabel('Conventional frequency (Hz)')

% === Average of power across QRS BW ()
subplot(4,5,[17:19])
% Double bottom x axis, based on:
% https://www.mathworks.com/matlabcentral/answers/344747-plotting-two-x-axis-in-one-plot-but-both-at-the-bottom
% First axes
plot(x_ms.freq_mod_axis, avg_power, '.-k');
xlim(mod_freq_lims);
hAX=gca;                  % first axes, save handle
pos=get(hAX,'position');  % get the position vector
pos1=pos(2);              % save the original bottom position
pos(2)=pos(2)+pos1; 
pos(4)=pos(4)-0.01;       % raise bottom/reduce height->same overall upper position
set(hAX,'position',pos)   % and resize first axes
pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
axes(hAX(1))
xlabel('Modulation frequency (Hz)')
set(colorbar,'visible','off'); 
ylabel('Power')
% Second axes
hAX(2)=axes('position',pos,'color','none');  % create the second axes
xlim(hAX(2), 60*mod_freq_lims)
axes(hAX(2))
xlabel('breathing rate (breaths per minute)')
set(colorbar,'visible','off'); 

