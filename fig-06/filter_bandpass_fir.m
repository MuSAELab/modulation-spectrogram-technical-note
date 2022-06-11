function Hd = filter_bandpass_fir(fc1, fc2 ,fs)
%FILTER_BANDPASS_FIR Returns a discrete-time filter object.

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs     = fs;  % Sampling Frequency
Fstop1 = max([fc1 - 1.2, 0]);             % First Stopband Frequency
Fpass1 = max([fc1 - 0.5, 0.01]);               % First Passband Frequency
Fpass2 = min([fc2 + 0.5, fs/2]);              % Second Passband Frequency
Fstop2 = min([fc2 + 1.2, fs/2]);            % Second Stopband Frequency
Astop1 = 40;      % First Stopband Attenuation (dB)
Apass  = 1;       % Passband Ripple (dB)
Astop2 = 40;      % Second Stopband Attenuation (dB)
dens   = 20;              % Density Factor


% Compute deviations
Dpass = (10^(Apass/20)-1)/(10^(Apass/20)+1);
Dstop1 = 10^(-Astop1/20);
Dstop2 = 10^(-Astop2/20);

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]