function [fxdm,fxdp,fpad] = basetran_mod(x,fhop)
% BASETRAN Base transform of the modulation codec
%   [fxdm,fxdp] = basetran(x,fhop)
%
%   By returning a matrix version of the framed signal we enable
%
%   x:     input signal
%   fhop:  hop length is samples
%   fpad:  padding information needed for OLA
%   fxdm:  magnitude part
%   fxdp:  phase part

% ------- basetran.m ---------------------------------------
% Marios Athineos, marios@ee.columbia.edu
% http://www.ee.columbia.edu/~marios/
% Copyright (c) 2003 by Columbia University.
% All rights reserved.
% ----------------------------------------------------------

% Window and its length
% (Other TDAC windows are 'lowin', 'rectwintdac', 'trapezwin')
%win  = 'kbdwin';
flen = 4*fhop; % for 75% overlap

% Frame it with 50% overlap
%[fx,fpad] = linframe(x,fhop,flen,'sym');%% commented out Tiago, Feb13/07
[fx,fpad] = linframe(x,fhop,flen,'pad'); %% modified by Tiago, Feb13/07

% If odd # of frames, zero pad with one more frame and rewindow to make even
% if mod(size(fx,2),2);
%     % Pad ...
%     x(end+fhop) = 0;
%     % ... and reframe
%     %[fx,fpad] = linframe(x,fhop,flen,'sym');%% commented out Tiago,
%     %Feb13/07
%     [fx,fpad] = linframe(x,fhop,flen,'pad');%% modified by Tiago, Feb13/07
% end
clear x;

% Window it
fx = winit(fx);

fxd = fft(fx);

% Split in odd and even shifts
% (Careful with odd and even since we are 1-based) 
% fxo = fx(:,2:2:end);
% fxe = fx(:,1:2:end);
% clear fx;

% Make the Cosine and Sine transformation matrices
% Ts  = mdxtmtx(flen,'f',@sin,0);
% Tc  = mdxtmtx(flen,'f',@cos,0);

% Take the modified sine and cosine transforms
% fxs = Ts*fxo;
% fxc = Tc*fxe;
% clear fxo fxe Ts Tc;

% Combine them into a complex transform
% fxd = fxc + i*fxs;
% clear fxc fxs;

% Now play with mag and phase (Tiago)
fxdm = abs(fxd); 
%fxdm = abs(fxd)/length(fx); 
fxdp = angle(fxd);


