
function y = winit(x,wintype)
% WINIT Applies a window to a signal.
%   y = winit(x,wintype)
%
%   Applies window to a signal. Includes extra windows that guarantee
%   perfect TDAC MDCT reconstruction.
%
%   x:       input signal
%   wintype: window name eg: 'hanning'(default), 'sinewin', 'kbdwin'
%            'lowin'
%   y:       windowed signal

% ------- winit.m ------------------------------------------
% Marios Athineos, marios@ee.columbia.edu
% http://www.ee.columbia.edu/~marios/
% Copyright (c) 2002 by Columbia University.
% All rights reserved.
% ----------------------------------------------------------

% Check arguments
if (nargin < 2)
    wintype = 'hann';
end

[flen,fnum] = size(x);
% Make column
if (flen==1)
    x = x(:);
end

switch wintype
case 'rectwintdac'
    % For TDAC in MDCT
    % Never used in practice only theoritical use
    w = rectwintdac(flen);
case 'trapezwin'
    % For TDAC in MDCT
    % Never used in practice only theoritical use
    w = trapezwin(flen);
case 'sinewin'
    % For TDAC in MDCT
    w = sinewin(flen);
case 'kbdwin'
    % For TDAC in MDCT
    w = kbdwin(flen,6.5);
case 'lowin'
    % For Low-Overlap TDAC in MDCT (transients)
    w = lowin(flen);
case 'hamming'
    w = hamming(flen);
    w_sqrt=w;
otherwise
    w = window(wintype,flen);
    w_sqrt = sqrt(w); %%% modified by Tiago ... Feb 13, 2007
end
% Tony's Trick produces a very large intermediate variable
% We use the sparse matrix trick instead (mathworks.com Tech-Note 1109)
% Multiplies each column of the signal x by w !!!
%y = diag(sparse(w)) * x; %%% commented out by Tiago ... Feb 13, 2007
y = diag(sparse(w_sqrt)) * x; %%% added by Tiago ... Feb 13, 2007
