function Hd = FDAToolOutput_ACheby2LowpassFilter2
%FDATOOLOUTPUT_ACHEBY2LOWPASSFILTER2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.6 and the Signal Processing Toolbox 7.1.
% Generated on: 11-Sep-2018 17:11:07

% Chebyshev Type II Lowpass filter designed using FDESIGN.LOWPASS.

global sosM;
global sV;

% All frequency values are in Hz.
Fs = 8000;  % Sampling Frequency

Fpass = 90;          % Passband Frequency
Fstop = 150;         % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 80;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);

sosM = Hd.sosMatrix; sV = Hd.ScaleValues;
end