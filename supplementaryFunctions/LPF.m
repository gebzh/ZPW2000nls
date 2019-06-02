function s_out = LPF( s )
% using a Chebyshev type II low-pass filter that has been designed using
% MATLAB FDA Tool to filter signal s (see 
% FDAToolOutput_ACheby2LowpassFilter()). For the convenience of 
% transplanting to embedded device, the low-pass filtering process is
% opened up which includes several solving processes of second order 
% difference equations.

%   s: the received signal; 

%   s_out: the filtered signal

% note that,
% realization of the filter function. In contrast to filtfilt that
% conducts zero-phase filtering, the filtered signal has some phase lags
% with the filter function.

global sosM;
global sV;

N = length(s);
s1 = s;
s2 = zeros(size(s));

if (isempty(sosM) || isempty(sV))
    FDAToolOutput_ACheby2LowpassFilter2();
end
TMP = 0;
for i = 1:size(sosM,1)
    s2(1) = sosM(i,1)*s1(1);  
    s2(2) = -sosM(i,5)*s2(1)+sosM(i,1)*s1(2)+sosM(i,2)*s1(1); 
    for j = 3:N
        s2(j) = -sosM(i,5)*s2(j-1)-sosM(i,6)*s2(j-2)+sosM(i,1)*s1(j)+...
            sosM(i,2)*s1(j-1)+sosM(i,3)*s1(j-2); 
    end
    s1 = s2;
end

Ak = 1;
for i = 1:length(sV)
    Ak = Ak * sV(i);
end

s_out = s1 * Ak;

end

