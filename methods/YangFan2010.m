function [ f0,f1 ] = YangFan2010( s,fs )
%Zoom-FFT is used to reduce the signal length. And the zero-padding after 
%the down-sampling processes of Zoom-FFT is the key process that EXACTLY 
%increases the frequency resolution.
%   s: the received signal; fs: the sampling frequency
%   f0: the carrier frequency; f1: the low frequency

global Timer;
Timer = 0;
tic;

f_shift = [1660 1960 2260 2560];

N = length(s);
reso_fft = fs/N;

% windowing
s = HannWindow(s);

S = fft(s); Y = conj(S).*S;

% coarse estimation
kl = round([1696,1996,2296,2596] / reso_fft)+1;
kr = round([1704,2004,2304,2604] / reso_fft)+1;
maximum = 0;
for i = 1:4
    [tmp] = max(Y(kl(i):kr(i)));
    if (tmp(1)>maximum)
        maximum = tmp(1);
        num_interest = i;
    end
end

% Zoom-FFT
n = (0:N-1)';  
i_n = downsample(LPF(s.*cos(2*pi*f_shift(num_interest)*n/fs)),26);
q_n = downsample(LPF(s.*sin(2*pi*f_shift(num_interest)*n/fs)),26);
s_zf = i_n - 1j*q_n;
 
% zero-padding
frp = 1e-4;
M = round(fs/(26*frp)); 
s_zf = [s_zf; zeros(M-length(s_zf),1)];
reso_after_zp = fs/26/M;

S = fft(s_zf); Y = conj(S).*S;

% centroid method applied to enhance the precision
kl = round(36/reso_after_zp)+1;
kr = round(44/reso_after_zp)+1;
[~,I] = max(Y(kl:kr));
m = kl + I(1) - 1;
mm = (m-6:m+4)';
f0 = reso_after_zp*sum(Y(m-5:m+5).*mm)/sum(Y(m-5:m+5));

kl = round([f0-32,f0+8] / reso_after_zp)+1;
kr = round([f0-8,f0+32] / reso_after_zp)+1;
[~,I1] = max(Y(kl(1):kr(1)));
[~,I2] = max(Y(kl(2):kr(2)));
m_l = I1(1)+kl(1)-1;
m_r = I2(1)+kl(2)-1;
mm = (m_l-6:m_l+4)';
fl = reso_after_zp*sum(Y(m_l-5:m_l+5).*mm)/sum(Y(m_l-5:m_l+5));
mm = (m_r-6:m_r+4)';
fr = reso_after_zp*sum(Y(m_r-5:m_r+5).*mm)/sum(Y(m_r-5:m_r+5));
f1 = (fr-fl)/2;
f0 = f0+f_shift(num_interest);

Timer = Timer + toc;
fprintf('time consuming YF %fs\n',Timer);

end

