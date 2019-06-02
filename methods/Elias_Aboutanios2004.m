function [ f0,f1 ] = Elias_Aboutanios2004( s,fs )
% Dichotomous search method for acquring the maximizer of the periodogram
% with modification that not requires the zero-padding 
%   s: the received signal; fs: the sampling frequency
%   f0: the carrier frequency; f1: the low frequency

N = length(s);

% user-defined parameters
frp = 1e-4;
m_rp = N*frp/fs;

fft_reso = fs/N; %the frequency resolution of FFT
s = BlackmanWindow(s);

global Timer;
Timer = 0;
tic;

S = fft(s); Y = conj(S).*S;

%% the carrier frequency estimation
% 1) coarse estimation
kl = round([1696,1996,2296,2596] / fft_reso);
kr = round([1704,2004,2304,2604] / fft_reso);
maximum = 0;
for i = 1:4
    [tmp,I] = max(Y(kl(i)+1:kr(i)+1));
    if (tmp(1)>maximum)
        maximum = tmp(1);
        m_hat = I(1)+kl(i)-1;
    end
end

f0 = binary_search(m_hat);

%% the low frequency estimation
% 1) coarse estimation
kl = round([f0-32,f0+8] / fft_reso);
kr = round([f0-8,f0+32] / fft_reso);
[~,I1] = max(Y(kl(1)+1:kr(1)+1));
[~,I2] = max(Y(kl(2)+1:kr(2)+1));
m_hat_l = I1(1)+kl(1)-1;
m_hat_r = I2(1)+kl(2)-1;
f1 = (binary_search(m_hat_r)-binary_search(m_hat_l))/2;

Timer = Timer + toc;
fprintf('time consuming EA %fs\n',Timer);

%% ========================* functions *==========================
    function f_hat = binary_search(m)        
        delta = 0.75;
        Y_m = Y(m);  Y_p = Y(m+2);
        if (Y_p>Y_m)
            Y_m = Y_lambda(m+1-2*delta);
            m = m + 1 - delta;
        else
            Y_p = Y_lambda(m-1+2*delta);
            m = m - 1 + delta;
        end
        
        while delta>m_rp
            Y_0 = Y_lambda(m); 
            delta=delta/2;
            if (Y_p>Y_m)
                Y_m = Y_0;
                m = m + delta;
            else
                Y_p = Y_0;
                m = m - delta;
            end
        end
        f_hat = m*fs/N;
    end
    function Y_L = Y_lambda(lambda)
        k = (0:N-1)';
        exp_vector = exp(-1j*2*pi*k*lambda/N);
        S_lambda = sum(s.*exp_vector);
        Y_L = conj(S_lambda)*S_lambda;
    end
end

