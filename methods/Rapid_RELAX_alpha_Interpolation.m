function [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_alpha_Interpolation( x_bold,fs )
% demodulating the ZPW-2000 signal using the RELAX algorithm (J Li 1996)
% /*rapid version 2*/ using the interpolation method to reduce the 
% computational burden.
%   inputs:
%   x_bold: the received signal;  fs: the sampling frequency
%   outputs:
%   fs_hat: the estimated frequencies of components
%   Amplis_hat: the estimated amplitudes of components
%   f0_hat: the carrier frequency est; f1_hat: the low frequency est
%   k_Emax: the numbers of the components belonging to the ZPW-2000 signal
% of the largest energy.

% =======================* initiation *=======================
% definition of the standard carrier freqs, the standard low freqs 
global f0_std;
global f1_std;
global f_std;

if isempty(f0_std) || isempty(f1_std) || isempty(f_std)
    f0_std = [1701.4;1698.7;2001.4;1998.7;2301.4;2298.7;2601.4;2598.7];
    f1_std = (10.3:1.1:29)';
    f_std = [-29:1.1:-10.3, 0, 10.3:1.1,29]';
end

% parameters
N = length(x_bold);
n = (0:N-1)';
Apslowerbound = 0.075; 
ii = N*0.05/fs;
f_rp = 0.0001; %(Hz)
i_rp = N*f_rp/fs;

% starting the timer
global Timer;
Timer = 0;
tic;

% =======================* RELAX algorithm *=======================
twopi_div_fs = 2*pi/fs;
twopi_div_N = 2*pi/N;
Nfreqs = 30;

% parameters for the zero-padding FFT
M_zpFFT = 4*N;

% iterations
alphas_hat = zeros(Nfreqs,1);
fs_hat = zeros(Nfreqs,1);
ms_hat = zeros(Nfreqs,1);
deltas_hat = zeros(Nfreqs,1);
components = zeros(N,Nfreqs);
has_been_deleted = zeros(Nfreqs,1);
xk = x_bold; 
for K = 1:Nfreqs
    C3a = -1000; C3b = -100;
    q = 0;
    while (abs((C3b-C3a)/C3a)>1e-3)  
        q = q + 1; k = K;
        while k>0
            % deleting other components
            for i = 1:K
                if i ~= k && ~has_been_deleted(i)
                    xk = xk - components(:,i);
                    has_been_deleted(i) = 1;
                end     
                if i == k && has_been_deleted(i)
                    xk = xk + components(:,i);
                    has_been_deleted(i) = 0;       
                end
            end

%             if q==1
                % estimating \m_k            
                FFT_result = fft(xk,M_zpFFT);
                J_FFT = FFT_result.*conj(FFT_result);
                [~,m_max] = max(J_FFT);
                ms_hat(k) = 0.25*(m_max(1)-1);
                deltas_hat(k) = 0;
%             end                     
            
            [fs_hat(k),FFT_result_fk] = Interpolation_Periodo(k);
            alphas_hat(k) = FFT_result_fk/N;

            components(:,k) = exp(1j*twopi_div_fs*fs_hat(k)*n)*alphas_hat(k);
            k = k - 1;
        end
        residual = xk - components(:,1);
        C3a = C3b; C3b = residual'*residual;
    end
    disp(q);
end

% finding out the carrier freq & the low-freq
Amplis_hat = abs(alphas_hat);
k_Emax = zeros(5,1);
deltaf_totalmax = 2;
Emax = 0;
f0_hat = -1; f1_hat = -1;
for p = 1:8
    for q = 1:18
        k_tmp = zeros(5,1);
        n_tmp = 0;
        deltaf_total = 0;
        for freq = f0_std(p)-2*f1_std(q):f1_std(q):f0_std(p)+2*f1_std(q)
            deltaf_min = fs;
            for k = 1:Nfreqs
                tmp = abs(freq-fs_hat(k));
                if deltaf_min > tmp
                    deltaf_min = tmp;
                    kmin = k;
                end
            end
            n_tmp = n_tmp + 1;
            k_tmp(n_tmp) = kmin;
            if n_tmp>1 && n_tmp<5
                deltaf_total = deltaf_total + deltaf_min;
            end
        end    
        E_tmp = Amplis_hat(k_tmp(1))^2 + Amplis_hat(k_tmp(2))^2 ...
            + Amplis_hat(k_tmp(3))^2 + Amplis_hat(k_tmp(4))^2 ...
            + Amplis_hat(k_tmp(5))^2; 
        if Amplis_hat(k_tmp(2)) > Apslowerbound &&...
                Amplis_hat(k_tmp(3)) > Apslowerbound &&...
                Amplis_hat(k_tmp(4)) > Apslowerbound &&...
                deltaf_total < deltaf_totalmax && E_tmp > Emax
            Emax = E_tmp;
            f0_hat = fs_hat(k_tmp(3));
            f1_hat = abs(fs_hat(k_tmp(2))-fs_hat(k_tmp(4)))/2;
            k_Emax = k_tmp(2:4);
        end
    end
end

% recording the time spent
Timer = Timer + toc;
fprintf('time consuming %fs\n',Timer);

%% ========================* functions *==========================
    function [f_hat,FFT_result_fhat] = Interpolation_Periodo(k)
        m = ms_hat(k)+deltas_hat(k);
        delta_now = m-floor(m);
        delta_add = 1;
        m = floor(m); 
        while (delta_add>i_rp)
            Xm = abs(CalcDTFT((m+delta_now-ii)*twopi_div_N));
            Xp = abs(CalcDTFT((m+delta_now+ii)*twopi_div_N));
            X0 = CalcDTFT((m+delta_now)*twopi_div_N);
            delta_add = N*atan((Xp-Xm)*sin(pi*ii/N)...
                /((Xp+Xm)*cos(pi*ii/N)-2*cos(pi*ii)*abs(X0)))/pi;
            delta_now = delta_now + delta_add;
        end
        f_hat = fs*(m+delta_now)/N;
        FFT_result_fhat = CalcDTFT((m+delta_now)*twopi_div_N);
    end
    function DTFT_result = CalcDTFT(omega)
        exp_vector = exp(-1j*omega*n);
        DTFT_result = sum(xk.*exp_vector);
    end
end