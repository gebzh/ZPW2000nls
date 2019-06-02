function [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_Based_Algorithm_wa( x_bold, fs )
% rapid RELAX method based ZPW-2000 signal demodulation alg
% according to ref[68]
%   inputs:
%   x_bold: the received signal;  fs: the sampling frequency
%   outputs:
%   fs_hat: the estimated frequencies of components
%   Amplis_hat: the estimated amplitudes of components
%   f0_hat: the carrier frequency est; f1_hat: the low frequency est
%   k_Emax: the numbers of the components belonging to the ZPW-2000 signal

N = length(x_bold);
n = (0:N-1)';
Apslowerbound = 0.075; 
% the lower bound of the amplitude of the carrier freq and the first-side
% freqs components

%========================* loading data *==========================
global f0_std;
global f1_std;

% definition of the standard carrier freqs, the standard low freqs and the
% modulation indexes 
if isempty(f0_std) || isempty(f1_std) 
    f0_std = [1701.4; 1698.7; 2001.4; 1998.7; 2301.4; 2298.7; 2601.4; 2598.7];
    f1_std = (10.3:1.1:29)';
end

% starting the timer
global Timer;
Timer = 0;
tic;

%========================* RELAX alg *==========================
Nfreqs = 24; % assuming there's L freq components of complex value
fs_hat = zeros(Nfreqs,1);
ms_hat = zeros(Nfreqs,1);
deltas_hat = zeros(Nfreqs,1); 
alphas_hat = zeros(Nfreqs,1);

q = 1;
x = x_bold;
L_now = Nfreqs;

while (q < 21)
    for i = 1:Nfreqs       
        % calculating \hat{m}_i
        if q == 1 
            temp = fft(x,N*4);
            periodogram = temp.*conj(temp); 
            [~,m_max] = max(periodogram);
            ms_hat(i) = 0.25 * (m_max(1) - 1);
            L_now = i;
        end
        
        % calculating leakages
        leakages_plus = 0; leakages_minus = 0; 
        for j = 1:L_now
            if j ~= i
                M_ji = ms_hat(j) - ms_hat(i);
                delta_ji = deltas_hat(j) - deltas_hat(i);
                leakages_plus = leakages_plus + ...
                    alphas_hat(j)*(1+exp(2j*pi*delta_ji))...
                    /(1-exp(2j*pi*(M_ji+delta_ji-0.5)/N))/N;
                leakages_minus = leakages_minus + ...
                    alphas_hat(j)*(1+exp(2j*pi*delta_ji))...
                    /(1-exp(2j*pi*(M_ji+delta_ji+0.5)/N))/N;
            end
        end
        
        % calculating leakage-free Fourier coefficient
        temp = exp(-1j*2*pi*(ms_hat(i)+deltas_hat(i)+0.5)*n/N);
        Xplus = sum(x_bold.*temp)/N;
        Splus = Xplus - leakages_plus;
        
        temp = exp(-1j*2*pi*(ms_hat(i)+deltas_hat(i)-0.5)*n/N);
        Xminus = sum(x_bold.*temp)/N;
        Sminus = Xminus - leakages_minus;  
        
        % updating the delta of ith component using AM interpolation
        h = real((Splus+Sminus)/(Splus-Sminus))/2;
        deltas_hat(i) = deltas_hat(i) + h;
        
        % updating the estimated amplitude of ith component 
        leakages_ampli = 0;
        for j = 1:L_now
            if j ~= i
                f_ji = ms_hat(j)+deltas_hat(j) - (ms_hat(i)+deltas_hat(i));  
                if abs(f_ji) > 1e-7
                    leakages_ampli = leakages_ampli + ...
                        alphas_hat(j)*(1-exp(2j*pi*f_ji))/(N*(1-exp(2j*pi*f_ji/N)));
                else
                    leakages_ampli = leakages_ampli - alphas_hat(j);  
                    alphas_hat(j) = 0;
                end
            end
        end
        temp = exp(-2j*pi*(ms_hat(i)+deltas_hat(i))*n/N);
        alphas_hat(i) = sum(temp.*x_bold)/N - leakages_ampli;
        
        % deleting the estimated component in the first iteration
        if q == 1
            temp = alphas_hat(i)*exp(2j*pi*(ms_hat(i)+deltas_hat(i))*n/N);
            x = x - temp;
        end
    end
    q = q + 1;
end

for i = 1:Nfreqs
    fs_hat(i) = fs*(ms_hat(i)+deltas_hat(i))/N;
end

% finding out the carrier freq & the low-freq
Amplis_hat = abs(alphas_hat);
k_Emax = zeros(3,1);
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

% record the time cost
Timer = Timer + toc;
fprintf('time consuming of the rapid RELAX %fs\n',Timer);

end
%========================* some notes *==========================
% 1.
% when the amplitudes of freqs occupying a large range of value, it means that, 
% the leakage of the larger-amplitude components left during the first iteration
% may higher than the amplitude of the lower-amplitude component and thus the 
% larger one would be regarded as 'many components' in the first iteration, and 
% the similar scenario also happens when the number of the components is not 
% a priori known, that the freqs of last few components produced by the first 
% iteration will be close and will be seemed as 'many components' at one freq 
% since the expected number of freqs is larger than the reality. Then along with
% the iterations, the freqs of the 'many components' converged together, the NaN 
% value would be produced cause f_ji would be zero in that case. 
% Hence, in order to avoid this situation, the alg may should be adjusted that
% the iteration process is adopted at first to confirm the convergence, and
% then the exptected number of the freq is able to increase.
% And an another choice to handle this situation is to merge 'many components'
% into one freq component during the implementation.
