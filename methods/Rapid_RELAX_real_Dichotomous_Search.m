function [ f0_hat, f1_hat, fs_hat, Amplis_hat, k_Emax ] = Rapid_RELAX_real_Dichotomous_Search( x_bold,fs )
% demodulating the ZPW-2000 signal using the RELAX algorithm (J Li 1996)
% /*rapid version 1*/ the using the dichotomous search method to reduce
% the computational burden.
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
Nfreqs = 12;
N = length(x_bold);
n0 = -(N-1)/2;
n = (n0:n0+N-1)';
twopi_div_fs = 2*pi/fs;
fs_div_4N = fs/(4*N);
fs_div_twopi = fs/(2*pi);
Apslowerbound = 0.15;
Alowestcomponent = 0.01;
f_rp = 0.0001; 
omega_dsrp = f_rp*twopi_div_fs;

% starting the timer
global Timer;
Timer = 0;
tic;

% =======================* freqs detection process *=======================

% parameters for the zero-padding FFT
M_zpFFT = 4*N;
M_half = round(M_zpFFT/2);
omega_zpFFT = 2*pi*(0:M_half-1)'/M_zpFFT;
% complementary parameters for the real-model RELAX algorithm 
temp = sin(N*omega_zpFFT) ./ sin(omega_zpFFT);
v1_FFT = 2 ./ (N + temp);
v2_FFT = 2 ./ (N - temp);

alphas_hat = zeros(2,Nfreqs);
fs_hat = -fs*ones(Nfreqs,1);
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
            
            % estimating \f_k & \alpha_k
            if fs_hat(k)==-fs
                FFT_result = fft(xk,M_zpFFT);
                temp = exp(-1j*n0*omega_zpFFT).*FFT_result(1:M_half);
                CTx_FFT = real(temp); % M/2 dimmension vector
                STx_FFT = -imag(temp); 
                J_FFT = v1_FFT.*CTx_FFT.^2 + v2_FFT.*STx_FFT.^2;
                kl = round(1660/fs_div_4N); kl = 1;
                [~,m_max] = max(J_FFT(kl:end));
                fs_hat(k) = fs_div_4N*(m_max(1)+kl-2);
            end
            
            omega0 = fs_hat(k)*twopi_div_fs;
            omegaminus = (fs_hat(k)-fs_div_4N)*twopi_div_fs;
            omegaplus = (fs_hat(k)+fs_div_4N)*twopi_div_fs;
            if CalcCostFunc(omegaminus)>CalcCostFunc(omegaplus)
                omegamid = omega0 + pi/(4*N);
                [fs_hat(k),alphas_hat(:,k)] = ...
                    Dichotomous_Search_of_Max_CostFunc(omegaminus,omegamid);
            else
                omegamid = omega0 - pi/(4*N);
                [fs_hat(k),alphas_hat(:,k)] = ...
                    Dichotomous_Search_of_Max_CostFunc(omegamid,omegaplus);
            end
            components(:,k) = [cos(fs_hat(k)*twopi_div_fs*n),...
                sin(fs_hat(k)*twopi_div_fs*n)]*alphas_hat(:,k);
            k = k - 1;
        end
        residual = xk - components(:,1);
        C3a = C3b; C3b = residual'*residual;
    end
    if sqrt(alphas_hat(1,K)^2+alphas_hat(2,K)^2)<Alowestcomponent
        break;
    end
%     disp(q);
end

% =======================* decoding process *=======================
% finding out the carrier freq & the low-freq
Amplis_hat = sqrt(alphas_hat(1,:).^2 + alphas_hat(2,:).^2);
k_Emax = zeros(5,1);
Emax = 0;
f0_hat = -1; f1_hat = -1;
for p = 1:8
    for q = 1:18
        k_tmp = zeros(5,1);
        n_tmp = 0;
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
        end    
        E_tmp = (Amplis_hat(k_tmp(1))^2 + Amplis_hat(k_tmp(2))^2 ...
            + Amplis_hat(k_tmp(3))^2 + Amplis_hat(k_tmp(4))^2 ...
            + Amplis_hat(k_tmp(5))^2)/2; 
        if Amplis_hat(k_tmp(2)) > Apslowerbound &&...
                Amplis_hat(k_tmp(3)) > Apslowerbound &&...
                Amplis_hat(k_tmp(4)) > Apslowerbound &&...
                k_tmp(2)~=k_tmp(3) && k_tmp(3)~=k_tmp(4) && k_tmp(2)~=k_tmp(4)...
                && E_tmp > Emax
            Emax = E_tmp;
            f0_hat = fs_hat(k_tmp(3));
            f1_hat = abs(fs_hat(k_tmp(2))-fs_hat(k_tmp(4)))/2;
            k_Emax = k_tmp(2:4);
        end
    end
end

% recording the time spent
Timer = Timer + toc;
fprintf('time consuming rRELAX %fs\n',Timer);

%% ========================* functions *==========================
    function [f_hat,alpha_fhat] = Dichotomous_Search_of_Max_CostFunc(omegal,omegar)
        % freqs in this function are in rad/sample 
        omegam = 0.5*(omegal+omegar);
        delta = 0.5*(omegar-omegal);
        [Jl,~] = CalcCostFunc(omegal);
        [Jr,~] = CalcCostFunc(omegar);
        while (delta>omega_dsrp)
            [Jm,alpha_fhat] = CalcCostFunc(omegam); 
            delta=delta*0.5;
            if (Jr>Jl)
                Jl = Jm;
                omegam = omegam + delta;
            else
                Jr = Jm;
                omegam = omegam - delta;
            end
        end
        f_hat = fs_div_twopi*omegam;
    end
    function [J,alpha] = CalcCostFunc(omega)
        Z = [cos(omega*n) sin(omega*n)];
        eta = Z'*xk;
        v1 = 2/(N+sin(omega*N)/sin(omega));
        v2 = 2/(N-sin(omega*N)/sin(omega));
        alpha = [eta(1)*v1; eta(2)*v2];
        J = eta'*alpha;
    end
end