%% simulation of Fig 5.6 the influences of the multiplicative interference and the deviations
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

global Timer;
global f0_std;
global f1_std;
TIMER_ALL = 0;

% loading data
if isempty(f0_std) || isempty(f1_std) 
    f0_std = [1701.4; 1698.7; 2001.4; 1998.7; 2301.4; 2298.7; 2601.4; 2598.7];
    f1_std = (10.3:1.1:29)';
end

fs = 8000;
Tds_multipliedeviated = (0.1:0.1:1);
snr = 30; 
times = 30;

f0_mse_multipliedeviated = zeros(length(Tds_multipliedeviated),4);
f1_mse_multipliedeviated = zeros(length(Tds_multipliedeviated),4);

i = 0;
for Td = Tds_multipliedeviated
    i = i + 1;
    for p = 1:8
        f0 = f0_std(p)+0.28;
        for q = 1:18
            f1 = f1_std(q)+0.03;

            y = Generate2000Signal(f0,f1,fs,Td); % noise-free signal 
            N = length(y);
            n = (0:N-1)';
            kn = 0.4051*n/fs + 1;     
            
            for k = 1:times
                x = awgn(kn.*Generate2000Signal(f0,f1,fs,Td),snr,'measured');
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                f0_mse_multipliedeviated(i,1) = f0_mse_multipliedeviated(i,1) + (f0_hat-f0)^2;
                f1_mse_multipliedeviated(i,1) = f1_mse_multipliedeviated(i,1) + (f1_hat-f1)^2;
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(Generate2000Signal(f0,f1,fs,Td),snr,'measured');
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                f0_mse_multipliedeviated(i,2) = f0_mse_multipliedeviated(i,2) + (f0_hat-f0)^2;
                f1_mse_multipliedeviated(i,2) = f1_mse_multipliedeviated(i,2) + (f1_hat-f1)^2;
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(Generate2000Signal(f0-0.28,f1-0.03,fs,Td),snr,'measured');
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                f0_mse_multipliedeviated(i,3) = f0_mse_multipliedeviated(i,3) + (f0_hat-(f0-0.28))^2;
                f1_mse_multipliedeviated(i,3) = f1_mse_multipliedeviated(i,3) + (f1_hat-(f1-0.03))^2;
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(kn.*Generate2000Signal(f0-0.28,f1-0.03,fs,Td),snr,'measured');
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                f0_mse_multipliedeviated(i,4) = f0_mse_multipliedeviated(i,4) + (f0_hat-(f0-0.28))^2;
                f1_mse_multipliedeviated(i,4) = f1_mse_multipliedeviated(i,4) + (f1_hat-(f1-0.03))^2;
                TIMER_ALL = TIMER_ALL + Timer;
            end
        end
    end
    f0_mse_multipliedeviated(i,:) = f0_mse_multipliedeviated(i,:)/times/144; 
    f1_mse_multipliedeviated(i,:) = f1_mse_multipliedeviated(i,:)/times/144;
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_6 simulation results.mat'),...
'Tds_multipliedeviated','f0_mse_multipliedeviated','f1_mse_multipliedeviated');
