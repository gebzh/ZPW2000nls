%% simulation of Fig 5.7 the comparison of accuracies of different algs varying with SNR
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
Td = 0.3;
SNRs = (-20:1:-1); 
times = 30;

f0_mse = zeros(length(SNRs),4);
f1_mse = zeros(length(SNRs),4);
errorcount = zeros(length(SNRs),4);
crbappf0 = zeros(length(SNRs),1);
crbappf1 = zeros(length(SNRs),1);
crbsimf0 = zeros(length(SNRs),1);
crbsimf1 = zeros(length(SNRs),1);

i = 0;
for snr = SNRs
    i = i + 1;
    for p = 1:8
        f0 = f0_std(p);
        for q = 1:18
            f1 = f1_std(q);

            A = 1;
            y = A * Generate2000Signal(f0,f1,fs,Td); % noise-free signal 
            N = length(y);
            P_ZPW = 10*log10((A^2)/2); 
            P_agwn = P_ZPW - snr;
            sigmaSquare = 10^(P_agwn/10);
            
            for k = 1:times
                x = awgn(y,snr,'measured');
                [f0_hat,f1_hat] = YangFan2010(x,fs);
                f0_mse(i,1) = f0_mse(i,1) + (f0_hat-f0)^2;
                f1_mse(i,1) = f1_mse(i,1) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('SNR = %fdB£¨‘ÿ∆µ°¢µÕ∆µŒ™%f Hz°¢%f Hz ±£¨YF≥ˆœ÷“Î¬Î¥ÌŒÛ£°\n',...
                        snr,f0,f1);
                    errorcount(i,1) = errorcount(i,1) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured');
                [f0_hat,f1_hat] = Elias_Aboutanios2004(x,fs);
                f0_mse(i,2) = f0_mse(i,2) + (f0_hat-f0)^2;
                f1_mse(i,2) = f1_mse(i,2) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('snr = %fdB£¨‘ÿ∆µ°¢µÕ∆µŒ™%f Hz°¢%f Hz ±£¨EA≥ˆœ÷“Î¬Î¥ÌŒÛ£°\n',...
                        snr,f0,f1);
                    errorcount(i,2) = errorcount(i,2) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured');
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                f0_mse(i,3) = f0_mse(i,3) + (f0_hat-f0)^2;
                f1_mse(i,3) = f1_mse(i,3) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('snr = %fdB£¨‘ÿ∆µ°¢µÕ∆µŒ™%f Hz°¢%f Hz ±£¨NLSM≥ˆœ÷“Î¬Î¥ÌŒÛ£°\n',...
                        snr,f0,f1);
                    errorcount(i,3) = errorcount(i,3) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured');
                [f0_hat,f1_hat] = Rapid_RELAX_real_Dichotomous_Search(x,fs);
                f0_mse(i,4) = f0_mse(i,4) + (f0_hat-f0)^2;
                f1_mse(i,4) = f1_mse(i,4) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('snr = %fdB£¨‘ÿ∆µ°¢µÕ∆µŒ™%f Hz°¢%f Hz ±£¨rRELAX≥ˆœ÷“Î¬Î¥ÌŒÛ£°\n',...
                        snr,f0,f1);
                    errorcount(i,4) = errorcount(i,4) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            [tmp1, tmp2] = AppCRBCalc(y,fs,f0,f1,sigmaSquare);
            crbappf0(i) = crbappf0(i) + tmp1;
            crbappf1(i) = crbappf1(i) + tmp2;
            [tmp1, tmp2] = SimCRBCalc(A,fs,N,f1,sigmaSquare);
            crbsimf0(i) = crbsimf0(i) + tmp1;
            crbsimf1(i) = crbsimf1(i) + tmp2;
        end
    end
    f0_mse(i,:) = f0_mse(i,:)/times/144; 
    f1_mse(i,:) = f1_mse(i,:)/times/144;
    crbappf0(i) = crbappf0(i) / 144;
    crbappf1(i) = crbappf1(i) / 144;
    crbsimf0(i) = crbsimf0(i) / 144;
    crbsimf1(i) = crbsimf1(i) / 144;
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_7 part1 simulation results.mat'),...
'SNRs','f0_mse','f1_mse','errorcount','crbappf0','crbappf1','crbsimf0','crbsimf1');
% costs about  50h in total