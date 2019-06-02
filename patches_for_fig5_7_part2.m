%% patch for Fig 5.7 part2
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
SNRs = (11:1:40); 
times = 30;

f0_mse_rRELAX = zeros(length(SNRs),1);
f1_mse_rRELAX = zeros(length(SNRs),1);
errorcount_rRELAX = zeros(length(SNRs),1);

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
                [f0_hat,f1_hat] = Rapid_RELAX_real_Dichotomous_Search(x,fs);
                f0_mse_rRELAX(i) = f0_mse_rRELAX(i) + (f0_hat-f0)^2;
                f1_mse_rRELAX(i) = f1_mse_rRELAX(i) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('snr = %fdB，载频、低频为%f Hz、%f Hz时，rRELAX出现译码错误！\n',...
                        snr,f0,f1);
                    errorcount_rRELAX(i) = errorcount_rRELAX(i) + 1;
                    breakpointForSectionImplementation( x,fs );
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
        end
    end
    f0_mse_rRELAX(i) = f0_mse_rRELAX(i)/times/144; 
    f1_mse_rRELAX(i) = f1_mse_rRELAX(i)/times/144;
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_7 simulation results patch part2.mat'),...
'SNRs','f0_mse_rRELAX','f1_mse_rRELAX','errorcount_rRELAX');
