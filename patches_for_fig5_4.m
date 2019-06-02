%% patch for Fig 5.4 for updating results of rRELAX
% adding a if in line 19th of rRELAX's pseudo-code to avoid infinite loops 
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
Tds = [0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 (0.2:0.1:1)];
snr = 30; 
times = 30;

f0_mse_rRELAX = zeros(length(Tds),1);
f1_mse_rRELAX = zeros(length(Tds),1);
errorcount_rRELAX = zeros(length(Tds),1);

i = 0;
for Td = Tds
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
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，rRELAX出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount_rRELAX(i) = errorcount_rRELAX(i) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end            
        end
    end
    f0_mse_rRELAX(i,:) = f0_mse_rRELAX(i,:)/times/144; 
    f1_mse_rRELAX(i,:) = f1_mse_rRELAX(i,:)/times/144;
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_4 simulation results patch.mat'),...
'Tds','f0_mse_rRELAX','f1_mse_rRELAX','errorcount_rRELAX');
% costs about  seconds in total
%% patch for Fig 5.4 for updating errorcount of the two periodogram-based algs
% the alg ref[20] may need more sampling duration than other algs, so
% adding some Tds here to test how much sampling duration is needed by
% alg ref[20]
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
Tds = (0.21:0.01:0.29);
snr = 30; 
times = 30;

errorcount_periodo = zeros(length(Tds),2);

i = 0;
for Td = Tds
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
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，YF出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount_periodo(i,1) = errorcount_periodo(i,1) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured');
                [f0_hat,f1_hat] = Elias_Aboutanios2004(x,fs);
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，EA出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount_periodo(i,2) = errorcount_periodo(i,2) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end    
        end
    end
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_4 simulation results patch2.mat'),...
'Tds_periodo','errorcount_periodo');