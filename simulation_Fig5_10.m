%% simulation of Fig 5.10 the comparison of accuracies of different algs under impulse interference
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

fs = 8000; Ts = 1/fs;
Tds = (0.1:0.1:1);
snr = 30; 
times = 30;

f0_mse = zeros(length(Tds),4);
f1_mse = zeros(length(Tds),4);
errorcount = zeros(length(Tds),4);

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
            
            t_impulse = (0:Ts:0.04)';
            t_traction = (0:Ts:(N-length(t_impulse)-401)*Ts)';

            A_im = 1111; a = 40; omega = 2*pi*300; 
            n_impulse = A_im*exp(-a*t_impulse).*sin(omega*t_impulse);

            n_traction = 222*sin(2*pi*50*t_traction);

            n_imp = [zeros(400,1);n_impulse;n_traction];
            
            for k = 1:times
                x = awgn(y,snr,'measured') + n_imp;
                [f0_hat,f1_hat] = YangFan2010(x,fs);
                f0_mse(i,1) = f0_mse(i,1) + (f0_hat-f0)^2;
                f1_mse(i,1) = f1_mse(i,1) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，YF出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,1) = errorcount(i,1) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured') + n_imp;
                [f0_hat,f1_hat] = Elias_Aboutanios2004(x,fs);
                f0_mse(i,2) = f0_mse(i,2) + (f0_hat-f0)^2;
                f1_mse(i,2) = f1_mse(i,2) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，EA出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,2) = errorcount(i,2) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured') + n_imp;
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                f0_mse(i,3) = f0_mse(i,3) + (f0_hat-f0)^2;
                f1_mse(i,3) = f1_mse(i,3) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，NLSM出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,3) = errorcount(i,3) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                x = awgn(y,snr,'measured') + n_imp;
                [f0_hat,f1_hat] = Rapid_RELAX_real_Dichotomous_Search(x,fs);
                f0_mse(i,4) = f0_mse(i,4) + (f0_hat-f0)^2;
                f1_mse(i,4) = f1_mse(i,4) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，rRELAX出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,4) = errorcount(i,4) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end           
        end
    end
    f0_mse(i,:) = f0_mse(i,:)/times/144; 
    f1_mse(i,:) = f1_mse(i,:)/times/144;
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_10 simulation results.mat'),...
'Tds','f0_mse','f1_mse','errorcount');
% costs about 116422 seconds in total