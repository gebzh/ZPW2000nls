%% simulation of Fig 5.9 the testing & results of traction harmonic immunity
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

global Timer;
TIMER_ALL = 0;

f0 = 2001.4; f1 = 29; fs = 8000; Td = 1;
s = Generate2000Signal(f0,f1,fs,Td);
snr = 5;
times = 30;

N = length(s);
n = (0:N-1)';

fhars = f0-40:1:f0+40;
Ahars = (0.2:0.2:7);
errorcount = zeros(length(Ahars),length(fhars));
f0_mse = zeros(length(Ahars),length(fhars));
f1_mse = zeros(length(Ahars),length(fhars));

i = 0; j = 0;
for Ahar = Ahars
    i = i + 1; j = 0;
    for fhar = fhars
        j = j + 1;
        for k = 1:times
            phihar = 2*pi*rand(1);
            n_har = Ahar * cos(2*pi*fhar*n/fs + phihar);
            x = awgn(s,snr,'measured') + n_har;
%             [f0_hat,f1_hat] = Rapid_RELAX_real_Dichotomous_Search(x,fs);
            [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
            TIMER_ALL = TIMER_ALL + Timer;
            f0_mse(i,j) = f0_mse(i,j) + (f0_hat-f0)^2;
            f1_mse(i,j) = f1_mse(i,j) + (f1_hat-f1)^2;
            if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                errorcount(i,j) = errorcount(i,j) + 1;
                fprintf('Ahar=%f A fhar=%f Hz时，出现译码错误！\n',Ahar,fhar);
            end
        end
        f0_mse(i,j) = f0_mse(i,j) / times;
        f1_mse(i,j) = f1_mse(i,j) / times;
    end
end
disp(TIMER_ALL);
msgbox('Operation Completed');
% save(strcat(pwd,'/data_files/Fig 5_9 simulation results_Td1 rRELAX.mat'),...
% 'Ahars','fhars','f0_mse','f1_mse','errorcount');
save(strcat(pwd,'/data_files/Fig 5_9 simulation results_Td1.mat'),...
    'Ahars','fhars','f0_mse','f1_mse','errorcount');
% costs about 116422 seconds in total

% for rRELAX, Ahar=1.400000 A fhar=2002.400000 Hz时，出现译码错误！
%% simulation of Fig 5.9 the testing & results of traction harmonic immunity 
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

global Timer;
TIMER_ALL = 0;

f0 = 2001.4; f1 = 29; fs = 8000; Td = 0.8;
s = Generate2000Signal(f0,f1,fs,Td);
snr = 5;
times = 30;

N = length(s);
n = (0:N-1)';

fhars = f0-40:1:f0+40;
Ahars = (0.2:0.2:7);
errorcount = zeros(length(Ahars),length(fhars));
f0_mse = zeros(length(Ahars),length(fhars));
f1_mse = zeros(length(Ahars),length(fhars));

i = 0; j = 0;
for Ahar = Ahars
    i = i + 1; j = 0;
    for fhar = fhars
        j = j + 1;
        for k = 1:times
            phihar = 2*pi*rand(1);
            n_har = Ahar * cos(2*pi*fhar*n/fs + phihar);
            x = awgn(s,snr,'measured') + n_har;
            [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
            TIMER_ALL = TIMER_ALL + Timer;
            f0_mse(i,j) = f0_mse(i,j) + (f0_hat-f0)^2;
            f1_mse(i,j) = f1_mse(i,j) + (f1_hat-f1)^2;
            if (abs(f0_hat-f0)>0.5 || abs(f1_hat-f1)>0.5)
                errorcount(i,j) = errorcount(i,j) + 1;
            end
        end
        f0_mse(i,j) = f0_mse(i,j) / times;
        f1_mse(i,j) = f1_mse(i,j) / times;
    end
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_9 simulation results_Td08.mat'),...
'Ahars','fhars','f0_mse','f1_mse','errorcount');
% costs about 116422 seconds in total