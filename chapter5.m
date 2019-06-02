%% simulation of Fig 5.1 the comparison of running times of different algs
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

global Timer;
TIMER_ALL = 0;

f0 = 2001.4; f1=29; fs = 8000; 
Tds = (0.1:0.1:1);
snr = 30; 
times = 1000;
i = 0;

runningTime = zeros(length(Tds),4);

for Td = Tds
    i = i + 1;
    yn = Generate2000Signal(f0,f1,fs,Td); 
            
    for k = 1:times
        [f0_hat,f1_hat] = YangFan2010(awgn(yn,snr,'measured'),fs);
        runningTime(i,1) = runningTime(i,1) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end

    for k = 1:times
        [f0_hat,f1_hat] = Elias_Aboutanios2004(awgn(yn,snr,'measured'),fs);
        runningTime(i,2) = runningTime(i,2) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end

    for k = 1:times
        [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(awgn(yn,snr,'measured'),fs);
        runningTime(i,3) = runningTime(i,3) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end
    
    for k = 1:times
        [f0_hat,f1_hat] = Rapid_RELAX_real_Dichotomous_Search(awgn(yn,snr,'measured'),fs);
        runningTime(i,4) = runningTime(i,4) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end
    
    runningTime(i,:) = runningTime(i,:) / times;
end
disp(TIMER_ALL); 
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_1 simulation results.mat'),...
'Tds','runningTime');
% costs about 2.8 hours in total
%% Fig 5.1 the comparison of running times of different algs
load(strcat(pwd,'/data_files/Fig 5_1 simulation results.mat'));

addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
plot(Tds, runningTime(:,1),'-*k','MarkerSize',3.5);
hold on;
plot(Tds, runningTime(:,2),'-ok','MarkerSize',3.5);
hold on;
plot(Tds, runningTime(:,3),'-^k','MarkerSize',3.5);
hold on;
plot(Tds, runningTime(:,4),'-dk','MarkerSize',3.5);
axis([0.1 1 0 max(runningTime(:,1))+0.02]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = (0.1:0.1:max(runningTime(:,1)));
axes.Position = [0.2 0.18 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel('\fontname{宋体}采样持续时间\fontname{Times new roman}/s','Fontsize',7.5);
ylabel('\fontname{宋体}平均运行时间\fontname{Times new roman}/s','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 237 with the tag '% for chapter 5' in columnlegend.m is uncommented 
%% simulation of Fig 5.2 comparison of the two kinds of CRB when the sampling duration varying
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
Tds = (0.1:0.1:1);

crbappf0 = zeros(length(Tds),1);
crbappf1 = zeros(length(Tds),1);
crbsimf0 = zeros(length(Tds),1);
crbsimf1 = zeros(length(Tds),1);

i = 0;
snr = 30; 
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
            
            [tmp1, tmp2] = AppCRBCalc(y,fs,f0,f1,sigmaSquare);
            crbappf0(i) = crbappf0(i) + tmp1;
            crbappf1(i) = crbappf1(i) + tmp2;
            [tmp1, tmp2] = SimCRBCalc(A,fs,N,f1,sigmaSquare);
            crbsimf0(i) = crbsimf0(i) + tmp1;
            crbsimf1(i) = crbsimf1(i) + tmp2;
        end
    end
    crbappf0(i) = crbappf0(i) / 144;
    crbappf1(i) = crbappf1(i) / 144;
    crbsimf0(i) = crbsimf0(i) / 144;
    crbsimf1(i) = crbsimf1(i) / 144;
end
save(strcat(pwd,'/data_files/Fig 5_2 simulation results.mat'),...
'Tds','crbappf0','crbappf1','crbsimf0','crbsimf1');
msgbox('Operation Completed');
%% Fig 5.2 comparison of the two kinds of CRB when the sampling duration varying
load(strcat(pwd,'/data_files/Fig 5_2 simulation results.mat'));
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(Tds, sqrt(crbappf0),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--k','MarkerSize',3.5);
axis([0.1 1 1e-4 1e-2]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-4 3e-4 1e-3 3e-3 1e-2];
axes.Position = [0.21 0.23 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}CRB/Hz','Fontsize',7.5);
columnlegend(2,{'$CRB_{\mathrm{appf0}}$','$CRB_{\mathrm{simf0}}$'},'Interpreter','latex',...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 237 with the tag '% for chapter 5' in columnlegend.m is uncommented 

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(Tds, sqrt(crbappf1),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--k','MarkerSize',3.5);
grid on; 
axis([0.1 1 3e-4 3e-2]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [3e-4 1e-3 3e-3 1e-2 3e-2 1e-1];
axes.Position = [0.21 0.23 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}CRB/Hz','Fontsize',7.5);
columnlegend(2,{'$CRB_{\mathrm{appf1}}$','$CRB_{\mathrm{simf1}}$'},'Interpreter','latex',...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 237 with the tag '% for chapter 5' in columnlegend.m is uncommented 
%% simulation of Fig 5.3 comparison of the two kinds of CRB when the SNR varying
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
SNRs = (-20:2:50); 

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
            
            [tmp1, tmp2] = AppCRBCalc(y,fs,f0,f1,sigmaSquare);
            crbappf0(i) = crbappf0(i) + tmp1;
            crbappf1(i) = crbappf1(i) + tmp2;
            [tmp1, tmp2] = SimCRBCalc(A,fs,N,f1,sigmaSquare);
            crbsimf0(i) = crbsimf0(i) + tmp1;
            crbsimf1(i) = crbsimf1(i) + tmp2;
        end
    end
    crbappf0(i) = crbappf0(i) / 144;
    crbappf1(i) = crbappf1(i) / 144;
    crbsimf0(i) = crbsimf0(i) / 144;
    crbsimf1(i) = crbsimf1(i) / 144;
end
disp(TIMER_ALL); 
save(strcat(pwd,'/data_files/Fig 5_3 simulation results.mat'),...
'SNRs','crbappf0','crbappf1','crbsimf0','crbsimf1');
msgbox('Operation Completed');
%% Fig 5.3 comparison of the two kinds of CRB when the SNR varying
load(strcat(pwd,'/data_files/Fig 5_3 simulation results.mat'));
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(SNRs, sqrt(crbappf0),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs, sqrt(crbsimf0),'--k','MarkerSize',3.5);
axis([-20 50 1e-4 1]);
axes = gca;
axes.XTick = (-20:10:50);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1];
axes.Position = [0.19 0.23 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{Times new roman}SNR/dB';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}CRB/Hz','Fontsize',7.5);
columnlegend(2,{'$CRB_{\mathrm{appf0}}$','$CRB_{\mathrm{simf0}}$'},'Interpreter','latex',...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 237 with the tag '% for chapter 5' in columnlegend.m is uncommented 

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(SNRs, sqrt(crbappf1),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs, sqrt(crbsimf1),'--k','MarkerSize',3.5);
grid on; 
axis([-20 50 1e-4 1]);
axes = gca;
axes.XTick = (-20:10:50);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10];
axes.Position = [0.19 0.23 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;
xlabel({'\fontname{Times new roman}SNR/dB';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}CRB/Hz','Fontsize',7.5);
columnlegend(2,{'$CRB_{\mathrm{appf1}}$','$CRB_{\mathrm{simf1}}$'},'Interpreter','latex',...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 237 with the tag '% for chapter 5' in columnlegend.m is uncommented 
%% Fig 5.4 the comparison of accuracies of different algs
load(strcat(pwd,'/data_files/Fig 5_4 simulation results.mat'));
load(strcat(pwd,'/data_files/Fig 5_4 simulation results patch.mat'));
f0_mse(:,4) = f0_mse_rRELAX; f1_mse(:,4) = f1_mse_rRELAX;
errorcount(:,4) = errorcount_rRELAX;

addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 6.3]);
% for the legend
semilogy(Tds, sqrt(f0_mse(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf0),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--k','MarkerSize',3.5);

% erasing
hold on;
semilogy(Tds, sqrt(f0_mse(:,1)),'-*w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-ow','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-^w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-dw','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf0),'-w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--w','MarkerSize',3.5);

% real plot
hold on; temp = sqrt(f0_mse(:,1))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'*k','MarkerSize',3.5);
hold on; temp = sqrt(f0_mse(:,2))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'ok','MarkerSize',3.5);
hold on; temp = sqrt(f0_mse(:,3))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'^k','MarkerSize',3.5);
hold on; temp = sqrt(f0_mse(:,4))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,1)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf0),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf0),'--k','MarkerSize',3.5);

axis([0.09 1 1e-4 max(sqrt(f0_mse(:,4)))]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.21 0.78 0.62];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX',...
    '{\itCRB}_{\fontname{Times new roman}appf0}','{\itCRB}_{\fontname{Times new roman}simf0}'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 238 with the tag '% for chapter 5 plus' in columnlegend.m is uncommented 
box on;
figure; 
set(gcf,'unit','centimeters','position',[10 5 6 6.3]);
% for the legend
semilogy(Tds, sqrt(f1_mse(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf1),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--k','MarkerSize',3.5);

% erasing
hold on;
semilogy(Tds, sqrt(f1_mse(:,1)),'-*w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-ow','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-^w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-dw','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf1),'-w','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--w','MarkerSize',3.5);

% real plot
hold on; temp = sqrt(f1_mse(:,1))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'*k','MarkerSize',3.5);
hold on; temp = sqrt(f1_mse(:,2))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'ok','MarkerSize',3.5);
hold on; temp = sqrt(f1_mse(:,3))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'^k','MarkerSize',3.5);
hold on; temp = sqrt(f1_mse(:,4))';
semilogy([Tds(2) Tds(12:end)], [temp(2) temp(12:end)],'dk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,1)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbappf1),'-k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(crbsimf1),'--k','MarkerSize',3.5);

axis([0.09 1 1e-4 max(sqrt(f1_mse(:,2)))]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.21 0.78 0.62];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX',...
    '{\itCRB}_{\fontname{Times new roman}appf1}','{\itCRB}_{\fontname{Times new roman}simf1}'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
box on;
%% Fig 5.5 the comparison of amounts of decoding errors of different algs
load(strcat(pwd,'/data_files/Fig 5_4 simulation results patch2.mat'));
load(strcat(pwd,'/data_files/Fig 5_4 simulation results.mat'));
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
temp1 = [Tds(1:12) Tds_periodo Tds(13:end)];
tmp = errorcount(:,1); temp2 = [tmp(1:12); errorcount_periodo(:,1);tmp(13:end)];
temp2 = temp2 + 1e-18;
semilogy(temp1, temp2,'-*k','MarkerSize',3.5);
hold on;
tmp = errorcount(:,2); temp2 = [tmp(1:12); errorcount_periodo(:,2);tmp(13:end)];
temp2 = temp2 + 1e-18;
semilogy(temp1, temp2,'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, errorcount(:,3)+1e-18,'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, errorcount(:,4)+1e-18,'-dk','MarkerSize',3.5);
axis([0.09 0.3 1e-1 max(max(errorcount))]);
axes = gca;
axes.XTick = [0.11 0.15 0.19 0.25 0.3];
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.19 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s'},'Fontsize',7.5);
ylabel('\fontname{宋体}译码错误次数','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');

%% Fig 5.6 the influences of the multiplicative interference and the deviations
load(strcat(pwd,'/data_files/Fig 5_6 simulation results.mat'));

addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6.5 6.3]);
semilogy(Tds_multipliedeviated, sqrt(f0_mse_multipliedeviated(:,1)),'-vk','MarkerSize',3.5);
hold on;
semilogy(Tds_multipliedeviated, sqrt(f0_mse_multipliedeviated(:,2)),'->k','MarkerSize',3.5);
hold on;
semilogy(Tds_multipliedeviated, sqrt(f0_mse_multipliedeviated(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds_multipliedeviated, sqrt(f0_mse_multipliedeviated(:,4)),'-<k','MarkerSize',3.5);
axis([0.1 1 1e-3 1e-1]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-3 3e-3 1e-2 3e-2 1e-1];
axes.Position = [0.2 0.22 0.78 0.58];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
legend({'\fontname{宋体}乘性干扰、频率参数存在偏移',...
    '\fontname{宋体}无乘性干扰、频率参数存在偏移',...
    '\fontname{宋体}无乘性干扰、频率参数未偏移',...
    '\fontname{宋体}乘性干扰、频率参数未偏移'}...
    ,'Fontsize',7.5);
legend boxoff;
% hint: note that the line 236 with the tag '% for chapter 4' in columnlegend.m is uncommented 

figure; 
set(gcf,'unit','centimeters','position',[10 5 6.5 6.3]);
semilogy(Tds_multipliedeviated, sqrt(f1_mse_multipliedeviated(:,1)),'-vk','MarkerSize',3.5);
hold on;
semilogy(Tds_multipliedeviated, sqrt(f1_mse_multipliedeviated(:,2)),'->k','MarkerSize',3.5);
hold on;
semilogy(Tds_multipliedeviated, sqrt(f1_mse_multipliedeviated(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds_multipliedeviated, sqrt(f1_mse_multipliedeviated(:,4)),'-<k','MarkerSize',3.5);
grid on; 
axis([0.1 1 3e-4 3e-1]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1];
axes.Position = [0.2 0.22 0.78 0.58];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
legend({'\fontname{宋体}乘性干扰、频率参数存在偏移',...
    '\fontname{宋体}无乘性干扰、频率参数存在偏移',...
    '\fontname{宋体}无乘性干扰、频率参数未偏移',...
    '\fontname{宋体}乘性干扰、频率参数未偏移'}...
    ,'Fontsize',7.5);
legend boxoff;
%% Fig 5.7 the comparison of accuracies of different algs varying with SNR
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

SNRs_all = zeros(61,1);
f0_mse_all =  zeros(61,4); f1_mse_all =  zeros(61,4);
crbappf0_all = zeros(61,1); crbappf1_all = zeros(61,1);
crbsimf0_all = zeros(61,1); crbsimf1_all = zeros(61,1);

head = 1;
load(strcat(pwd,'/data_files/Fig 5_7 part1 simulation results.mat'));
SNRs_all(head:head+length(SNRs)-1) = SNRs;
f0_mse_all(head:head+length(SNRs)-1,:) = f0_mse;
f1_mse_all(head:head+length(SNRs)-1,:) = f1_mse;
crbappf0_all(head:head+length(SNRs)-1) = crbappf0;
crbappf1_all(head:head+length(SNRs)-1) = crbappf1;
crbsimf0_all(head:head+length(SNRs)-1) = crbsimf0;
crbsimf1_all(head:head+length(SNRs)-1) = crbsimf1;

head = head + length(SNRs);
load(strcat(pwd,'/data_files/Fig 5_7 part2 simulation results.mat'));
SNRs_all(head:head+length(SNRs)-1) = SNRs;
f0_mse_all(head:head+length(SNRs)-1,:) = f0_mse;
f1_mse_all(head:head+length(SNRs)-1,:) = f1_mse;
crbappf0_all(head:head+length(SNRs)-1) = crbappf0;
crbappf1_all(head:head+length(SNRs)-1) = crbappf1;
crbsimf0_all(head:head+length(SNRs)-1) = crbsimf0;
crbsimf1_all(head:head+length(SNRs)-1) = crbsimf1;
head = head + length(SNRs);
load(strcat(pwd,'/data_files/Fig 5_7 part3 simulation results.mat'));
SNRs_all(head:head+length(SNRs)-1) = SNRs;
f0_mse_all(head:head+length(SNRs)-1,:) = f0_mse;
f1_mse_all(head:head+length(SNRs)-1,:) = f1_mse;
crbappf0_all(head:head+length(SNRs)-1) = crbappf0;
crbappf1_all(head:head+length(SNRs)-1) = crbappf1;
crbsimf0_all(head:head+length(SNRs)-1) = crbsimf0;
crbsimf1_all(head:head+length(SNRs)-1) = crbsimf1;

load(strcat(pwd,'/data_files/Fig 5_7 simulation results patch part1.mat'));
f0_mse_all(1:31,4) = f0_mse_rRELAX; f1_mse_all(1:31,4) = f1_mse_rRELAX; 
load(strcat(pwd,'/data_files/Fig 5_7 simulation results patch part2.mat'));
f0_mse_all(32:end,4) = f0_mse_rRELAX; f1_mse_all(32:end,4) = f1_mse_rRELAX; 

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 6.3]);
% for the legend
semilogy(SNRs_all, sqrt(f0_mse_all(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,4)),'-dk','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbappf0_all),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbsimf0_all),'--k','MarkerSize',3.5);

% erasing
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,1)),'-*w','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,2)),'-ow','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,3)),'-^w','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,4)),'-dw','MarkerSize',3.5);

% real plot
temp1 = SNRs_all(1:4:end); 
hold on; tmp = sqrt(f0_mse_all(:,1)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'*k','MarkerSize',3.5);
hold on; tmp = sqrt(f0_mse_all(:,2)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'ok','MarkerSize',3.5);
hold on; tmp = sqrt(f0_mse_all(:,3)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'^k','MarkerSize',3.5);
hold on; tmp = sqrt(f0_mse_all(:,4)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'dk','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,1)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,2)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,3)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f0_mse_all(:,4)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbappf0_all),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbsimf0_all),'--k','MarkerSize',3.5);

axis([-20 40 1e-4 1e4]);
axes = gca;
axes.XTick = (-20:10:40);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
axes.Position = [0.18 0.21 0.78 0.6];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{Times new roman}SNR/dB';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX',...
    '{\itCRB}_{\fontname{Times new roman}appf0}','{\itCRB}_{\fontname{Times new roman}simf0}'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 238 with the tag '% for chapter 5 plus' in columnlegend.m is uncommented 
box on;
figure; 
set(gcf,'unit','centimeters','position',[10 5 6 6.3]);
% for the legend
semilogy(SNRs_all, sqrt(f1_mse_all(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,4)),'-dk','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbappf1_all),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbsimf1_all),'--k','MarkerSize',3.5);
% erasing
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,1)),'-*w','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,2)),'-ow','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,3)),'-^w','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,4)),'-dw','MarkerSize',3.5);

% real plot
temp1 = SNRs_all(1:4:end); 
hold on; tmp = sqrt(f1_mse_all(:,1)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'*k','MarkerSize',3.5);
hold on; tmp = sqrt(f1_mse_all(:,2)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'ok','MarkerSize',3.5);
hold on; tmp = sqrt(f1_mse_all(:,3)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'^k','MarkerSize',3.5);
hold on; tmp = sqrt(f1_mse_all(:,4)); temp2 = tmp(1:4:end);
semilogy(temp1,temp2,'dk','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,1)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,2)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,3)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(f1_mse_all(:,4)),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbappf1_all),'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, sqrt(crbsimf1_all),'--k','MarkerSize',3.5);

axis([-20 40 1e-4 1e2]);
axes = gca;
axes.XTick = (-20:10:40);
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 10 100 1000];
axes.Position = [0.18 0.21 0.78 0.6];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{Times new roman}SNR/dB';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX',...
    '{\itCRB}_{\fontname{Times new roman}appf1}','{\itCRB}_{\fontname{Times new roman}simf1}'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
box on;
%% Fig 5.8 the comparison of amounts of decoding errors of different algs
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

SNRs_all = zeros(61,1);
errorcount_all =  zeros(61,4);

head = 1;
load(strcat(pwd,'/data_files/Fig 5_7 part1 simulation results.mat'));
SNRs_all(head:head+length(SNRs)-1) = SNRs;
errorcount_all(head:head+length(SNRs)-1,:) = errorcount;

head = head + length(SNRs);
load(strcat(pwd,'/data_files/Fig 5_7 part2 simulation results.mat'));
SNRs_all(head:head+length(SNRs)-1) = SNRs;
errorcount_all(head:head+length(SNRs)-1,:) = errorcount;

head = head + length(SNRs);
load(strcat(pwd,'/data_files/Fig 5_7 part3 simulation results.mat'));
SNRs_all(head:head+length(SNRs)-1) = SNRs;
errorcount_all(head:head+length(SNRs)-1,:) = errorcount;

load(strcat(pwd,'/data_files/Fig 5_7 simulation results patch part1.mat'));
errorcount_all(1:31,4) = errorcount_rRELAX; 
load(strcat(pwd,'/data_files/Fig 5_7 simulation results patch part2.mat'));
errorcount_all(32:end,4) = errorcount_rRELAX; 

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(SNRs_all, errorcount_all(:,1)+1e-18,'-*k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,2)+1e-18,'-ok','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,3)+1e-18,'-^r','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,4)+1e-18,'-dr','MarkerSize',3.5);
% erasing
hold on;
semilogy(SNRs_all, errorcount_all(:,1)+1e-18,'-*w','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,2)+1e-18,'-ow','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,3)+1e-18,'-^w','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,4)+1e-18,'-dw','MarkerSize',3.5);
% real plot
temp1 = SNRs_all(2:2:end); 
hold on;
tmp = errorcount_all(:,1)+1e-18; temp2 = tmp(2:2:end);
semilogy(temp1, temp2,'*k','MarkerSize',3.5);
hold on;
tmp = errorcount_all(:,2)+1e-18; temp2 = tmp(2:2:end);
semilogy(temp1, temp2,'ok','MarkerSize',3.5);
hold on;
tmp = errorcount_all(:,3)+1e-18; temp2 = tmp(2:2:end);
semilogy(temp1, temp2,'^r','MarkerSize',3.5);
hold on;
tmp = errorcount_all(:,4)+1e-18; temp2 = tmp(2:2:end);
semilogy(temp1, temp2,'dr','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,1)+1e-18,'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,2)+1e-18,'-k','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,3)+1e-18,'-r','MarkerSize',3.5);
hold on;
semilogy(SNRs_all, errorcount_all(:,4)+1e-18,'-r','MarkerSize',3.5);
axis([-20 5 1e-1 1e4]);
axes = gca;
axes.XTick = [-20 -15 -10 -2 3 5];
axes.YTick = [1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];
axes.Position = [0.18 0.19 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{Times new roman}SNR/dB'},'Fontsize',7.5);
ylabel('\fontname{宋体}译码错误次数','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
%% Fig 5.9 the results of traction harmonic immunity of NLSM
%(a)
load(strcat(pwd,'/data_files/Fig 5_9 simulation results_Td08.mat'));
% load(strcat(pwd,'/data_files/Fig 5_9 simulation results_Td1 rRELAX.mat'));

f0 = 2001.4;

figure;
set(gcf,'unit','centimeters','position',[10 5 6.8 5.8]);
tmp = errorcount; size1 = size(tmp,1);
temp = zeros(size(tmp));
for i = 1:size1
    temp(i,:) = tmp(size1-i+1,:);
end
img = imagesc(temp);
xticklabels = f0-40:20:f0+40;
xticks = linspace(1, size(temp,2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
set(gca,'Position',[0.16 0.32 0.69 0.63],'Linewidth',0.283,...
    'XTickLabelRotation',45,'Fontname','Times new roman','Fontsize',7.5);
xlabel({'{\it\fontname{Times new roman}f}_{\fontname{Times new roman}har}\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 采样持续时间为\fontname{Times new roman}0.8s'}...
    ,'Fontsize',7.5);

yticklabels = (0.2:1.7:7);
yticks = linspace(1, size(temp,1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)));
ylabel('$A_{\mathrm{har}}$/V','Interpreter','latex','Fontsize',7.5);

colormap(flipud(gray));
cb = colorbar;
set(cb,'position',[0.88 0.32 0.05 0.63]);

%(b)
load(strcat(pwd,'/data_files/Fig 5_9 simulation results_Td1.mat'));

figure;
set(gcf,'unit','centimeters','position',[10 5 6.8 5.8]);
tmp = errorcount; size1 = size(tmp,1);
temp = zeros(size(tmp));
for i = 1:size1
    temp(i,:) = tmp(size1-i+1,:);
end
img = imagesc(temp);
xticklabels = f0-40:20:f0+40;
xticks = linspace(1, size(temp,2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
set(gca,'Position',[0.16 0.32 0.69 0.63],'Linewidth',0.283,...
    'XTickLabelRotation',45,'Fontname','Times new roman','Fontsize',7.5);
xlabel({'{\it\fontname{Times new roman}f}_{\fontname{Times new roman}har}\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 采样持续时间为\fontname{Times new roman}1s'}...
    ,'Fontsize',7.5);

yticklabels = (0.2:1.7:7);
yticks = linspace(1, size(temp,1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)));
ylabel('$A_{\mathrm{har}}$/V','Interpreter','latex','Fontsize',7.5);

colormap(flipud(gray));
cb = colorbar;
set(cb,'position',[0.88 0.32 0.05 0.63]);
%% Fig 5.10 the comparison of accuracies of different algs under impulse interference
load(strcat(pwd,'/data_files/Fig 5_10 simulation results.mat'));
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(Tds, sqrt(f0_mse(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,4)),'-dk','MarkerSize',3.5);
axis([0.1 1 1e-3 1e3]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-3 1e-2 1e-1 1 1e1 1e2 1e3];
axes.Position = [0.2 0.22 0.78 0.64];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
semilogy(Tds, sqrt(f1_mse(:,1)),'-*k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-ok','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,4)),'-dk','MarkerSize',3.5);
grid on; 
axis([0.1 1 1e-3 3e1]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1 3 1e1 3e1];
axes.Position = [0.2 0.22 0.78 0.64];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法',...
    '\fontname{Times new roman}NLSM','\fontname{Times new roman}rRELAX'},...
    'Location','northoutside','Fontsize',7.5,'boxoff');
%% Fig 5.11 the periodogram of a ZPW-2000 signal along with a impulse interference 
load(strcat(pwd,'/data_files/Fig 5_11 signal.mat'));

N = length(x);

periodo = abs(fft(x,20*N)).^2/N; 

reso_periodo = fs/(20*N);
f_right_margin = 1800;

kl = 0;
kr = round(f_right_margin/reso_periodo);
horizontaxis = reso_periodo * (kl:kr)';
verticaxis = 10*log10(periodo(kl+1:kr+1)/max(periodo(kl+1:kr+1)));

figure;
set(gcf,'unit','centimeters','position',[10 5 10 6.5]);
plot(horizontaxis,verticaxis,'k');
axis([0 f_right_margin min(verticaxis) 0]);
axes = gca;
axes.XTick = [50 (300:200:1700) ];
axes.YTick = (-100:10:0);
axes.Position = [0.13 0.23 0.82 0.75];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.XTickLabelRotation = 45;
axes.GridAlpha = 0.2;
box on;
xlabel({'\fontname{Times new roman}{\itf}/Hz'},'Fontsize',7.5);
ylabel('$10log_{10}[W_{\mathrm{x4}}(f)/max\{W_{\mathrm{x4}}(f)\}]/\mathrm{dB}$',...
    'Interpreter','latex','Fontsize',7.5);
grid on;
%% Fig 5.12 the periodogram of a windowed signal
load(strcat(pwd,'/data_files/Fig 5_11 signal.mat'));

N = length(x);
x = BlackmanWindow(x);

periodo = abs(fft(x,20*N)).^2/N; 

reso_periodo = fs/(20*N);
f_right_margin = 1800;

kl = 0;
kr = round(f_right_margin/reso_periodo);
horizontaxis = reso_periodo * (kl:kr)';
verticaxis = 10*log10(periodo(kl+1:kr+1)/max(periodo(kl+1:kr+1)));

figure;
set(gcf,'unit','centimeters','position',[10 5 10 6.5]);
plot(horizontaxis,verticaxis,'k');
axis([0 f_right_margin min(verticaxis) 0]);
axes = gca;
axes.XTick = [50 (300:200:1700) ];
axes.YTick = (-100:10:0);
axes.Position = [0.13 0.23 0.82 0.75];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.XTickLabelRotation = 45;
axes.GridAlpha = 0.2;
box on;
xlabel({'\fontname{Times new roman}{\itf}/Hz'},'Fontsize',7.5);
ylabel('$10log_{10}[W_{\mathrm{x5}}(f)/max\{W_{\mathrm{x5}}(f)\}]/\mathrm{dB}$',...
    'Interpreter','latex','Fontsize',7.5);
grid on;
%% Fig 5.13 
load(strcat(pwd,'/data_files/Fig 5_13 simulation results part2.mat'));
errorcount_all = zeros(size(errorcount,1),3);
errorcount_nosa_all = zeros(size(errorcount,1),3);
errorcount_all(:,3) = errorcount(:,1);
errorcount_nosa_all(:,3) = errorcount_nosa(:,1);

load(strcat(pwd,'/data_files/Fig 5_13 simulation results.mat'));
errorcount_all(:,1:2) = errorcount;
errorcount_nosa_all(:,1:2) = errorcount_nosa;


figure;
set(gcf,'unit','centimeters','position',[10 5 11 6.5]);
h = bar(Aadjs,[errorcount_all(:,1) errorcount_nosa_all(:,1)...
    errorcount_all(:,2) errorcount_nosa_all(:,2)...
    errorcount_all(:,3) errorcount_nosa_all(:,3)]);
axis([0.35 0.95 0 7000]);
set(gca,'XTick',(0.4:0.1:0.9));
set(gca,'YTick',(0:1000:7000));
set(gca,'GridAlpha',0.2);
set(gca,'Position',[0.13 0.15 0.85 0.82]);
colormap(summer(6));
grid on;
xlabel('$\overleftarrow{A}$/V','Interpreter','latex','Fontsize',7.5);
ylabel('\fontname{宋体}译码错误次数','Fontsize',7.5);
legend(h,{'\fontname{Times new roman}NLSM\fontname{宋体}总译码错误次数',...
    '\fontname{宋体}非同载频邻线干扰下\fontname{Times new roman}NLSM\fontname{宋体}译码错误次数',...
    '\fontname{宋体}文献\fontname{Times new roman}[19]\fontname{宋体}算法总译码错误次数',...
    '\fontname{宋体}非同载频邻线干扰下文献\fontname{Times new roman}[19]\fontname{宋体}算法译码错误次数 ',...
    '\fontname{宋体}文献\fontname{Times new roman}[20]\fontname{宋体}算法总译码错误次数',...
    '\fontname{宋体}非同载频邻线干扰下文献\fontname{Times new roman}[20]\fontname{宋体}算法译码错误次数 '},...
    'Fontsize',7.5);