%% Fig 4.1 cost function grid map of the ZPW-2000 signal
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f0 = 2001.4; f1=29; fs = 8000; Td = 0.3;
s = Generate2000Signal(f0,f1,fs,Td);

[~,~,~,~,~,~,~,~,~,omega0s, omega1s, Jmat] = NLSM_Based_Algorithm_2( s,fs,1 );

figure;
set(gcf,'unit','centimeters','position',[10 5 8 4.7]);
x = fs*omega1s/2/pi; y = fs*omega0s/2/pi; z = Jmat;
mesh(x,y,z);
axis([10 31 1994.4 2008.4 0 max(max(Jmat))]);
grid on;
set(gca,'LineWidth',0.283);
set(gca,'position',[0.2 0.2 0.7 0.8]);
set(get(gca,'xlabel'),'Rotation',14);
set(gca,'XTick',[10 15 20 25 29]);
set(get(gca,'ylabel'),'Rotation',-24);
set(gca,'Fontsize',7.5,'Fontname','Times new roman');
xlabel('\fontname{宋体}低频频率\fontname{Times new roman}/Hz','Fontsize',7.5); 
ylabel('\fontname{宋体}载频频率\fontname{Times new roman}/Hz','Fontsize',7.5); 
zlabel('\fontname{宋体}价值函数','Fontsize',7.5);
% set(gcf,'WindowStyle','normal');
%% simulation of Fig 4.2 the comparison of accuracies of different algs
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
snr = 30; 
times = 30;
i = 0;

f0_mse = zeros(length(Tds),3);
f1_mse = zeros(length(Tds),3);
errorcount = zeros(length(Tds),3);

for Td = Tds
    i = i + 1;
    for p = 1:8
        f0 = f0_std(p);
        for q = 1:18
            f1 = f1_std(q);

            yn = Generate2000Signal(f0,f1,fs,Td);  
            
            for k = 1:times
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_1(awgn(yn,snr,'measured'),fs);
                f0_mse(i,1) = f0_mse(i,1) + (f0_hat-f0)^2;
                f1_mse(i,1) = f1_mse(i,1) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.4 || abs(f1_hat-f1)>0.4)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，alg1出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,1) = errorcount(i,1) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_2(awgn(yn,snr,'measured'),fs);
                f0_mse(i,2) = f0_mse(i,2) + (f0_hat-f0)^2;
                f1_mse(i,2) = f1_mse(i,2) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.4 || abs(f1_hat-f1)>0.4)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，alg1出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,2) = errorcount(i,2) + 1;
                end
                TIMER_ALL = TIMER_ALL + Timer;
            end
            
            for k = 1:times
                [f0_hat,f1_hat] = NLSM_Based_Algorithm_3(awgn(yn,snr,'measured'),fs);
                f0_mse(i,3) = f0_mse(i,3) + (f0_hat-f0)^2;
                f1_mse(i,3) = f1_mse(i,3) + (f1_hat-f1)^2;
                if (abs(f0_hat-f0)>0.4 || abs(f1_hat-f1)>0.4)
                    fprintf('Td = %fs，载频、低频为%f Hz、%f Hz时，alg1出现译码错误！\n',...
                        Td,f0,f1);
                    errorcount(i,3) = errorcount(i,3) + 1;
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
% save(strcat(pwd,'/data_files/Fig 4_2 simulation results.mat'),...
% 'Tds','f0_mse','f1_mse','errorcount');
% costs about 61400 seconds in total
%% Fig 4.2 the comparison of accuracies of different algs
load(strcat(pwd,'/data_files/Fig 4_2 simulation results.mat'));

addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6.2 5.8]);
semilogy(Tds, sqrt(f0_mse(:,1)),'-sk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,2)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f0_mse(:,3)),'-xk','MarkerSize',3.5);
axis([0.1 1 0.003 1e-1]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [3e-3 1e-2 3e-2 1e-1];
axes.Position = [0.2 0.22 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 载频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}方案一','\fontname{宋体}方案二',...
    '\fontname{宋体}方案三'},'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 236 with the tag '% for chapter 4' in columnlegend.m is uncommented 

figure; 
set(gcf,'unit','centimeters','position',[10 5 6.2 5.8]);
semilogy(Tds, sqrt(f1_mse(:,1)),'-sk','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,2)),'-^k','MarkerSize',3.5);
hold on;
semilogy(Tds, sqrt(f1_mse(:,3)),'-xk','MarkerSize',3.5);
grid on; 
axis([0.1 1 0.0005 0.3]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = [1e-4 1e-3 1e-2 1e-1];
axes.Position = [0.2 0.22 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;
xlabel({'\fontname{宋体}采样持续时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 低频估计'},'Fontsize',7.5);
ylabel('\fontname{Times new roman}RMSE/Hz','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}方案一','\fontname{宋体}方案二',...
    '\fontname{宋体}方案三'},'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 236 with the tag '% for chapter 4' in columnlegend.m is uncommented 
%% simulation of Fig 4.3 the comparison of running times of different algs
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

f0 = 2001.4; f1=29; fs = 8000; 
Tds = (0.1:0.1:1);
snr = 30; 
times = 1000;
i = 0;

runningTime = zeros(length(Tds),3);

for Td = Tds
    i = i + 1;
    yn = Generate2000Signal(f0,f1,fs,Td);  
            
    for k = 1:times
        [f0_hat,f1_hat] = Proposed_Algorithm_1(awgn(yn,snr,'measured'),fs);
        runningTime(i,1) = runningTime(i,1) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end

    for k = 1:times
        [f0_hat,f1_hat] = Proposed_Algorithm_2(awgn(yn,snr,'measured'),fs);
        runningTime(i,2) = runningTime(i,2) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end

    for k = 1:times
        [f0_hat,f1_hat] = Proposed_Algorithm_3(awgn(yn,snr,'measured'),fs);
        runningTime(i,3) = runningTime(i,3) + Timer;
        TIMER_ALL = TIMER_ALL + Timer;
    end
    
    runningTime(i,:) = runningTime(i,:) / times;
end
disp(TIMER_ALL); 
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 4_3 simulation results.mat'),...
'Tds','runningTime');
% costs about 13745 seconds in total
%% Fig 4.3 the comparison of running times of different algs
load(strcat(pwd,'/data_files/Fig 4_3 simulation results.mat'));

addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'));

figure; 
set(gcf,'unit','centimeters','position',[10 5 6 5.8]);
plot(Tds, runningTime(:,1),'-sk','MarkerSize',3.5);
hold on;
plot(Tds, runningTime(:,2),'-^k','MarkerSize',3.5);
hold on;
plot(Tds, runningTime(:,3),'-xk','MarkerSize',3.5);
axis([0.1 1 0 max(runningTime(:,1))]);
axes = gca;
axes.XTick = (0.1:0.1:1);
axes.YTick = (0.1:0.1:max(runningTime(:,1)));
axes.Position = [0.2 0.21 0.78 0.65];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on; 
xlabel('\fontname{宋体}采样持续时间\fontname{Times new roman}/s','Fontsize',7.5);
ylabel('\fontname{宋体}平均运行时间\fontname{Times new roman}/s','Fontsize',7.5);
columnlegend(2,{'\fontname{宋体}方案一','\fontname{宋体}方案二',...
    '\fontname{宋体}方案三'},'Location','northoutside','Fontsize',7.5,'boxoff');
% hint: note that the line 236 with the tag '% for chapter 4' in columnlegend.m is uncommented 
%% Fig 4.4 the spectrums estimated by RELAX alg and periodogram
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f_c = [1701.4 1698.7 2001.4 1998.7 2301.4 2298.7 2601.4 2598.7];
f_l = 10.3:1.1:29;

pp(1) = 1; qq(1) = 11; pp(2) = 2; qq(2) = 1;
% pp(1) = 2; qq(1) = 18; pp(2) = 1; qq(2) = 3;

AA(1) = 1; AA(2) = 0.7;

Td = 0.3;

f0 = f_c(pp(1)); f1 = f_l(qq(1)); fs = 8000; 
s1 = AA(1)*Generate2000Signal(f0,f1,fs,Td);
snr = 30;
s1 = awgn(s1,snr,'measured');
f0 = f_c(pp(2)); f1 = f_l(qq(2)); fs = 8000; 
s2 = AA(2)*Generate2000Signal(f0,f1,fs,Td);
s = s1 + s2;

[ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = RELAX_alpha_complex_model(s,fs);
% [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_complex_Dichotomous_Search(s,fs);
% [f0_hat, f1_hat, fs_hat, Amplis_hat, k_Emax ] = Rapid_RELAX_real_Dichotomous_Search(s1,fs);
% Amplis_hat = Amplis_hat / 2;
N = length(s);
M = 100 * N;

periodogram = abs(fft(s,M))/N;

reso_FFT = fs / M;
kl = round(1660/reso_FFT);
kr = round(1740/reso_FFT);

figure;
set(gcf,'unit','centimeters','position',[10 5 10 6.5]);
% plot the standard spectrum
f = zeros(10,1);
A_true = zeros(10,1);
for i = 1:2
    f(1+5*(i-1)) = f_c(pp(i)) - 2 * f_l(qq(i));
    f(2+5*(i-1)) = f_c(pp(i)) - f_l(qq(i));
    f(3+5*(i-1)) = f_c(pp(i));
    f(4+5*(i-1)) = f_c(pp(i)) + f_l(qq(i));
    f(5+5*(i-1)) = f_c(pp(i)) + 2 * f_l(qq(i));
    m = 11/f_l(qq(i));
    A_true(1+5*(i-1)) = AA(i)*abs( 2*m*sin(m*pi/2)/((m^2-4)*pi) );
    A_true(2+5*(i-1)) = AA(i)*abs( 2*m*cos(m*pi/2)/((m^2-1)*pi) );
    A_true(3+5*(i-1)) = AA(i)*abs( 2*sin(m*pi/2)/(m*pi) );
    A_true(4+5*(i-1)) = A_true(2+5*(i-1));
    A_true(5+5*(i-1)) = A_true(1+5*(i-1));
end
hold on;
plot(f, A_true, '*k','MarkerSize',3.5);
% plot estimation results
hold on;
stem(fs_hat, 2*Amplis_hat, '-k','Marker','none','Linewidth',0.3);
% plot the periodogram
hold on;
plot(kl*reso_FFT:reso_FFT:kr*reso_FFT,2*periodogram(kl+1:kr+1),'-.k','Linewidth',0.3);
if k_Emax(1) > 0
    hold on;
    stem([fs_hat(k_Emax(1)),fs_hat(k_Emax(2)),fs_hat(k_Emax(3))],...
        2*[Amplis_hat(k_Emax(1)),Amplis_hat(k_Emax(2)),Amplis_hat(k_Emax(3))]...
        ,'k','Linewidth',1.5,'Marker','none');
end
grid on;
axis([1670 1730 0 1]);
axes = gca;
axes.XTick = sort([f_c(pp(1)) f_c(pp(1))-f_l(qq(1)) f_c(pp(1))+f_l(qq(1))...
    f_c(pp(2)) f_c(pp(2))-f_l(qq(2)) f_c(pp(2))+f_l(qq(2))...
    f_c(pp(2))-2*f_l(qq(2)) f_c(pp(2))+2*f_l(qq(2)) ]);
axes.YTick = (0:0.2:1);
axes.Position = [0.13 0.23 0.82 0.76];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.XTickLabelRotation = 45;
axes.GridAlpha = 0.2;
box on;
xlabel({'\fontname{宋体}频率\fontname{Times new roman}/Hz'},'Fontsize',7.5);
ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
legend({'\fontname{宋体}实际频谱',...
    '\fontname{Times new roman}RELAX\fontname{宋体}法所得频谱',...
    '\fontname{宋体}周期图'},'Fontsize',7.5,'Location','northoutside',...
    'Orientation','horizontal');
legend boxoff;
%% Fig 4.5 the spectrums estimated by ref 67&68 and periodogram
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f_c = [1701.4 1698.7 2001.4 1998.7 2301.4 2298.7 2601.4 2598.7];
f_l = 10.3:1.1:29;

pp(1) = 1; qq(1) = 11; pp(2) = 2; qq(2) = 1;
% pp(1) = 2; qq(1) = 18; pp(2) = 1; qq(2) = 3;

AA(1) = 1; AA(2) = 0.7;

Td = 0.3;

f0 = f_c(pp(1)); f1 = f_l(qq(1)); fs = 8000; 
s1 = AA(1)*Generate2000Signal(f0,f1,fs,Td);
snr = 30;
s1 = awgn(s1,snr,'measured');
f0 = f_c(pp(2)); f1 = f_l(qq(2)); fs = 8000; 
s2 = AA(2)*Generate2000Signal(f0,f1,fs,Td);
s = s1 + s2;

% [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_alpha_Dichotomous_Search(s,fs);
% [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_real_Dichotomous_Search(s,fs);
[ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_Based_Algorithm_wa(s,fs);
% [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_alpha(s,fs);

N = length(s);
M = 100 * N;

figure;
set(gcf,'unit','centimeters','position',[10 5 10 6.5]);
% plot the standard spectrum
f = zeros(10,1);
A_true = zeros(10,1);
for i = 1:2
    f(1+5*(i-1)) = f_c(pp(i)) - 2 * f_l(qq(i));
    f(2+5*(i-1)) = f_c(pp(i)) - f_l(qq(i));
    f(3+5*(i-1)) = f_c(pp(i));
    f(4+5*(i-1)) = f_c(pp(i)) + f_l(qq(i));
    f(5+5*(i-1)) = f_c(pp(i)) + 2 * f_l(qq(i));
    m = 11/f_l(qq(i));
    A_true(1+5*(i-1)) = AA(i)*abs( 2*m*sin(m*pi/2)/((m^2-4)*pi) );
    A_true(2+5*(i-1)) = AA(i)*abs( 2*m*cos(m*pi/2)/((m^2-1)*pi) );
    A_true(3+5*(i-1)) = AA(i)*abs( 2*sin(m*pi/2)/(m*pi) );
    A_true(4+5*(i-1)) = A_true(2+5*(i-1));
    A_true(5+5*(i-1)) = A_true(1+5*(i-1));
end
hold on;
plot(f, A_true, '*k','MarkerSize',3.5);
% plot estimation results
hold on;
stem(fs_hat, 2*Amplis_hat, '-k','Marker','none','Linewidth',0.3);
if k_Emax(1) > 0
    hold on;
    stem([fs_hat(k_Emax(1)),fs_hat(k_Emax(2)),fs_hat(k_Emax(3))],...
        2*[Amplis_hat(k_Emax(1)),Amplis_hat(k_Emax(2)),Amplis_hat(k_Emax(3))]...
        ,'k','Linewidth',1.5,'Marker','none');
end
grid on;
axis([1670 1730 0 1]);
axes = gca;
axes.XTick = sort([f_c(pp(1)) f_c(pp(1))-f_l(qq(1)) f_c(pp(1))+f_l(qq(1))...
    f_c(pp(2)) f_c(pp(2))-f_l(qq(2)) f_c(pp(2))+f_l(qq(2))...
    f_c(pp(2))-2*f_l(qq(2)) f_c(pp(2))+2*f_l(qq(2)) ]);
axes.YTick = (0:0.2:1);
axes.Position = [0.13 0.23 0.82 0.76];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.XTickLabelRotation = 45;
axes.GridAlpha = 0.2;
box on;
xlabel({'\fontname{宋体}频率\fontname{Times new roman}/Hz'},'Fontsize',7.5);
ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
legend({'\fontname{宋体}实际频谱',...
    '\fontname{宋体}文献\fontname{Times new roman}[68]\fontname{宋体}算法所得频谱'...
    },'Fontsize',7.5,'Location','northoutside',...
    'Orientation','horizontal');
legend boxoff;
%% Fig 4.6 the spectrums estimated by RELAX alg with a constant K
% note: in RELAX_alpha_complex_model, replace line 52nd by K = Nfreqs;
% and delete the line 85th 
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f_c = [1701.4 1698.7 2001.4 1998.7 2301.4 2298.7 2601.4 2598.7];
f_l = 10.3:1.1:29;

pp(1) = 1; qq(1) = 11; pp(2) = 2; qq(2) = 1;
% pp(1) = 2; qq(1) = 18; pp(2) = 1; qq(2) = 3;

AA(1) = 1; AA(2) = 0.7;

Td = 0.3;

f0 = f_c(pp(1)); f1 = f_l(qq(1)); fs = 8000; 
s1 = AA(1)*Generate2000Signal(f0,f1,fs,Td);
snr = 30;
s1 = awgn(s1,snr,'measured');
f0 = f_c(pp(2)); f1 = f_l(qq(2)); fs = 8000; 
s2 = AA(2)*Generate2000Signal(f0,f1,fs,Td);
s = s1 + s2;

[ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = RELAX_constantK_complex_model(s,fs);

N = length(s);
M = 100 * N;

periodogram = abs(fft(s,M))/N;

reso_FFT = fs / M;
kl = round(1660/reso_FFT);
kr = round(1740/reso_FFT);

figure;
set(gcf,'unit','centimeters','position',[10 5 10 6.5]);
% plot the standard spectrum
f = zeros(10,1);
A_true = zeros(10,1);
for i = 1:2
    f(1+5*(i-1)) = f_c(pp(i)) - 2 * f_l(qq(i));
    f(2+5*(i-1)) = f_c(pp(i)) - f_l(qq(i));
    f(3+5*(i-1)) = f_c(pp(i));
    f(4+5*(i-1)) = f_c(pp(i)) + f_l(qq(i));
    f(5+5*(i-1)) = f_c(pp(i)) + 2 * f_l(qq(i));
    m = 11/f_l(qq(i));
    A_true(1+5*(i-1)) = AA(i)*abs( 2*m*sin(m*pi/2)/((m^2-4)*pi) );
    A_true(2+5*(i-1)) = AA(i)*abs( 2*m*cos(m*pi/2)/((m^2-1)*pi) );
    A_true(3+5*(i-1)) = AA(i)*abs( 2*sin(m*pi/2)/(m*pi) );
    A_true(4+5*(i-1)) = A_true(2+5*(i-1));
    A_true(5+5*(i-1)) = A_true(1+5*(i-1));
end
hold on;
plot(f, A_true, '*k','MarkerSize',3.5);
% plot estimation results
hold on;
stem(fs_hat, 2*Amplis_hat, '-k','Marker','none','Linewidth',0.3);
if k_Emax(1) > 0
    hold on;
    stem([fs_hat(k_Emax(1)),fs_hat(k_Emax(2)),fs_hat(k_Emax(3))],...
        2*[Amplis_hat(k_Emax(1)),Amplis_hat(k_Emax(2)),Amplis_hat(k_Emax(3))]...
        ,'k','Linewidth',1.5,'Marker','none');
end
grid on;
axis([1670 1730 0 1]);
axes = gca;
axes.XTick = sort([f_c(pp(1)) f_c(pp(1))-f_l(qq(1)) f_c(pp(1))+f_l(qq(1))...
    f_c(pp(2)) f_c(pp(2))-f_l(qq(2)) f_c(pp(2))+f_l(qq(2))...
    f_c(pp(2))-2*f_l(qq(2)) f_c(pp(2))+2*f_l(qq(2)) ]);
axes.YTick = (0:0.2:1);
axes.Position = [0.13 0.23 0.82 0.76];
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.XTickLabelRotation = 45;
axes.GridAlpha = 0.2;
box on;
xlabel({'\fontname{宋体}频率\fontname{Times new roman}/Hz'},'Fontsize',7.5);
ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
legend({'\fontname{宋体}实际频谱',...
    '\fontname{Times new roman}RELAX\fontname{宋体}法所得频谱'},...
    'Fontsize',7.5,'Location','northoutside',...
    'Orientation','horizontal');
legend boxoff;
%% Fig 4.7 ploting the periodogram in two-tone situation
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));
N = 2400;
n = (0:N-1)';
fs = 8000;
f1 = 1701.4; f2 = [1703 1703.5 1703.6 1705];
xlabels = cell(2,4);
xlabels(:,1) = {'\fontsize{7.5}\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 频差为\fontname{Times new roman}1.6Hz'};
xlabels(:,2) = {'\fontsize{7.5}\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 频差为\fontname{Times new roman}2.1Hz'};
xlabels(:,3) = {'\fontsize{7.5}\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(c)\fontname{宋体} 频差为\fontname{Times new roman}2.2Hz'};
xlabels(:,4) = {'\fontsize{7.5}\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(d)\fontname{宋体} 频差为\fontname{Times new roman}3.6Hz'};
figure;
set(gcf,'unit','centimeters','position',[5 5 10.5 8]);
for i = 1:4
    x = cos(2*pi*f1*n/fs) + cos(2*pi*f2(i)*n/fs);
    
    % obtaining periodogram
    fft_result = fft(x,20*N);
    fft_reso = fs/(20*N);
    periodogram = abs(fft_result)/N;
    % obtaining cost func of the ML estimator for single real tone
    M_half = 10*N; M_zpFFT = 20*N; n0 = -(N-1)/2;
    omega_zpFFT = 2*pi*(0:M_half-1)'/M_zpFFT;
    temp1 = sin(N*omega_zpFFT) ./ sin(omega_zpFFT);
    v1_FFT = 2 ./ (N + temp1);
    v2_FFT = 2 ./ (N - temp1);
    temp2 = exp(-1j*n0*omega_zpFFT).*fft_result(1:M_half);
    CTx_FFT = real(temp2); % M/2 dimmension vector
    STx_FFT = -imag(temp2); 
    J_FFT = v1_FFT.*CTx_FFT.^2 + v2_FFT.*STx_FFT.^2;

    kl = round(1695/fft_reso);
    kr = round(1712/fft_reso);

    subplot(2,2,i);
    plot(kl*fft_reso:fft_reso:kr*fft_reso, periodogram(kl+1:kr+1)/max(periodogram(kl+1:kr+1)),'-.k');
    hold on;
    plot(kl*fft_reso:fft_reso:kr*fft_reso, J_FFT(kl+1:kr+1)/max(J_FFT(kl+1:kr+1)),'-k');

    xlabel(xlabels(:,i),'Fontsize',7.5);
    ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
    axis([kl*fft_reso kr*fft_reso 0 1]);
    axes = gca;
    axes.XTick = round((0:fs/(2*N):fs-fs/(2*N))*10)/10;
    axes.XTickLabelRotation = 45;
    axes.YTick = 0:0.2:1;
    axes.LineWidth = 0.283;
    axes.FontSize = 7.5;
    axes.FontName = 'Times new roman';
    axes.GridAlpha = 0.2;
    grid on;
end
legend({'\fontname{宋体}周期图',...
    '\fontname{宋体}实单频信号极大似然估计价值函数'},...
    'Fontsize',7.5,'Orientation','horizontal','box','off');
% set(gcf,'WindowStyle','normal');

%% Fig 4.8 ploting the estimation results in two-tone situation
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));
N = 2400;
n = (0:N-1)';
fs = 8000;
f1 = 1701.4; f2 = [1702.1 1702.2 1702.3 1702.4];
xlabels = cell(2,4);
xlabels(:,1) = {'\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 频差为\fontname{Times new roman}0.7Hz'};
xlabels(:,2) = {'\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 频差为\fontname{Times new roman}0.8Hz'};
xlabels(:,3) = {'\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(c)\fontname{宋体} 频差为\fontname{Times new roman}0.9Hz'};
xlabels(:,4) = {'\fontname{宋体}频率\fontname{Times new roman}/Hz';...
    '\fontsize{10.5}\fontname{Times new roman}(d)\fontname{宋体} 频差为\fontname{Times new roman}1.0Hz'};
figure;
set(gcf,'unit','centimeters','position',[5 5 10.5 8]);
for i = 1:4
    x = cos(2*pi*f1*n/fs) + cos(2*pi*f2(i)*n/fs);

    [f0_hat, f1_hat, fs_hat, Amplis_hat, k_Emax ] = Rapid_RELAX_real_Dichotomous_Search(x,fs);
    
    fprintf('freqs are %f  %f  %f  %f  %f  %f\n',fs_hat(1),fs_hat(2),fs_hat(3),fs_hat(4),fs_hat(5),fs_hat(6));
    
    fft_result = fft(x,20*N);
    fft_reso = fs/(20*N);
    periodogram = abs(fft_result)/N;

    kl = round(1695/fft_reso);
    kr = round(1712/fft_reso);

    subplot(2,2,i);
    plot(kl*fft_reso:fft_reso:kr*fft_reso, periodogram(kl+1:kr+1),'k','Linewidth',0.3);
    hold on;
    stem(fs_hat,Amplis_hat,'k','Marker','none','Linewidth',1.5);
    hold on;
    plot([f1 f2(i)],[1 1],'*k');
    
    xlabel(xlabels(:,i),'Fontsize',7.5);
    ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
    tmp = max(periodogram(kl+1:kr+1));
    axis([kl*fft_reso kr*fft_reso 0 1]);
    axes = gca;
    axes.XTick = round((0:fs/(2*N):fs-fs/(2*N))*10)/10;
    axes.XTickLabelRotation = 45;
    axes.YTick = 0.1:0.2:1;
    axes.LineWidth = 0.283;
    axes.FontSize = 7.5;
    axes.FontName = 'Times new roman';
    axes.GridAlpha = 0.2;
    grid on;
end
legend({'\fontname{宋体}周期图',...
    '\fontname{宋体}快速\fontname{Times new roman}RELAX\fontname{宋体}法所得频谱',...
    '\fontname{宋体}实际频谱'},...
    'Fontsize',7.5,'Orientation','horizontal','box','off');
% set(gcf,'WindowStyle','normal');
