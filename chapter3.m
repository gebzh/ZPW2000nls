%% Fig 3.1 comparison of the ideal model and the approximate model
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f0 = 2001.4; f1=29; fs = 8000; Td = 0.3;
x = Generate2000Signal(f0,f1,fs,Td,0);

[f0_h11,f1_h11,alpha] = NLSM_Based_Algorithm_plus( x,fs );
N = length(x); n0 = 0; n = (n0:n0+N-1)';
omg0 = 2*pi*f0_h11/fs; omg1 = 2*pi*f1_h11/fs;
Zl = [cos((omg0-2*omg1)*n),cos((omg0-omg1)*n),cos(omg0*n),...
    cos((omg0+omg1)*n),cos((omg0+2*omg1)*n),...
    sin((omg0-2*omg1)*n),sin((omg0-omg1)*n),sin(omg0*n),...
    sin((omg0+omg1)*n),sin((omg0+2*omg1)*n)];
x_tilde = Zl*alpha;

figure;
set(gcf,'unit','centimeters','position',[10 5 7.5 8]);
subplot(2,1,1);
TSegment = 0.015:1/fs:0.02;
nStart = round(0.015*fs);
plot(TSegment, x( nStart:nStart+length(TSegment)-1 ),'-r','LineWidth',0.567);
axis([0.015,0.02,-1,1]);
axes = gca;
axes.LineWidth = 0.283;
axes.Position = [0.16 0.64 0.77 0.28];
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
hold on; 
plot(TSegment, x_tilde(nStart:nStart+length(TSegment)-1),'--k','LineWidth',0.567);
grid on;
legend({'$\boldmath{x_1}$','$\boldmath{\tilde{x}_1}$'},'Interpreter','latex',...
    'Position',[0.05 0.95 0.94 0.03],'Orientation','horizontal','Box','off',...
    'Fontsize',7.5);
xlabel({'\fontname{宋体}时间/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)\fontname{宋体} 时域波形对比'},'Fontsize',7.5); 
ylabel('\fontname{宋体}信号\fontname{Times new roman}/V','Fontsize',7.5);

subplot(2,1,2);
plot(0:1/fs:Td-1/fs,abs(x-x_tilde),'k','LineWidth',0.567);
grid on; 
axis([0,Td,0,1]);
axes = gca;
axes.LineWidth = 0.283;
axes.Position = [0.16 0.164 0.77 0.28];
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
xlabel({'\fontname{宋体}时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)\fontname{宋体} 绝对误差随时间变化情况'},'Fontsize',7.5); 
ylabel('\fontname{宋体}绝对误差\fontname{Time s new roman}/V','Fontsize',7.5);
%% Fig 3.2 periodogram of a sinosoid of frequency equaling 1701.4 Hz
fs = 8000;
t = 0:1/fs:2399/fs;
x = awgn(exp(1j*2*pi*1701.4*t),30,'measured');
N = length(x);

temp = fft(x,length(x)*100);
periodogram = temp.*conj(temp)/N;
P0 = max(periodogram);
periodogram_dB = 10*log10(periodogram/P0);

N = length(x);
reso_FFT = fs/(100*N);
kl = round((1701.4-30)/reso_FFT);
kr = round((1701.4+30)/reso_FFT);

figure;
set(gcf,'unit','centimeters','position',[5 5 8 6]);
plot(reso_FFT*(kl:kr),periodogram_dB(kl+1:kr+1),'k','LineWidth',0.3);
axis([kl*reso_FFT (kr+1)*reso_FFT -60 0]);
axes = gca;
axes.Position = [0.15 0.22 0.8 0.75];
axes.XTick = round((1701.4 + fs*(-9:2:9)/N)*10)/10;
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
axes.XTickLabelRotation = 45;
xlabel('$f$/Hz','Interpreter','latex','Fontsize',7.5); 
ylabel('$10log_{10}[W_{\mathrm{f0}}(f)/max\{W_{\mathrm{f0}}(f)\}]/\mathrm{dB}$',...
    'Interpreter','latex','Fontsize',7.5);
grid on;
%% Fig 3.3 periodogram of a ZPW-2000 signal
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f0 = 1701.4; f1 = 10.3; fs = 8000; Td = 0.3; 
s = Generate2000Signal(f0,f1,fs,Td);
x = awgn(s,30,'measured');
N = length(x);

% obtaining periodogram
temp = fft(x,N*100);
periodogram = temp.*conj(temp)/N;
P0 = max(periodogram);
periodogram_dB = 10*log10(periodogram/P0);

% setting the freq range of the plot
reso_FFT = fs/(100*N);
kl = round((1701.4-25)/reso_FFT);
kr = round((1701.4+25)/reso_FFT);

% obtaining the local maximum points of the sinusoidal components 
mainlobe_width = 2*fs/N;
freq_components = 1701.4 + f1*(-2:2)';
kkl = round((freq_components-mainlobe_width/4)/reso_FFT); 
kkr = round((freq_components+mainlobe_width/4)/reso_FFT); 
pos_maximum = zeros(5,1); 
W_maximum = zeros(5,1); 
for i = 1:5
    [W_maximum(i),pos_maximum(i)] = max(periodogram_dB(kkl(i)+1:kkr(i)+1));    
end
pos_maximum = pos_maximum + kkl;

figure;
set(gcf,'unit','centimeters','position',[5 5 8 6]);
plot(reso_FFT*(kl:kr),periodogram_dB(kl+1:kr+1),'k','LineWidth',0.3);
hold on;
plot((pos_maximum-1)*reso_FFT,W_maximum,'.k','MarkerSize',6);
axis([kl*reso_FFT (kr+1)*reso_FFT -40 0]);
axes = gca;
axes.Position = [0.15 0.22 0.8 0.75];
axes.XTick = round((1701.4 + f1*(-2:2))*10)/10;
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
axes.XTickLabelRotation = 45;
xlabel('$f$/Hz','Interpreter','latex','Fontsize',7.5); 
ylabel('$10log_{10}[W_{\mathrm{x2}}(f)/max\{W_{\mathrm{x2}}(f)\}]/\mathrm{dB}$',...
    'Interpreter','latex','Fontsize',7.5);
grid on;
%% Fig 3.4 periodogram of the 4-compoents-reduced ZPW-2000 signal
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f0 = 1701.4; f1 = 10.3; fs = 8000; Td = 0.3; phi0 = 0;
s = Generate2000Signal(f0,f1,fs,Td);
x = awgn(s,30,'measured');
N = length(x);

[f0_hat,f1_hat,alpha_hat] = NLSM_Based_Algorithm_plus(x,fs);

% reducing the other 4 freq components
n = (0:N-1)';
omg = 2*pi*(f0_hat-f1_hat*2:f1_hat:f0_hat+f1_hat*2)'/fs;
Z_four_components = [cos(omg(2)*n),cos(omg(3)*n),cos(omg(4)*n),cos(omg(5)*n),...
                    sin(omg(2)*n),sin(omg(3)*n),sin(omg(4)*n),sin(omg(5)*n)];
x = x - Z_four_components*[alpha_hat(2:5);alpha_hat(7:10)];

% obtaining periodogram
temp = fft(x,N*100);
periodogram = temp.*conj(temp)/N;
P0 = max(periodogram);
periodogram_dB = 10*log10(periodogram/P0);

% setting the freq range of the plot
reso_FFT = fs/(100*N);
kl = round((1701.4-25)/reso_FFT);
kr = round((1701.4+25)/reso_FFT);

% obtaining the local maximum points of the sinusoidal components 
mainlobe_width = 2*fs/N;
freq_components = 1701.4 + f1*(-2:2)';
kkl = round((freq_components-mainlobe_width/4)/reso_FFT); 
kkr = round((freq_components+mainlobe_width/4)/reso_FFT); 
pos_maximum = zeros(5,1); 
W_maximum = zeros(5,1); 
for i = 1:5
    [W_maximum(i),pos_maximum(i)] = max(periodogram_dB(kkl(i)+1:kkr(i)+1));    
end
pos_maximum = pos_maximum + kkl;

figure;
set(gcf,'unit','centimeters','position',[5 5 8 6]);
plot(reso_FFT*(kl:kr),periodogram_dB(kl+1:kr+1),'k','LineWidth',0.3);
hold on;
plot((pos_maximum(1)-1)*reso_FFT,W_maximum(1),'.k','MarkerSize',6);
axis([kl*reso_FFT (kr+1)*reso_FFT -40 0]);
axes = gca;
axes.Position = [0.15 0.22 0.8 0.75];
axes.XTick = round((1701.4 + f1*(-2:2))*10)/10;
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
axes.XTickLabelRotation = 45;
xlabel('$f$/Hz','Interpreter','latex','Fontsize',7.5); 
ylabel('$10log_{10}[W_{\mathrm{x3}}(f)/max\{W_{\mathrm{x3}}(f)\}]/\mathrm{dB}$',...
    'Interpreter','latex','Fontsize',7.5);
grid on;