%% tests for RELAX alg
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

f0 = 2301.4; f1=10.3; fs = 8000; Td = 1;

m = 11/f1;
A_true(1) = abs( 2*m*sin(m*pi/2)/((m^2-16)*pi) );
A_true(2) = abs( 2*m*cos(m*pi/2)/((m^2-9)*pi) );
A_true(3) = abs( 2*m*sin(m*pi/2)/((m^2-4)*pi) );
A_true(4) = abs( 2*m*cos(m*pi/2)/((m^2-1)*pi) );
A_true(5) = abs( 2*sin(m*pi/2)/(m*pi) );
A_true(6) = A_true(4);
A_true(7) = A_true(3);
A_true(8) = A_true(2);
A_true(9) = A_true(1);
f_true = [f0-4*f1, f0-3*f1, f0-2*f1, f0-f1, f0, ...
    f0+f1, f0+2*f1, f0+3*f1, f0+4*f1]';

s = Generate2000Signal(f0,f1,fs,Td,0);
N = length(s);
n = (0:N-1)';
snr = 5;
s = awgn(s,snr,'measured') + 0.8*cos(2*pi*(f0-0.8)/fs*n);
[ fs_hat, Amplis_hat, f0_hat, f1_hat, k_compo ] = RELAX_alpha_complex_model(s,fs);
% [ fs_hat, Alphas_hat ] = RELAX_alpha(s,fs);
% Amplis_hat = sqrt(Alphas_hat(1,:).^2 + Alphas_hat(2,:).^2);

Periodo1 = abs(fft(s,100*length(s)))/N;

reso_FFT = fs/N/100;
kl = round(2250/reso_FFT)+1;
kr = round(2350/reso_FFT)+1;

figure;
set(gcf,'unit','centimeters','position',[10 5 13 8]);
stem(fs_hat,2*Amplis_hat,'k','Linewidth',0.9,'Marker','none');
hold on;
plot(f_true,A_true,'*k');
hold on;
plot((kl-1)*reso_FFT:reso_FFT:(kr-1)*reso_FFT,2*Periodo1(kl:kr),...
    '--k','Linewidth',0.45);
hold on;
stem([fs_hat(k_compo(1)),fs_hat(k_compo(2)),fs_hat(k_compo(3)),...
    fs_hat(k_compo(4)),fs_hat(k_compo(5))],...
    2*[Amplis_hat(k_compo(1)),Amplis_hat(k_compo(2)),Amplis_hat(k_compo(3))...
    ,Amplis_hat(k_compo(4)),Amplis_hat(k_compo(5))]...
    ,'r','Linewidth',0.9,'Marker','none');
axis([2250 2350 0 0.9]);
axes = gca;
axes.XTick = [2250 2275 2300 2325 2350];
axes.YTick = (0.1:0.1:0.9);
axes.LineWidth = 0.283;
grid on;
xlabel('频率/Hz'); ylabel('幅值/V');
legend({'RELAX算法频谱','ZPW-2000理想频谱','周期图'});

%% testing the immunity of interferences from the neighboring track
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));
addpath(strcat(pwd,'/pkgFromWeb/simonhenin-columnlegend-8883602/'))

f_c = [1701.4 1698.7 2001.4 1998.7 2301.4 2298.7 2601.4 2598.7];
f_l = 10.3:1.1:29;

% pp(1) = 2; qq(1) = 18; pp(2) = 1; qq(2) = 6;
pp(1) = 1; qq(1) = 11; pp(2) = 2; qq(2) = 1;
AA(1) = 1; AA(2) = 0.7;

Td = 0.3;

f0 = f_c(pp(1)); f1 = f_l(qq(1)); fs = 8000; 
s1 = AA(1)*Generate2000Signal(f0,f1,fs,Td,0);
snr = 30;
s1 = awgn(s1,snr,'measured');
f0 = f_c(pp(2)); f1 = f_l(qq(2)); fs = 8000; 
s2 = AA(2)*Generate2000Signal(f0,f1,fs,Td,0);
s = s1 + s2;


[ fs_hat, Amplis_hat, f0_hat, f1_hat, k_compo ] = RELAX_alpha_complex_model(s,fs);

N = length(s);
M = 100 * N;

periodogram = abs(fft(s,M))/N;
% periodogram = abs(fft(s,M))/N;

reso_FFT = fs / M;
kl = round(1660/reso_FFT);
kr = round(1740/reso_FFT);

figure;
set(gcf,'unit','centimeters','position',[10 5 5.8 5.4]);
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
stem(fs_hat, 2*Amplis_hat, '-k','Marker','none');
% plot the periodogram
hold on;
plot(kl*reso_FFT:reso_FFT:kr*reso_FFT,2*periodogram(kl+1:kr+1),'-.k','Linewidth',0.5);
hold on;
stem([fs_hat(k_compo(1)),fs_hat(k_compo(2)),fs_hat(k_compo(3)),...
    fs_hat(k_compo(4)),fs_hat(k_compo(5))],...
    2*[Amplis_hat(k_compo(1)),Amplis_hat(k_compo(2)),Amplis_hat(k_compo(3))...
    ,Amplis_hat(k_compo(4)),Amplis_hat(k_compo(5))]...
    ,'r','Linewidth',0.9,'Marker','none');
box off;
grid on;
axis([1660 1740 0 1]);
axes = gca;
axes.XTick = sort([1698.7-29 1698.7 1701.4+15.8 1698.7+29]);
axes.YTick = (0:0.2:1);
axes.Position = [0.15 0.23 0.82 0.61];
axes.LineWidth = 0.283;
axes.FontSize = 7;
axes.FontName = 'Times new roman';
axes.XTickLabelRotation = 45;
xlabel({'\fontname{宋体}频率\fontname{Times new roman}(Hz)'},'Fontsize',7);
ylabel('\fontname{宋体}幅值','Fontsize',7);
legend({'\fontname{宋体}实际频谱',...
    '\fontname{宋体}粗估计所得频谱',...
    '\fontname{宋体}汉宁窗截取信号周期图'},'Fontsize',7);
legend boxoff;
[f0_h,f1_h] = YangFan2010(s,8000);
fprintf('now the estimated f0 and f1 of yf method are %f %f respectively\n',...
    f0_h, f1_h);
%% ploting the periodogram in two-tone situation
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));
N = 2400;
n = (0:N-1)';
fs = 8000;
f1 = 1701.4; f2 = 1703.6;
x = cos(2*pi*f1*n/fs) + cos(2*pi*f2*n/fs);

% [ fs_hat, Amplis_hat, f0_hat, f1_hat, k_Emax ] = Rapid_RELAX_alpha_Dichotomous_Search(x,fs);

fft_result = fft(x,20*N);
fft_reso = fs/(20*N);
periodogram = fft_result.*conj(fft_result);

kl = round(1695/fft_reso);
kr = round(1712/fft_reso);

figure;
plot(kl*fft_reso:fft_reso:kr*fft_reso, periodogram(kl+1:kr+1));
hold on;
stem([1703.33333333333          1701.66666666667          1704.16666666667],ones(3,1)*max(periodogram(kl+1:kr+1)),'k','Marker','none');
axis([kl*fft_reso kr*fft_reso 0 max(periodogram(kl+1:kr+1))+10]);
axes = gca;
axes.XTick = round((0:fs/N:fs-fs/N)*10)/10;
grid on;
%% ploting the time complexities of different algs
frps = 0.0001:0.0005:0.1;
Ns = 1000:500:10000;
fs = 8000;

TCs = zeros(length(frps),length(Ns),4);
% ref 19
for i = 1:length(frps)
    frp = frps(i);
    for j = 1:length(Ns)
        N = Ns(j);
        M19 = round(fs/(26*frp)); 
        TCs(i,j,1) = 25*N + M19*log2(M19);
    end
end
% ref20
for i = 1:length(frps)
    frp = frps(i);
    for j = 1:length(Ns)
        N = Ns(j);
        Q20 = round(log2(fs/(N*frp))); 
        TCs(i,j,2) = N*log2(N) + Q20*N;
    end
end
% NLSM
for i = 1:length(frps)
    frp = frps(i);
    for j = 1:length(Ns)
        N = Ns(j);
        N_c = 144; L = 10; M_NLSM = 25*N; 
        Q_NLSM = round(log2(fs/(M_NLSM*frp)));
        TCs(i,j,3) = N_c*L*N + M_NLSM*log2(M_NLSM) + (Q_NLSM^2)*(L^2)*N;
    end
end
% rRELAX
for i = 1:length(frps)
    frp = frps(i);
    for j = 1:length(Ns)
        N = Ns(j);
        M_rRELAX = 4*N; 
        Q_k = 4; L_bar = 12;
        TCs(i,j,4) = M_rRELAX*log2(M_rRELAX)*Q_k*(1+L_bar)*L_bar/2;
    end
end
figure;
subplot(2,2,1),mesh(Ns,frps,TCs(:,:,1));

subplot(2,2,2),mesh(Ns,frps,TCs(:,:,2));
subplot(2,2,3),mesh(Ns,frps,TCs(:,:,3));
subplot(2,2,4),mesh(Ns,frps,TCs(:,:,4));
