%% Fig 2.2 the generation of the ZPW-2000 signal
f0 = 1701.4; f1 = 29; fs = 8000; Td = 0.07; phi0 = 0;
f_delta = 11; N = round(Td*fs); t = ( 0:1/fs:(N-1)/fs )';

g = square(2*pi*f1*t);
phi = cumsum(g)/fs;
out = cos(2*pi*f0*t + 2*pi*f_delta*phi + phi0); 

figure;
set(gcf,'unit','centimeters','position',[5 5 10.5 9]);
subplot(3,1,1),plot(0:1/fs:Td-1/fs,g,'k'); axis([0 Td -1.5 1.5]);
xlabel({'\fontname{宋体}时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(a)'},'Fontsize',7.5); 
ylabel('$g(t)$','interpreter','latex','Fontsize',7.5);
axes = gca;
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;

subplot(3,1,2),plot(0:1/fs:Td-1/fs,phi,'k'); 
axis([0 Td min(phi)-0.005 max(phi)+0.005]);
xlabel({'\fontname{宋体}时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(b)'},'Fontsize',7.5); 
ylabel('$\phi(t)$','interpreter','latex','Fontsize',7.5);
axes = gca;
axes.LineWidth = 0.283;
axes.YTick = [0 0.01 0.02];
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;

subplot(3,1,3),plot(0:1/fs:Td-1/fs,out,'k'); axis([0 Td -1.5 1.5]);
xlabel({'\fontname{宋体}时间\fontname{Times new roman}/s';...
    '\fontsize{10.5}\fontname{Times new roman}(c)'},'Fontsize',7.5); 
ylabel('$s(t)$','interpreter','latex','Fontsize',7.5);
axes = gca;
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;
%% Fig 2.3 the ideal periodogram of the ZPW-2000 signal
f0 = 1701.4; f1 = 29; fs = 8000; Td = 20; phi0 = 0;

% generating signal
s = Generate2000Signal(f0,f1,fs,Td,phi0);
N = length(s);

% acquiring the periodogram
periodogram = 2*abs(fft(s))/N;
reso_FFT = fs/N;
kl = round(1635/reso_FFT); kr = round(1765/reso_FFT);

% obtaining the standard amplitudes
m = 11 / f1; A = 1;
A_0 = abs(2*A*sin(m*pi/2)/(m*pi));
A_pm1 = abs(2*m*A*cos(m*pi/2)/((m^2-1)*pi));
A_pm2 = abs(2*m*A*sin(m*pi/2)/((m^2-4)*pi));
A_pm3 = abs(2*m*A*cos(m*pi/2)/((m^2-9)*pi));

figure;
set(gcf,'unit','centimeters','position',[5 5 8 5]);
stem(kl*reso_FFT:reso_FFT:kr*reso_FFT,periodogram(kl+1:kr+1),'k',...
    'Marker','none');
axis([kl*reso_FFT kr*reso_FFT 0 1]);
xlabel('\fontname{宋体}频率\fontname{Times new roman}/Hz','Fontsize',7.5); 
ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
axes = gca;
axes.LineWidth = 0.283;
axes.XTick = [f0-3*f1 f0-2*f1 f0-f1 f0 f0+f1 f0+2*f1 f0+3*f1];
axes.YTick = [A_pm2 A_pm1 A_0];
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;

%% Fig 2.4 the impulse interference
Ts = 1/8000; 
Td = 0.3;
N = round(Td/Ts);
t_impulse = (0:Ts:0.04)';
t_traction = (0:Ts:(N-length(t_impulse)-401)*Ts)';

A_im = 1111; a = 40; omega = 2*pi*300; 
y_impulse = A_im*exp(-a*t_impulse).*sin(omega*t_impulse);

y_traction = 222*sin(2*pi*50*t_traction);

y = [zeros(400,1);y_impulse;y_traction];
t = (0:Ts:Ts*(length(y)-1));

figure;
set(gcf,'unit','centimeters','position',[5 5 8 5]);
plot(t,y,'k');
xlabel('\fontname{宋体}时间\fontname{Times new roman}/s','Fontsize',7.5); 
ylabel('\fontname{宋体}幅值\fontname{Times new roman}/V','Fontsize',7.5);
axis([0 Td -1200 1200]);
axes = gca;
axes.XTick = (0:0.05:0.3);
axes.YTick = (-1200:200:1200);
axes.LineWidth = 0.283;
axes.FontSize = 7.5;
axes.FontName = 'Times new roman';
axes.GridAlpha = 0.2;
grid on;