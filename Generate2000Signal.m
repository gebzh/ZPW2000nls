function out = Generate2000Signal(f0,f1,fs,Td,phi0)
% 生成振幅为1的ZPW-2000仿真信号
%   输入参数:
%   f0:载频频率; f1:低频频率; fs:采样频率; Td:采样持续时间; phi0:初相位
%   输出参数:
%   out:ZPW-2000仿真信号采样序列
    f_delta = 11;                   %频偏
    N = round(Td*fs);               %信号长度
    t = ( 0:1/fs:(N-1)/fs )';       %采样时刻定义
    if nargin < 5
        phi0 = 2*pi*rand(1);            %初相位
    end
    
    % 产生频率为f1的方波信号,幅值为+-1
    g = square(2*pi*f1*t);                      
    % 对方波积分，得到相对于时间，斜率为+-1的连续的三角波
    phi = cumsum(g)/fs;                       
    % 三角波信号乘上频偏的角频率，与载波信号相位及初相位相加，即可得信号相位
    % 将信号相位代入cos函数，即可得时域振幅为1的ZPW-2000仿真信号采样序列
    out = cos(2*pi*f0*t + 2*pi*f_delta*phi + phi0);      
end
