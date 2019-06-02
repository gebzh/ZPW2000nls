function [ ] = Load_Data_plus( N,fs )
% 1. calculating the parameters which are irrelevant to the received signal 
% 2. loading these parameters into the memory as global variables  
% hint: this function is written in accordance to the patent, thus the
% names of variables are similar to the patent.
%   inputs:
%   N: the length of the received signal; fs: the sampling frequency

% definition of the standard carrier freqs and the low freqs
f0_std = [1701.4 1698.7 2001.4 1998.7 2301.4 2298.7 2601.4  2598.7];
f1_std = 10.3:1.1:29;

global Z_AI;
global ZTZ_AI;
global ZTZ_AII;
global ks_AII;
global gridparas_AII;
global Amplis_std;
global lengthOfSignal;
global n0iszero;

lengthOfSignal = N;
n0iszero = 1;

n = (0:N-1)';
L = 5;
omega = zeros(L,1);

% (A.1)
Z_AI = zeros(N,10,8,18);
ZTZ_AI = zeros(10,10,8,18);
Amplis_std = zeros(5,18);
for p = 1:8
    f0 = f0_std(p);
    for q = 1:18
        f1 = f1_std(q);
        l = (1:L)';
        omega = 2*pi*(f0+(l-3)*f1)/fs;
        Z = [cos(omega(1)*n),cos(omega(2)*n),cos(omega(3)*n),...
             cos(omega(4)*n),cos(omega(5)*n),...
             sin(omega(1)*n),sin(omega(2)*n),sin(omega(3)*n),...
             sin(omega(4)*n),sin(omega(5)*n)];
        Z_AI(:,:,p,q) = Z;
        ZTZ_AI(:,:,p,q) = Z'*Z;
        if p==1
            m = 11/f1;
            Amplis_std(1,q) = abs(2*m*sin(m*pi/2)/((m^2-4)*pi));
            Amplis_std(2,q) = abs(2*m*cos(m*pi/2)/((m^2-1)*pi));
            Amplis_std(3,q) = abs(2*sin(m*pi/2)/(m*pi));
            Amplis_std(4,q) = Amplis_std(2,q);
            Amplis_std(5,q) = Amplis_std(1,q);
        end
    end
end

% (A.2)
Deltaomega_d = Grid_Size_Calc();
M = round(2*pi/Deltaomega_d);
Deltaomega = 2*pi/M;
L_h = round(0.6*pi/(fs*Deltaomega));
gridparas_AII = [M,L_h];

ZTZ_AII = zeros(10,10,2*L_h+1,2*L_h+1,8,18);
ks_AII = zeros(3,8,18);
for p = 1:8
    f0 = f0_std(p);
    for q = 1:18
        f1 = f1_std(q);
        k0 = round(2*pi*f0/(fs*Deltaomega)); 
        k1 = round(2*pi*f1/(fs*Deltaomega)); 
        k2 = round(4*pi*f1/(fs*Deltaomega));
        ks_AII(:,p,q) = [k0;k1;k2];
        for i = -L_h : L_h
            omega0 = (k0+i)*Deltaomega;
            for j = -L_h : L_h 
                omega(1) = omega0 - (k2+2*j)*Deltaomega;
                omega(2) = omega0 - (k1+j)*Deltaomega;
                omega(3) = omega0;
                omega(4) = omega0 + (k1+j)*Deltaomega;
                omega(5) = omega0 + (k2+2*j)*Deltaomega;
                Z = [cos(omega(1)*n),cos(omega(2)*n),cos(omega(3)*n),...
                     cos(omega(4)*n),cos(omega(5)*n),...
                     sin(omega(1)*n),sin(omega(2)*n),sin(omega(3)*n),...
                     sin(omega(4)*n),sin(omega(5)*n)];
                ZTZ_AII(:,:,i+L_h+1,j+L_h+1,p,q) = Z'*Z;
            end
        end
    end
end

    function delta_omega = Grid_Size_Calc()
        % calculating the grid size of the grid search algorithm (see sec4)
        % the model order of the ZPW-2000 signal is set as 3.
        gr = 1.15;  
        delta_omega = sqrt(12*(gr-1)/gr)/N/L;
    end


end

