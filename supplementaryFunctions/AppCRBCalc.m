function [ CRBf0, CRBf1 ] = AppCRBCalc( noise_free_signal,fs,f0,f1,sigmaSquare )
%Calculating exact CRB for the approximate model
%   inputs:
%   noise_free_signal: noise free ZPW-2000 signal
%   fs: the sampling frequency
%   f0: the carrier freq; f1: the low freq
%   sigmaSquare: the variance of the white Gaussian noise
%   outputs:
%   CRBf0: the CRB of the carrier freq
%   CRBf1: the CRB of the low freq

N = length(noise_free_signal);
n0 = -(N-1)/2;
n = (n0:n0+N-1)';
[ ~, ~, alpha_hat] = NLSM_Based_Algorithm_plus( noise_free_signal,fs );
omega0 = 2*pi*f0/fs;
omega1 = 2*pi*f1/fs;

% ==========* loading data *==========
global f0_std;
global f1_std;
global mies;
if isempty(f0_std) || isempty(f1_std) || isempty(mies)
    f0_std = [1701.4; 1698.7; 2001.4; 1998.7; 2301.4; 2298.7; 2601.4; 2598.7];
    f1_std = (10.3:1.1:29)';
    mies = 11./f1_std;
end

% ==========* calculating the Fisher information matrix *==========

% calculating the partial derivative with respect to omega_0 i.e. xi0
pd_c_1 = -sin((omega0-2*omega1)*n).*n;
pd_c_2 = -sin((omega0-omega1)*n).*n;
pd_c_3 = -sin(omega0*n).*n;
pd_c_4 = -sin((omega0+omega1)*n).*n;
pd_c_5 = -sin((omega0+2*omega1)*n).*n;
pd_s_1 = cos((omega0-2*omega1)*n).*n;
pd_s_2 = cos((omega0-omega1)*n).*n;
pd_s_3 = cos(omega0*n).*n;
pd_s_4 = cos((omega0+omega1)*n).*n;
pd_s_5 = cos((omega0+2*omega1)*n).*n;
xi0 = [pd_c_1,pd_c_2,pd_c_3,pd_c_4,pd_c_5,pd_s_1,pd_s_2,pd_s_3,pd_s_4,pd_s_5]...
    *alpha_hat;
% calculating the partial derivative with respect to omega_2 i.e. xi1
pd_c_1 = -sin((omega0-2*omega1)*n).*(-2*n);
pd_c_2 = -sin((omega0-omega1)*n).*(-n);
pd_c_3 = zeros(N,1);
pd_c_4 = -sin((omega0+omega1)*n).*n;
pd_c_5 = -sin((omega0+2*omega1)*n).*(2*n);
pd_s_1 = cos((omega0-2*omega1)*n).*(-2*n);
pd_s_2 = cos((omega0-omega1)*n).*(-n);
pd_s_3 = zeros(N,1);
pd_s_4 = cos((omega0+omega1)*n).*n;
pd_s_5 = cos((omega0+2*omega1)*n).*(2*n);
xi1 = [pd_c_1,pd_c_2,pd_c_3,pd_c_4,pd_c_5,pd_s_1,pd_s_2,pd_s_3,pd_s_4,pd_s_5]...
    *alpha_hat;
% calculating the partial derivative with respect to a_l i.e. psi
Z = [cos((omega0-2*omega1)*n) cos((omega0-omega1)*n) cos(omega0*n)...
    cos((omega0+omega1)*n) cos((omega0+2*omega1)*n)...
    sin((omega0-2*omega1)*n) sin((omega0-omega1)*n) sin(omega0*n)...
    sin((omega0+omega1)*n) sin((omega0+2*omega1)*n)];
psi = zeros(N,5);
psi(:,1) = Z*[1;0;0;0;0; 0;0;0;0;0];
psi(:,2) = Z*[0;1;0;0;0; 0;0;0;0;0];
psi(:,3) = Z*[0;0;1;0;0; 0;0;0;0;0];
psi(:,4) = Z*[0;0;0;1;0; 0;0;0;0;0];
psi(:,5) = Z*[0;0;0;0;1; 0;0;0;0;0];
% calculating the partial derivative with respect to b_l i.e. zeta
zeta = zeros(N,5);
zeta(:,1) = Z*[0;0;0;0;0; -1;0;0;0;0];
zeta(:,2) = Z*[0;0;0;0;0; 0;-1;0;0;0];
zeta(:,3) = Z*[0;0;0;0;0; 0;0;-1;0;0];
zeta(:,4) = Z*[0;0;0;0;0; 0;0;0;-1;0];
zeta(:,5) = Z*[0;0;0;0;0; 0;0;0;0;-1];
% obtaining the Fisher info matrix 
I_theta = [xi0'*xi0, xi0'*xi1, xi0'*psi(:,1), xi0'*psi(:,2), xi0'*psi(:,3),...
    xi0'*psi(:,4), xi0'*psi(:,5), xi0'*zeta(:,1), xi0'*zeta(:,2), ...
    xi0'*zeta(:,3), xi0'*zeta(:,4), xi0'*zeta(:,5);
    
    xi1'*xi0, xi1'*xi1, xi1'*psi(:,1), xi1'*psi(:,2), xi1'*psi(:,3),...
    xi1'*psi(:,4), xi1'*psi(:,5), xi1'*zeta(:,1), xi1'*zeta(:,2), ...
    xi1'*zeta(:,3), xi1'*zeta(:,4), xi1'*zeta(:,5);
    
    psi(:,1)'*xi0, psi(:,1)'*xi1, psi(:,1)'*psi(:,1), ...
    psi(:,1)'*psi(:,2), psi(:,1)'*psi(:,3), psi(:,1)'*psi(:,4), ...
    psi(:,1)'*psi(:,5), psi(:,1)'*zeta(:,1), psi(:,1)'*zeta(:,2), ...
    psi(:,1)'*zeta(:,3), psi(:,1)'*zeta(:,4), psi(:,1)'*zeta(:,5);
    
    psi(:,2)'*xi0, psi(:,2)'*xi1, psi(:,2)'*psi(:,1), ...
    psi(:,2)'*psi(:,2), psi(:,2)'*psi(:,3), psi(:,2)'*psi(:,4), ...
    psi(:,2)'*psi(:,5), psi(:,2)'*zeta(:,1), psi(:,2)'*zeta(:,2), ...
    psi(:,2)'*zeta(:,3), psi(:,2)'*zeta(:,4), psi(:,2)'*zeta(:,5);
    
    psi(:,3)'*xi0, psi(:,3)'*xi1, psi(:,3)'*psi(:,1), ...
    psi(:,3)'*psi(:,2), psi(:,3)'*psi(:,3), psi(:,3)'*psi(:,4), ...
    psi(:,3)'*psi(:,5), psi(:,3)'*zeta(:,1), psi(:,3)'*zeta(:,2), ...
    psi(:,3)'*zeta(:,3), psi(:,3)'*zeta(:,4), psi(:,3)'*zeta(:,5);
    
    psi(:,4)'*xi0, psi(:,4)'*xi1, psi(:,4)'*psi(:,1), ...
    psi(:,4)'*psi(:,2), psi(:,4)'*psi(:,3), psi(:,4)'*psi(:,4), ...
    psi(:,4)'*psi(:,5), psi(:,4)'*zeta(:,1), psi(:,4)'*zeta(:,2), ...
    psi(:,4)'*zeta(:,3), psi(:,4)'*zeta(:,4), psi(:,4)'*zeta(:,5);
    
    psi(:,5)'*xi0, psi(:,5)'*xi1, psi(:,5)'*psi(:,1), ...
    psi(:,5)'*psi(:,2), psi(:,5)'*psi(:,3), psi(:,5)'*psi(:,4), ...
    psi(:,5)'*psi(:,5), psi(:,5)'*zeta(:,1), psi(:,5)'*zeta(:,2), ...
    psi(:,5)'*zeta(:,3), psi(:,5)'*zeta(:,4), psi(:,5)'*zeta(:,5);    
    
    zeta(:,1)'*xi0, zeta(:,1)'*xi1, zeta(:,1)'*psi(:,1), zeta(:,1)'*psi(:,2),...
    zeta(:,1)'*psi(:,3), zeta(:,1)'*psi(:,4), zeta(:,1)'*psi(:,5), ...
    zeta(:,1)'*zeta(:,1), zeta(:,1)'*zeta(:,2), zeta(:,1)'*zeta(:,3), ...
    zeta(:,1)'*zeta(:,4), zeta(:,1)'*zeta(:,5);    
    
    zeta(:,2)'*xi0, zeta(:,2)'*xi1, zeta(:,2)'*psi(:,1), zeta(:,2)'*psi(:,2),...
    zeta(:,2)'*psi(:,3), zeta(:,2)'*psi(:,4), zeta(:,2)'*psi(:,5), ...
    zeta(:,2)'*zeta(:,1), zeta(:,2)'*zeta(:,2), zeta(:,2)'*zeta(:,3), ...
    zeta(:,2)'*zeta(:,4), zeta(:,2)'*zeta(:,5);    
    
    zeta(:,3)'*xi0, zeta(:,3)'*xi1, zeta(:,3)'*psi(:,1), zeta(:,3)'*psi(:,2),...
    zeta(:,3)'*psi(:,3), zeta(:,3)'*psi(:,4), zeta(:,3)'*psi(:,5), ...
    zeta(:,3)'*zeta(:,1), zeta(:,3)'*zeta(:,2), zeta(:,3)'*zeta(:,3), ...
    zeta(:,3)'*zeta(:,4), zeta(:,3)'*zeta(:,5);    
    
    zeta(:,4)'*xi0, zeta(:,4)'*xi1, zeta(:,4)'*psi(:,1), zeta(:,4)'*psi(:,2),...
    zeta(:,4)'*psi(:,3), zeta(:,4)'*psi(:,4), zeta(:,4)'*psi(:,5), ...
    zeta(:,4)'*zeta(:,1), zeta(:,4)'*zeta(:,2), zeta(:,4)'*zeta(:,3), ...
    zeta(:,4)'*zeta(:,4), zeta(:,4)'*zeta(:,5);    
    
    zeta(:,5)'*xi0, zeta(:,5)'*xi1, zeta(:,5)'*psi(:,1), zeta(:,5)'*psi(:,2),...
    zeta(:,5)'*psi(:,3), zeta(:,5)'*psi(:,4), zeta(:,5)'*psi(:,5), ...
    zeta(:,5)'*zeta(:,1), zeta(:,5)'*zeta(:,2), zeta(:,5)'*zeta(:,3), ...
    zeta(:,5)'*zeta(:,4), zeta(:,5)'*zeta(:,5);];
        
I_theta = I_theta / sigmaSquare;

% ==========* calculating CRB of f0 & f1 *==========

% obtaining the inverse of the Fisher info matrix 
I = eye(size(I_theta,1));
I_theta_inv = zeros(size(I_theta));
for i = 1:size(I_theta,1)
    I_theta_inv(:,i) = GaussianElimination(I_theta,I(:,i));
end

% obtaining the CRB of f0 & f1
CRBf0 = (fs^2)*I_theta_inv(1,1)/(4*pi^2);
CRBf1 = (fs^2)*I_theta_inv(2,2)/(4*pi^2);

%% ========================* functions *==========================
    function x = GaussianElimination(A,b)
        % Solving the linear system Ax = b via Gaussian elimination
        % Designed based on Gilbert Strang. Linear Algebra and Its
        % Application fourth edition Chapter 1
        % Matrix A must be non-singular that the system can be solved,
        % thus if A is singular x in return will be NaN.
        % Cause Z_l'*Z_l is certainly non-singular, the issue about the
        % singularity mentioned above can be ignored.
        
        Ab = [A,b];
        N_rows = size(A,1);
        % forward eliminaton
        for kk = 1:N_rows-1
            % searching for the largest pivot to reduce the round-off error
            % named /emph{partial pivoting}
            maxim = abs(Ab(kk,kk)); r = kk;
            for ii = kk+1:N_rows
                temp = abs(Ab(ii,kk));
                if temp>maxim
                    maxim = temp;
                    r = ii;
                end
            end            
            if (maxim < 1e-7) 
                %either A is singular or the round-off error may be large
                x = nan;
                return;
            end
            if kk~=r
                % exchanging the rows
                Ab([kk r],:) = Ab([r kk],:);                
            end
            
            % elimination of the lower triangle elements of the kth column 
            pivot = Ab(kk,kk); 
            for ii = kk+1:N_rows
                if (Ab(ii,kk)~=0)
                    multiple = Ab(ii,kk)/pivot;
                    Ab(ii,:) = Ab(ii,:) - multiple * Ab(kk,:); 
                end
            end
        end
        
        % backward substitution
        ii = N_rows;        
        if abs(Ab(ii,ii))<1e-7 %checking the last pivot
            %either A is singular or the round-off error may be large
            x = nan;
            return;
        end        
        x = zeros(N_rows,1);
        while (ii>0)
            temp = Ab(ii,1:end-1)*x;
            x(ii) = (Ab(ii,end)-temp)/Ab(ii,ii);
            ii = ii -1;
        end
    end
end

