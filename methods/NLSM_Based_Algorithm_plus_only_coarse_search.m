function [pc, qc, n_all_ZPW, p_all, q_all, Amplis_all,cos_correlation_all] ...
    = NLSM_Based_Algorithm_plus_only_coarse_search( x_bold, fs )
% NLS method based ZPW-2000 signal demodulation method

N = length(x_bold);
Apslowerbound = 0.15; 
% the lower bound of the amplitude of the carrier freq and the first-side
% freqs components

%========================* loading data *==========================
global Z_AI;
global ZTZ_AI;
global ZTZ_AII;
global SupplementaryParas_AII;
global lengthOfSignal;
global Amplis_std;
global f0_std;
global f1_std;
global mies;

if (isempty(Z_AI) || isempty(ZTZ_AI) || isempty(ZTZ_AII) || ...
        isempty(SupplementaryParas_AII) || isempty(Amplis_std) ||  ...
        lengthOfSignal ~= N)
    Load_Data_plus(N,fs);
end

% definition of the standard carrier freqs, the standard low freqs and the
% modulation indexes 
if isempty(f0_std) || isempty(f1_std) || isempty(mies)
    f0_std = [1701.4; 1698.7; 2001.4; 1998.7; 2301.4; 2298.7; 2601.4; 2598.7];
    f1_std = (10.3:1.1:29)';
    mies = 11./f1_std;
end

% starting the timer
global Timer;
Timer = 0;
tic;
%========================* coarse estimation *==========================
n_all_ZPW = 0;
p_all = zeros(10,1);
q_all = zeros(10,1);
Amplis_all = zeros(5,10);
cos_correlation_all = zeros(10,1);

Jmax = 0; pc = -1; qc = -1;
for p = 1:8
    for q = 1:18
        [J, alph_hat, cc] = CostCalcForCoarseEst(p,q);
        if (J > Jmax)
            Jmax = J;
            pc = p; qc = q;
        end
        if (J > 0)
            n_all_ZPW = n_all_ZPW + 1;
            p_all(n_all_ZPW) = p;
            q_all(n_all_ZPW) = q;
            Amplis_all(:,n_all_ZPW) = sqrt_BinarySearchpro(alph_hat(1:5).^2 ...
                + alph_hat(6:10).^2);
            cos_correlation_all(n_all_ZPW) = cc;
        end
    end
end

% record the time cost
Timer = Timer + toc;
fprintf('time consuming%fs\n',Timer);

%% ========================* functions *==========================
    function [J, lph_hat, c_c] = CostCalcForCoarseEst(pp,qq)
        % If the Amplitudes of the three main freq components are larger than
        % Aslowerbound, the corresponding carrier freq and the low freq will 
        % be choosed as candidates, and the corresponding alpha_hat and
        % the r square will be recorded as references for fault diagnosis
        eta = Z_AI(:,:,pp,qq)'*x_bold;
        lph_hat = GaussianElimination(ZTZ_AI(:,:,pp,qq),eta);
        J = eta'*lph_hat;
        Apls_hat = sqrt_BinarySearchpro(lph_hat(1:5).^2 + lph_hat(6:10).^2);
        c_c = CalcCosCorrelation(Apls_hat,Amplis_std(:,qq));
        if ( Apls_hat(2)<Apslowerbound || Apls_hat(3)<Apslowerbound...
                || Apls_hat(4)<Apslowerbound || (qq <= 9 ...
                && Apls_hat(1)>Apls_hat(2) && Apls_hat(5)>Apls_hat(4)) )
            J = -1;
            c_c = 0;
        end
    end
    function c_c = CalcCosCorrelation(y,x)
        c_c = x'*y/sqrt_BinarySearchpro(sum(x.^2))...
            /sqrt_BinarySearchpro(sum(y.^2));
        c_c = 1-2*arccos_CORDICpro(c_c)/pi;
    end
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
