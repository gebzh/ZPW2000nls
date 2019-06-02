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