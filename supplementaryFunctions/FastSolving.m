function a = FastSolving(t,h,qL,L,omega_L)
LlpLlT = zeros(3,3,3);
for i = 1:3
    for j = 2:i
        LlpLlT(j-1,j,i) = 1;
    end
    for j = 1:i-1
        LlpLlT(j+1,j,i) = 1;
    end
end

    % Solving the linear system Ax = b with time complexity O(L) 
    % according to properties of the Toeplitz-plus-Hankel system 
    gamma = zeros(L,1);  phi = zeros(L,1);  psi = zeros(L,1); 
    r_bold = zeros(L,1);  rho = zeros(L,1);  beta = zeros(L,1); 
    lambda = zeros(L,1);  a = zeros(L,1); 
    t = [t;0]; h = [h;0]; %for situation that L==1

    r_little = t(0+1) + h(2);
    gamma(1) = 1/r_little;
    phi(1) = qL(1)/r_little;
    psi(1) = 1/r_little;
    r_bold(1) = t(1+1) + h(3);
    rho(1) = -r_bold(1)/r_little;
    beta(1) = rho(1) * gamma(1);
    lambda(1) = omega_L(1);
    a(1) = lambda(1)*gamma(1);
    for ll = 1:L-1 
        r_little = t(0+1) + h(2*ll+2);
        gamma_ll = gamma(1:ll);
        gamma(ll+1) = 1/(r_little+r_bold(1:ll)'*beta(1:ll)/gamma_ll(ll));
        gamma(1:ll) = gamma(ll+1) * beta(1:ll) /gamma_ll(ll);
        a(1:ll+1) = [a(1:ll);0] + gamma(1:ll+1)*(omega_L(ll+1)-r_bold(1:ll)'*a(1:ll));
        if ll == L -1
            break;
        end
        phi(1:ll+1) = [phi(1:ll);0] + (qL(ll+1) - r_bold(1:ll)'*phi(1:ll))*gamma(1:ll+1);
        psi(1:ll+1) = [psi(1:ll);0] - r_bold(1:ll)'*psi(1:ll)*gamma(1:ll+1);
        r_bold(1:ll+1) = t(ll+2:-1:2) + h(ll+3:2*ll+3);
        rho(ll+1) = -r_bold(1:ll+1)'*gamma(1:ll+1);
        beta(1:ll+1) = ((rho(ll+1)-rho(ll))*eye(ll+1)+LlpLlT(1:ll+1,1:ll+1,ll+1))*...
            gamma(1:ll+1) - [gamma_ll;0] + psi(ll+1)*phi(1:ll+1) - phi(ll+1)*...
            psi(1:ll+1);
    end
end