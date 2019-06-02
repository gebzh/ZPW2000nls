function b = CalcThreeOrderNLSCostFunc(omega0,omega1,s)
    N = length(s);
    n0 = -(N-1)/2;
    n = (n0:n0+N-1)';
    Z_l = [cos((omega0-omega1)*n),cos(omega0*n),...
         cos((omega0+omega1)*n),sin((omega0-omega1)*n),...
         sin(omega0*n),sin((omega0+omega1)*n)]; 
    omega_l = Z_l'*s;
    [t,h,qL] = load_data(3,Z_l);
    a_l = FastSolving(t,h,qL,3,omega_l(1:3));
    qL = [-t(2);-t(3)+h(2) ];
    b_l = FastSolving(-t,h,qL,3,omega_l(4:6));
    J = omega_l'*[a_l;-b_l];
    b = [a_l;-b_l];
    A_fl = sqrt(a_l(1)^2 + b_l(1)^2); 
    A_fr = sqrt(a_l(3)^2 + b_l(3)^2); 
    % With the criterion: the both side-frequency components should
    % relatively have the same amplitude, the ZPW-2000 signal can be
    % distinguished from noises. Thus here a comparation is adopted  
    % with respect to the relative amplitude difference of the both 
    % side-freqs 
    if abs(A_fl-A_fr)/A_fl > 0.2
        J = J;
    end
    
    function [t,h,qL] = load_data(i,Z_l)
        if i==3
            t = zeros(3,1); h = zeros(6,1); 
            ZZ = Z_l'*Z_l;
            t = ((ZZ(1,1:3)+ZZ(4,4:6))/2)';
            h(2:4) =  ((ZZ(1,1:3)-ZZ(4,4:6))/2)';
            h(5:6) = ((ZZ(2:3,3) - ZZ(5:6,6))/2)';
%             t(1) = N/2; omg = [omega1; 2*omega1];
%             t(2:3) = cos(omg*(n0+(N-1)/2)).*sin(omg*N/2)./sin(omg/2)/2;
%             omg = 2*omega0-2*omega1:omega1:2*omega0+2*omega1;
%             h(2:6) = cos(omg*(n0+(N-1)/2)).*sin(omg*N/2)./sin(omg/2)/2;
            qL = [t(2); t(3)+h(2) ];
        end
    end
end