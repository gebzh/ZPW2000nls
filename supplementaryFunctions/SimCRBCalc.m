function [ CRBf0, CRBf1 ] = SimCRBCalc( A,fs,N,f1,sigmaSquare )
%Calculating exact CRB for the simplified model
%   inputs:
%   A: the amplitude of the noise&interference-free ZPW-2000 signal
%   fs: the sampling freq
%   N: the length of the signal 
%   f1: the low freq
%   sigmaSquare: the variance of the white Gaussian noise
%   outputs:
%   CRBf0: the CRB of the carrier freq
%   CRBf1: the CRB of the low freq

m = 11/f1;
Asquare = zeros(5,1);
Asquare(1) = (2*m*A*sin(m*pi/2)/((m^2-4)*pi))^2;
Asquare(2) = (2*m*A*cos(m*pi/2)/((m^2-1)*pi))^2;
Asquare(3) = (2*A*sin(m*pi/2)/(m*pi))^2;
Asquare(4) = Asquare(2);
Asquare(5) = Asquare(1);

CRBf0 = 6*sigmaSquare*(fs^2)/((pi^2)*Asquare(3)*(N^3));
CRBf1 = 3*sigmaSquare*(fs^2)/(2*(pi^2)*Asquare(4)*(N^3))...
    + 3*sigmaSquare*(fs^2)/(2*(pi^2)*Asquare(2)*(N^3));

end

