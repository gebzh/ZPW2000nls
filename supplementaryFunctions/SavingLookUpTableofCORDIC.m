function [] = SavingLookUpTableofCORDIC()
% saving the look-up table of the CORDIC method into CORDIC algorithm data.mat
% file
% note that the amplis_mult is according to the hyperbolic rotation for
% square root, it has no use when i=1
powers = zeros(301,1);
arctan = zeros(301,1);
amplis_mult = zeros(301,1);
sqrts_for_scaling = zeros(11,1);

powers(1) = 1;
arctan(1) = atan(1);
amplis_mult(1) = 1;
temp = 1;
k = 4;
for i = 1:300
    powers(i+1) = powers(i)/2;
    arctan(i+1) = atan(powers(i+1));
    temp = temp / 4;
    amplis_mult(i+1) = amplis_mult(i)*sqrt(1-temp);
    if i==k
        % for ensuring the convergence, 
        % see https://www.mathworks.com/help/fixedpoint/examples/compute-square-root-using-cordic.html
        amplis_mult(i+1) = amplis_mult(i+1)*sqrt(1-temp);
        k = 3*k+1;
    end
end

for i = 1:10
    sqrts_for_scaling(i) = sqrt(i);
end
sqrts_for_scaling(11) = sqrt(0.1);

save('CORDIC algorithm data.mat','powers','arctan','amplis_mult',...
    'sqrts_for_scaling');
end