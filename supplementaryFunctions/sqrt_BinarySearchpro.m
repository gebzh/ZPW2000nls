function [ out ] = sqrt_BinarySearchpro( in )
% calculating the square root of the input using binary search alg
%   input: 
%   in: the input value 
%   output:
%   out: the square root of x as the output
%========================* load data *==========================
global sqrts_for_scaling;

if (isempty(sqrts_for_scaling))
    SavingLookUpTableofCORDIC();
    load('CORDIC algorithm data.mat','sqrts_for_scaling');
end
%========================* scaling *==========================
% this scaling process makes the value of in to be within [1,2)
one_of_nine = ones(size(in));
in_sc = in;
factor_of_ten = zeros(size(in));
out = factor_of_ten;
for i = 1:length(in)
    if in(i)<0
        out(i) = -1;
        continue;
    end
    if in(i) <= 1e-7
        out(i) = 0;
        continue;
    end
    if in(i)>2 || in(i)<1
        % convert the number to the form of scientific notation
        while (in_sc(i)>=10 || in_sc(i)<1)
            if (in_sc(i) >= 10)
                in_sc(i) = in_sc(i) * 0.1;
                factor_of_ten(i) = factor_of_ten(i) + 1;
            end
            if (in_sc(i) < 1)
                in_sc(i) = in_sc(i) * 10;
                factor_of_ten(i) = factor_of_ten(i) - 1;
            end
        end
        % dividing the value with scientific notation with a single digit,
        % resulting in a value within the interval [1,2)
        one_of_nine(i) = floor(in_sc(i));
        in_sc(i) = in_sc(i) / one_of_nine(i);
    end
%========================* binary search alg *==========================
    l = 1; r = 2;
    while (r-l>1e-9)
        m = (l+r)/2;
        if m*m > in_sc(i)
            r = m;
        else
            l = m;
        end
    end
    out(i) = l;
%========================* scaling back *==========================

    if factor_of_ten(i) > 0
        factor = 10;
        idx = 10;
    else
        factor = 0.1;
        factor_of_ten(i) = -factor_of_ten(i);
        idx = 11;
    end
    for j = 1:floor(factor_of_ten(i)/2)
        out(i) = out(i) * factor;
    end
    if mod(factor_of_ten(i),2)
        out(i) = out(i) * sqrts_for_scaling(idx);
    end
    out(i) = out(i) * sqrts_for_scaling(one_of_nine(i));
end
end

