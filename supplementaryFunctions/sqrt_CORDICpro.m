function [ out ] = sqrt_CORDICpro( in )
% calculating the square root of the input using CORDIC alg with hyperbolic 
% rotation in vectoring mode
% the hyperbolic rotation means the vector (x,y) is rotated along with a
% hyperbolic funtion, such that the length of the vector is varying
% comparing with the elementary rotation whose vector rotates along with
% circle
%   input: 
%   in: the input value 
%   output:
%   out: the square root of x as the output
%========================* load data *==========================
global amplis_mult;
global sqrts_for_scaling;
global powers;

if (isempty(amplis_mult) || isempty(sqrts_for_scaling) || isempty(powers))
    SavingLookUpTableofCORDIC();
    load('CORDIC algorithm data.mat','amplis_mult','sqrts_for_scaling',...
        'powers');
end
%========================* scaling *==========================
% refer to https://stackoverflow.com/questions/18646105/cordic-for-square-roots
% when in is out of the range [0.5,2), the alg cannot converage, which has been
% confirmed by simulaiton. Hence here a scaling process is introduced.
one_of_nine = ones(size(in));
in_sc = in;
factor_of_ten = zeros(size(in));
out = factor_of_ten;
for i = 1:length(in)
    if in(i)>2 || in(i)<0.5
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


%========================* CORDIC alg *==========================
    x = in_sc(i) + 0.25; 
    y = in_sc(i) - 0.25; 
    % the angle is within -pi/4 and pi/4, thus the rotation angle can be atan(2^(-1)) 
    % in the first iteration
    idx = 1;
    k = 4;
    while (abs(y)>1e-5)
        xtmp = x*powers(idx+1); %bitsra(x,idx); % multiply by 2^(-idx)
        ytmp = y*powers(idx+1); %bitsra(y,idx); % multiply by 2^(-idx)
        if y < 0
            x = x + ytmp;
            y = y + xtmp;
        else
            x = x - ytmp;
            y = y - xtmp;
        end
        if idx == k
            xtmp = bitsra(x,idx); % multiply by 2^(-idx)
            ytmp = bitsra(y,idx); % multiply by 2^(-idx)
            if y < 0
                x = x + ytmp;
                y = y + xtmp;
            else
                x = x - ytmp;
                y = y - xtmp;
            end
            k = 3*k+1;
        end
        idx = idx + 1;
    end
    out(i) = x/amplis_mult(idx);
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

