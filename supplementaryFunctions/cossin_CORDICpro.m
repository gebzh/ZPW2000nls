function out = cossin_CORDICpro( theta , k )
    %calculating the value of cos(theta) using CORDIC algorithm
    %   inputs:
    %   theta: the angle value 
    %   k: k=1 refers to cos, while k=2 refers to sin
    %   outputs:
    %   out: the result of cos(theta) or sin(theta)
%========================* load data *==========================
global arctan;
global powers;

if isempty(arctan) || isempty(powers)
    SavingLookUpTableofCORDIC();
    load('CORDIC algorithm data.mat','arctan','powers');
end
%==============* adjusting angle to be within [-pi/2,pi/2] *=============
    theta = theta - 2*pi*round(theta / (2*pi));
    minu = 1;
    if theta > pi/2 
        theta = theta - pi;
        minu = -1;
    end
    if theta < -pi/2
        theta = theta + pi;
        minu = -1;
    end
%========================* CORDIC alg *========================
    z = theta;
    x = 1; y = 0;
    idx = 0;

    while abs(z)>1e-9
        xtmp = x*powers(idx+1); %bitsra(x,idx); % multiply by 2^(-idx)
        ytmp = y*powers(idx+1); %bitsra(y,idx); % multiply by 2^(-idx)
        if z < 0
            x = x + ytmp;
            y = y - xtmp;
            z = z + arctan(idx+1);
        else
            x = x - ytmp;
            y = y + xtmp;
            z = z - arctan(idx+1);
        end
        idx = idx + 1;
    end
    out = -100000;
    if k == 1
        out = minu * x / sqrt(x^2+y^2);
    end
    if k == 2
        out = minu * y / sqrt(x^2+y^2);
    end
end