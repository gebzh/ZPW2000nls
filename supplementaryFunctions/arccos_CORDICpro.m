function out = arccos_CORDICpro( r )
    %calculating the value of arccos(r) using CORDIC algorithm
    %   inputs:
    %   r: the cosine value
    %   outputs:
    %   out: the result of arccos(theta) 
%========================* load data *==========================
global arctan;
global powers;

if isempty(arctan) || isempty(powers)
    SavingLookUpTableofCORDIC();
    load('CORDIC algorithm data.mat','arctan','powers');
end

%==============* adjusting input value to be within [0,1] *=============
if (r<0)
    r = -r;
    out_pluspie = 1;
else
    out_pluspie = 0;
end
if (r>1 && abs(1-r)<1e-5)
    r = 1;
end
%========================* CORDIC alg *========================
    z = 0;
    x = r; y = sqrt_BinarySearchpro(1-r^2);
    idx = 0;

    while abs(y)>1e-9
        xtmp = x*powers(idx+1); %bitsra(x,idx); % multiply by 2^(-idx)
        ytmp = y*powers(idx+1); %bitsra(y,idx); % multiply by 2^(-idx)
        if y < 0
            x = x - ytmp;
            y = y + xtmp;
            z = z - arctan(idx+1);
        else
            x = x + ytmp;
            y = y - xtmp;
            z = z + arctan(idx+1);
        end
        idx = idx + 1;
    end
    out = z;
    if out_pluspie
        out = out + pi;        
    end
end