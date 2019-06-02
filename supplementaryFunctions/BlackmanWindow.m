function [ s ] = BlackmanWindow( signal )
    N1 = length(signal);
    n1 = (0:(N1-1))';
    win = 0.42 - 0.5*cos(2*pi*n1/(N1-1)) + 0.08*cos(4*pi*n1/(N1-1));
    tmp = sum(win);
    s = win .*signal;
end

