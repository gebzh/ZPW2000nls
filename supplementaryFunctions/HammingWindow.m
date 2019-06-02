function s = HammingWindow(signal)
    N1 = length(signal);
    n1 = 0:(N1-1);
    win = 0.54-0.46*cos(2*pi*n1/(N1-1));
    tmp = sum(win);
    s = win'.*signal;
end