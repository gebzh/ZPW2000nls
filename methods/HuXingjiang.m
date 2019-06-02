function [ f0,f1 ] = HuXingjiang( in,fs ) %胡幸江方法
%- sub-Nyquist sampling and zero-padding are applied to a windowed ZPW-2000
%  signal
%- when the sampling length M is 2048，the precision of this method is 
%  about +-0.05 
%- for all of the carrier frequency, the downsampling scales are tuned to
%  fit for the sampling frequency 8000 Hz and to make K in Hu2012 always be
%  9 in this program.
%   in: the received signal; fs: the sampling frequency
%   f0: the estimated carrier frequency; f1: the estimated low frequency
    
global Timer;
global Timer1;
Timer = 0;
Timer1 = 0;
tic;

    if (fs~=8000)
        fprintf('错误！采样频率应为8000Hz\n'); 
        return;
    end

    % windowing
    s = HammingWindow(in); 
    
    dsscale = [23 19 17 15]; %the downsampling scale
    
    FreqCalc(s);
    
    function FreqCalc(signal)
        N = length(signal);
        S = fft(signal); Y = conj(S).*S;
        fft_reso = fs/N;

        % coarse estimation
        kl = round([1696,1996,2296,2596] / fft_reso)+1;
        kr = round([1704,2004,2304,2604] / fft_reso)+1;
        maximum = 0;
        for i = 1:4
            [tmp] = max(Y(kl(i):kr(i)));
            if (tmp(1)>maximum)
                maximum = tmp(1);
                num_interest = i;
            end
        end
        
        % down-sampling        
        M = round(1000*fs/dsscale(num_interest)); 
        %in order to bring out the same freq resolution with the grid
        %search of the NLS real model based method.
        sig = zeros(M,1);
        tmp = downsample(signal,dsscale(num_interest));
        sig(1:length(tmp)) = tmp;
        f_s = fs/dsscale(num_interest);
        freq_resolution = f_s/M;
        
        % caculating the exact carrier frequency
        S = fft(sig);      
        spectra = S.*conj(S);        
        [~,kcf] = max(spectra);
        fb_prime = (kcf(1)-1)*freq_resolution;
        f0 = f_s/2 - (fb_prime - 9*f_s/2);

        %fprintf('载频为%12.1f\n',f0);
        right_margin = kcf(1) - round(9 / freq_resolution);
        [~,kfl] = max(spectra(1:max([right_margin,1])));
        left_margin = kcf(1) + round(9 / freq_resolution);
        [~,tmp] = max(spectra(min(left_margin,round(end/2)):round(end/2)));
        kfr = left_margin + tmp(1) - 1;
        f1 = freq_resolution*(kfr-kfl(1))/2;
        %fprintf('低频为%12.1f\n',f1);
    end

Timer = Timer + toc;

end