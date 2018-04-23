function [Y,f,YLog]=mypsd(data,Fs)
[R,C] = size(data);

Y = zeros(length(0:Fs/R:Fs/2),C);

N = R;
for i = 1:C
    x = data(:,i);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    Y(:,i) = psdx;
    YLog(:,i)= 10*log10(psdx);
end

f= 0:Fs/length(x):Fs/2;