function [Y,f]=myfft(y,Fs,Fc)
%MYFFT    Calcola lo spettro di un segnale
%   [Y,f]=myfft(y,Fs)
%   Y: spettro del segnale
%   f: frequenze
%   y: segnale
%   Fs: frequenza di sampling del segnale
[R,C] = size(y);


for i=1:C
    y2 = y(:,i);
    yy = detrend(y2);
    T = 1/Fs;
    L = length(yy);
    NFFT = 2^nextpow2(L);
    tmp = fft(yy,NFFT)/L;
    Y(:,i) = 2*abs(tmp(1:NFFT/2));

end

    f = Fs/2*linspace(0,1,NFFT/2);

% filtro lo spettro
if(Fc ~= -1)
    T = NFFT/Fs;
    [b,a] = besself(2,2*pi*Fc);
    [b,a] = bilinear(b,a,T);
    Y = filtfilt(b,a,Y);
end


