function [ Y,f ] = PlotMyFFT( y,Fs,Fc,plotType )
%PLOTMYFFT Summary of this function goes here
%   Detailed explanation goes here

if (~exist('plotType','var'))
    plotType = 'plot';
end


[Y,f] = myfft(y,Fs,Fc);

if(strcmp(plotType,'stem'))
    stem(f,Y);
else
    plot(f,Y);
end

end

