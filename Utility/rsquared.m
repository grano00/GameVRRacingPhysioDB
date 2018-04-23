function [ r2,rss,tss ] = rsquared( x,y )
%RSQUARED Summary of this function goes here
%   Detailed explanation goes here
    rss = sum((x-mean(x)).^2);
    tss = sum((x - y).^2);
    
    r2 = 1-(rss/tss);
    
end

