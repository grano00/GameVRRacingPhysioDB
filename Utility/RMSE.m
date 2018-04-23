function [ rmse ] = RMSE( a, b )
%RMSE calculate the rmse between a and b
%   Detailed explanation goes here

err = immse(a,b);
rmse = sqrt(err);

return;
end