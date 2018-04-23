function [ y ] = ResampleTo( array, obj_sample, functionType, interpType )
%RESAMPLETO Summary of this function goes here
%   Array is the array to resample
%   obj_sample is obj value to reach
%   functionType define the type of function to use: 
%       'r' in order to use resample
%       'd' in order to use interp1 and decimal -> Not recommended
%       otherwise it uses interp1

[p,q] = rat(obj_sample / length(array));

if(p/q == 1)
%    error('the sample rate is yet the same');
end

x =  linspace(1, length(array), obj_sample);

if(nargin == 2)
     y = interp1( 1:length(array), array,x);
elseif(nargin == 4)
    y = interp1( 1:length(array), array,x,interpType);
elseif(functionType == 'r')
   %Use Resample 
    y = resample(array,p,q);
elseif(functionType == 'd')
   %Use Interp and Decimation
   if(p/q < 1)
       %decimal
       
       %interp to prox multiple value
       m = ceil(q/p);
       nextVal = obj_sample * m;
       y = interp1(1:length(array),array,linspace(1, length(array), nextVal));
       y = decimate(y,m);
   else
       %interp1
       y = interp1( 1:length(array), array, x);
   end
else
    %%Simple interp1
     y = interp1( 1:length(array), array,x);
end

end