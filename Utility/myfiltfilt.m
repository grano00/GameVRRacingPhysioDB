%This function try to solve the problem of filter a small ammout of data
function [data] = myfiltfilt(myfilter,var1,var2)
    
    if(isscalar(var1))
        data = var2;
        try
            data = filtfilt(myfilter,var1,var2);
        catch ME
        end
    else
        data = var1;
        try
            data = filtfilt(myfilter,data);
        catch ME
        end
    end
end