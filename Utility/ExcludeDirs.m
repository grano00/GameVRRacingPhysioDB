function [ dirs ] = ExcludeDirs( dirs,dirExclude,mytype )
%EXCLUDEDIRS This function remove from the list the element of dirExclude
%   Detailed explanation goes here

    todelete = [];
    
    for i=1:length(dirExclude)
        element = dirExclude{i};
        
        for j=1:length(dirs)
            if(dirs(j).name == element)
                todelete = [todelete j];
            end
        end
        
    end
if(exist('mytype') > 0 && contains(mytype,'save'))
    ndirs = dirs(todelete);
    dirs = ndirs;
else    
    dirs(todelete) = [];
end
end

