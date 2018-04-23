function [ path ] = GetNameP( p, howmany, until )
%GETPATH Summary of this function goes here
%   Detailed explanation goes here
if(p(end) == '\'), p(end) = []; end

myFolders = split(p,'\');
path = '';


myFolders = myFolders((end-howmany-until+1):(end-until));

path = strjoin(myFolders,'\');
path = [path '\'];
end

