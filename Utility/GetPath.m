function [ path ] = GetPath( p )
%GETPATH Summary of this function goes here
%   Detailed explanation goes here
if(p(end) == '\'), p(end) = []; end
myFolders = split(p,'\');
path = '';

%for i=1:(length(myFolders)-1)
%    path = ([path myFolders(i) '/']);
%end
myFolders(end) = [];
path = strjoin(myFolders,'\');
path = [path '\'];
end

