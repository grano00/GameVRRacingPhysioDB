function [ M ] = cell2deepmat( mycell )
%CELL2NUM Summary of this function goes here
%   Detailed explanation goes here
M = [];
for i=1:length(mycell)
    t = mycell{i};
    for p=1:length(t)
        M = [M; t{p}];
    end
end

end

