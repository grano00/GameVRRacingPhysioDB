function [  ] = savetable( tableOfFiles, dir )
%SAVETABLE Summary of this function goes here
%   Detailed explanation goes here
save([dir '/tableOfFiles.mat'], 'tableOfFiles');

end

