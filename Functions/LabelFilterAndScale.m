function [ newTableOfFile ] = LabelFilterAndScale( tableOfFiles,fps,sr )
%LABELFILTERANDSCALE This function loop for each user and filter the
%labeling data
%   Detailed explanation goes here

vrInGame = tableOfFiles.VRInGame;
novrInGame = tableOfFiles.NOVRInGame;

%In the following function there is a loop where the data will be selected, filtered, resampled, and
%restored.
newInGVR = AllLabelDataFilter(vrInGame,fps,sr);
newInGNoVR = AllLabelDataFilter(novrInGame,fps,sr);

%create the new table
newTableOfFile = tableOfFiles;
newTableOfFile.VRFilteredLabel = newInGVR;
newTableOfFile.NOVRFilteredLabel = newInGNoVR;

end

function [newInGame] = AllLabelDataFilter(inGame,fps,sr)
%This function loop for all files and call the filter function.

newInG = cell(size(inGame));
for i = 1:length(inGame)
    data = inGame{i};
    load(data);
    gpig = gamePartsTable.InGame;
    labelDataFilt = cell(size(labelData));
    secondLabeling = cell(size(labelData));
    for j=1:length(labelData)
        [labelDataFilt{j}, secondLabeling{j}]= labelingFilter(labelData{j},fps,sr,length(gpig{j}));
    end
    newInG{i} = [GetPath(data) 'labelFilt.mat'];
    save(newInG{i},'labelDataFilt','secondLabeling');
end
newInGame = newInG;
end

function [filtData, secondData] = labelingFilter(M,fps,sr,newLength)

    %Scale and center
    filtData = (M(:,2:3) - 50) / 50;
   
 
    %Length Resample
    filtData = ResampleTo(filtData,newLength);    
    
    %Moving Average Filter
    coeff = ones(1, sr)/sr;
    filtData(:,1) = filtfilt(coeff, 1, filtData(:,1));
    filtData(:,2) = filtfilt(coeff, 1, filtData(:,2));
    
    secondData = ResampleTo(filtData,newLength/sr);
    
end