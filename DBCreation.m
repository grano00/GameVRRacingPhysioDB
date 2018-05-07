%This function is used in order to crate a Readable DB
%It was structured in directories with the player ID with subdirectories
%with the game name
%For each participant will be created  a Mat file that cointains the
%labeling data, the physiological data, and a table with the combination of
%both
%Futhermore, we store these three variables as csv files with the name
%labeling.csv, physio.csv, physioWithLabeling.csv

clear all
close all

%The the function available in the subdir
addpath(genpath('Functions'));
addpath(genpath('Utility'));

dirNameDB = 'ReadableDataSet/';
mkdir(dirNameDB);

load('./MatlabDataset/Examples/tableOfFiles.mat');
inVR = tableOfFiles.VRInGame;
noVR = tableOfFiles.NOVRInGame;

for i = 1:length(inVR)
    %Get User ID
    id = GetNameP(inVr{1},1,2);
    path = [dirNameDB id];
    mkdir(path);
    
    %NOVR
    load(noVR{i});
    
    %Data of PCs
    subpath = [path 'PCs/']
    [LABEL,PHYSIO,MERGED] = extractInfo(labelData,gamePartsTable.InGame,1);
    save([subpath 'dataset.mat'],'LABEL','PHYSIO','MERGED');
    writetable(LABEL,[subpath 'label.csv'],'Delimiter',',');
    writetable(PHYSIO,[subpath 'physio.csv'],'Delimiter',',');
    writetable(MERGED,[subpath 'physioWithLabeling.csv'],'Delimiter',',');
    
    %Data of RO
    subpath = [path 'RO/']
    load(noVR{i});
    [LABEL,PHYSIO,MERGED] = extractInfo(labelData,gamePartsTable.InGame,2);
    save([subpath 'dataset.mat'],'LABEL','PHYSIO','MERGED');
    writetable(LABEL,[subpath 'label.csv'],'Delimiter',',');
    writetable(PHYSIO,[subpath 'physio.csv'],'Delimiter',',');
    writetable(MERGED,[subpath 'physioWithLabeling.csv'],'Delimiter',',');
    
    
    %VR
    load(inVR{i});
    
    %Data of PCs
    subpath = [path 'PCsInVR/']
    [LABEL,PHYSIO,MERGED] = extractInfo(labelData,gamePartsTable.InGame,1);
    save([subpath 'dataset.mat'],'LABEL','PHYSIO','MERGED');
    writetable(LABEL,[subpath 'label.csv'],'Delimiter',',');
    writetable(PHYSIO,[subpath 'physio.csv'],'Delimiter',',');
    writetable(MERGED,[subpath 'physioWithLabeling.csv'],'Delimiter',',');
    
    %Data of RO
    subpath = [path 'ROInVR/']
    load(noVR{i});
    [LABEL,PHYSIO,MERGED] = extractInfo(labelData,gamePartsTable.InGame,2);
    save([subpath 'dataset.mat'],'LABEL','PHYSIO','MERGED');
    writetable(LABEL,[subpath 'label.csv'],'Delimiter',',');
    writetable(PHYSIO,[subpath 'physio.csv'],'Delimiter',',');
    writetable(MERGED,[subpath 'physioWithLabeling.csv'],'Delimiter',',');
    
end

function [L,P,T] = extractInfo(lab, gp,which)

l = lab{which};
p = gp{which};

%Remove Frame Counter
l = l(:,[2,3]);

%Transform in table
L = table(l(:,1),l(:,2),'VariableNames',{'Arousal' 'Valence'});

%Remove inGame id
p = p(:,2:end);
names = {'ECG', 'EMG_Zygomatic', 'EMG_Supercilii', 'EMG_Nasalis', 'EMG_OrbicularisOculi', 'EMG_Temporalis', ...
    'GSR', 'Light', 'Respiration', 'SampleRate'};
P = array2table(p,'VariableNames',names);

%Interp label 
N = length(p(:,1));
l = ResampleTo(l,N);

T = P;
T.Arousal = l(:,1);
T.Valence = l(:,2);
    

end
