function [ newTableOfFiles ] = AutoRegFeatureSelection( tableOfFiles,type)
%AUTOFEATURESELECTION Summary of this function goes here
%   Detailed explanation goes here

newTableOfFiles = tableOfFiles;

if(strcmp(type,'mergeGame'))
    
    disp('Start merged game');
    vrFMatrix = tableOfFiles.VRfeaturesMatrix;
    novrFMatrix = tableOfFiles.NOVRfeaturesMatrix;
    vrLabel = tableOfFiles.VRFilteredLabel;
    novrLabel = tableOfFiles.NOVRFilteredLabel;
    
    %In the following function there is a loop where the data will be selected, filtered, resampled, and
    %restored.
    [newGamesFeaturesTable] = FGameLoop(vrFMatrix,novrFMatrix,vrLabel, novrLabel);
    
    newTableOfFiles.GamesFeaturesSel = newGamesFeaturesTable;
    
elseif(strcmp(type,'vr'))
    
    disp('Start vr');
    vrFMatrix = tableOfFiles.VRfeaturesMatrix;
    vrLabel = tableOfFiles.VRFilteredLabel;
    
    %In the following function there is a loop where the data will be selected, filtered, resampled, and
    %restored.
    newVRFeaturesTable = FSelLoop(vrFMatrix,vrLabel);
   
    %create the new table
    newTableOfFiles.VRFeaturesSel = newVRFeaturesTable;
    
elseif(strcmp(type,'novr'))
    
    disp('Start novr');
    novrFMatrix = tableOfFiles.NOVRfeaturesMatrix;
    novrLabel = tableOfFiles.NOVRFilteredLabel;
    
    %In the following function there is a loop where the data will be selected, filtered, resampled, and
    %restored.
    newNOVRFeaturesTable = FSelLoop(novrFMatrix,novrLabel);
    
    
    %create the new table
    newTableOfFiles.NOVRFeaturesSel = newNOVRFeaturesTable;
end


end


function [newTable] = FGameLoop(DataMat1,DataMat2,Label1,Label2)

newTable = cell(size(DataMat1));
for i = 1:length(newTable)
    
    
    load(DataMat1{i});
    load(Label1{i});
    MLMatrixes1 = MLMatrixes;
    secondLabeling1 = secondLabeling;
    
    load(DataMat2{i});
    load(Label2{i});
    MLMatrixes2 = MLMatrixes;
    secondLabeling2 = secondLabeling;
    
    disp(['Start user: ' GetNameP(DataMat1{i},2,1)]);
    
    TAB1 = [MLMatrixes1{1};MLMatrixes2{1}];
    TAB2 = [MLMatrixes1{2};MLMatrixes2{2}];
    LAB1 = [secondLabeling1{1};secondLabeling2{1}];
    LAB2 = [secondLabeling1{2};secondLabeling2{2}];
    
    clear 'secondLabeling1' 'secondLabeling2' 'secondLabeling' ...
        'MLMatrixes1' 'MLMatrixes2' 'MLMatrixes'
    
    FeatureMatrix = cell(1,2);
    
    disp('Start FS 1/2');
    FeatureMatrix{1} = FSelection(TAB1,LAB1);
    disp('Start FS 2/2');
    FeatureMatrix{2} = FSelection(TAB2,LAB2);
    
    disp('Saving');
    newTable{i} = [GetPath(GetPath(DataMat1{i})) 'GamesFeatureSel.mat'];
    save(newTable{i},'FeatureMatrix');
    
end

end

function [newTable] = FSelLoop(DataMat,Labeling)

%Set the new empty path table
newTable = cell(size(DataMat));

for i = 1:length(newTable)
    
    load(DataMat{i});
    load(Labeling{i});
    
    disp(['Start user: ' GetNameP(DataMat{i},2,1)]);
    
    MLMatrixes = JoinTable(MLMatrixes,table());
    secondLabeling = JoinTable(secondLabeling,[]);
    
    FeatureMatrix = cell(size(MLMatrixes));
    
    for p = 1:length(FeatureMatrix)
        disp(['Start FS ' num2str(p) '/' num2str(length(FeatureMatrix))]);
        FeatureMatrix{p} = FSelection(MLMatrixes{p},secondLabeling{p});
        save(['test' num2str(p) 'InFSelLoop.mat']);
    end
    
    disp('Saving...');
    newTable{i} = [GetPath(DataMat{i}) 'FeatureMatrix.mat'];
    save(newTable{i},'FeatureMatrix');
end

end


function [FM] = FSelection(MLMat,lab)

FM = struct;
%Get the correlation of all dataFeatures
% -> Arousal
[corrTableArousal, corrArouNames,arouCorrIndexes] = ...
    GenerateCorrelationTable(MLMat,lab(:,1),0.05);

% -> Valence
[corrTableValence, corrValNames,valCorrIndexes] = ...
    GenerateCorrelationTable(MLMat,lab(:,2),0.05);


fun = @myRandomForestRegressionErr;

funVariables = {'Err' 'Pred' 'TrainY' 'TrainX' 'TestY' 'TestX' 'Fun'};

arousal = lab(:,1);
valence = lab(:,2);
X = table2array(MLMat);
[XR,XC] = size(X);
XArou = X(:,logical(arouCorrIndexes));
XVal = X(:,logical(valCorrIndexes));
seed = 320;
niteration = 1;
cva = cvpartition(arousal,'kfold',5);
cvv = cvpartition(valence,'kfold',5);
log = 1;

disp('Start Arousal at');
dt = datetime
[arouSFBSIndexes,arouHist,arouErr,arouSelected,cvArouInfo] = mysffs( XArou,arousal,fun,...
    log,niteration,cva,funVariables,seed);
disp('End at');
datetime
disp('Duration:');
datetime-dt
disp(['N of feature selected = ' num2str(sum(arouSFBSIndexes))]);
disp('-----------');
disp('Start Valence at');
dt=datetime
[valSFBSIndexes,valHist,valErr,valSelected,cvValInfo] = mysffs( XVal,valence,fun,...
    log,niteration,cvv,funVariables,seed);
disp('End at');
datetime
disp('Duration:');
datetime-dt
disp(['N of feature selected = ' num2str(sum(valSFBSIndexes))]);

myBase = zeros(1,XC);
myBase(logical(arouCorrIndexes)) = arouSFBSIndexes;
arouSFBSIndexes = myBase;

myBase = zeros(1,XC);
myBase(logical(valCorrIndexes)) = valSFBSIndexes;
valSFBSIndexes = myBase;

save('testBefFM.mat');

myrownames = {'pval','pindex','sffsindex'};

arousalFeatures = array2table([corrTableArousal;arouCorrIndexes;...
    arouSFBSIndexes],...
    'VariableNames',corrArouNames,'RowNames',myrownames);
valenceFeatures = array2table([corrTableValence;valCorrIndexes;...
    valSFBSIndexes],...
    'VariableNames',corrValNames,'RowNames',myrownames);

FM.arousalFeatures = arousalFeatures;
FM.arousalCV = cvArouInfo;
FM.valenceFeatures = valenceFeatures;
FM.valenceCV = cvValInfo;

save('testAftFM.mat');
end


%Get the correlation table for each features
function [pcorrelation,names,indexes] = GenerateCorrelationTable(MLMat,targ,bound)
pcorrelation = zeros(1,width(MLMat));

for i = 1:width(MLMat)
    [~,pcorrelation(i)] = corr(MLMat{:,i},targ);
end
indexes = zeros(1,width(MLMat));
indexes(pcorrelation<=bound) = 1;
names = MLMat.Properties.VariableNames;

end


function newM = JoinTable(M,emptyMat)
newM = M;
for i = 1:length(M)
    emptyMat = [emptyMat; M{i}];
end
newM{i+1} = emptyMat;

end


%Functions for sfs/sbs/sffs
function [err] = myCrossValidation(XT,yT,xt,yt)
X = [XT;xt];
Y = [yT;yt];
fun = @myRandomForestRegressionErr;
FoldCVP = cvpartition(length(Y),'KFold',10);
rescv = crossval(fun,X,Y,'partition',resubCVP)/...
    resubCVP.TestSize;
err = rescv;
end


function [ret] = myRandomForestRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
tb = TreeBagger(100,XT,yT,'Method','regression','Options',ops);
pr = predict(tb,xt);
err = sqrt(immse(yt,pr));
ret = {{err},{pr},{yT},{XT},{yt},{xt},{tb}};
end


function [err] = mySVMQRegressionErr(XT,yT,xt,yt)
%ops = statset('UseParallel',true);
mysvm = fitrsvm(XT,yT,'KernelFunction','gaussian');
pr = predict(mysvm,xt);
err = sqrt(immse(yt,pr));
end
