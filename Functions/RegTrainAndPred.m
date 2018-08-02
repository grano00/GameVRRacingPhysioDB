function [ newTableOfFiles ] = RegTrainAndPred( tableOfFiles,subrootPath )
%REGTRAINANDPRED This function collect the data with the selected features
%and train them with a ML Function/5 cv.
%   The results will be achieved into results subpath. In this folder will
%   be find the plot and a 4 text file = "*id*(VAL/AROU)TrainData.csv" that
%   cointains the final data of prediction/target; "*id*(VAL/AROU)Infos.txt" with
%   the info of the results (e.g. the missregression indexes).
%   Furthermore it create "FinalIndexes.mat" into subrootPath/indexes with the
%   index ranking, and the connected plots
%   Lastly, return the tableOfFiles with added the reference of
%   "RegResults.mat" with the results.


%TODO FINIRE FINALINDEXES + AGGIUNGERE INFO E LEGGENDA AL PLOT + MODIFICARE
%KENEL GAUSSIAN REGRESSION

newTableOfFiles = tableOfFiles;
subrootPath = [subrootPath '/indexes/'];
mkdir(subrootPath);


%Disable the warning of directory created
mywarning = 'MATLAB:MKDIR:DirectoryExists';
warning('off',mywarning)

%In the following function there is a loop where the data will be selected,
%filtered, resampled, and restored.
[newTableOfFiles.VRbyUserResults,arouTab1,valTab1]= ...
    TrainLoop(tableOfFiles.VRFeaturesSel,{'VRGame1','VRGame2','VR'});
[newTableOfFiles.NOVRbyUserResults,arouTab2,valTab2] = ...
    TrainLoop(tableOfFiles.NOVRFeaturesSel,{'NOVRGame1','NOVRGame2','NOVR'});
[newTableOfFiles.GAMESbyUserResults,arouTab3,valTab3] = ...
    TrainLoop(tableOfFiles.GamesFeaturesSel,{'Game1','Game2'});

warning('on',mywarning)

fun = @collapse;
ArousalIndexes = varfun(fun,[arouTab1,arouTab2,arouTab3]);
ValenceIndexes = varfun(fun,[valTab1,valTab2,valTab3]);

save([subrootPath 'FinalIndexes.mat'],'ArousalIndexes','ValenceIndexes');

plotBar(ArousalIndexes,subrootPath,'arousal');
plotBar(ValenceIndexes,subrootPath,'valence');




end

function [results,arouIndex,valIndex] = TrainLoop(mylist,names)
results = cell(size(mylist));


%Preset the indexes
load(mylist{1});

arouIndex = cell(length(mylist),length(FeatureMatrix));
valIndex = cell(length(mylist),length(FeatureMatrix));

clearvars FeatureMatrix

for i = 1:length(mylist)
    [results{i},arouIndex(i,:),valIndex(i,:)] = ...
        MyTraining(mylist{i},names);
end

arouIndex = cell2table(arouIndex,'VariableNames',names);
valIndex = cell2table(valIndex,'VariableNames',names);


end

function [npath,arouIndexes,valIndexes] = MyTraining(data,names)

load(data);

id = GetNameP(data,1,2);
id(end)  = [];

rpath = [GetPath(data) 'Results\'];
mkdir(rpath);

if(length(names) ~= length(FeatureMatrix))
    error(['Error in Mytraining of ' GetPath(data) ...
        '.. The lengths of var names and FeatureMatrix do not fit.']);
end


arouIndexes = cell(1,length(FeatureMatrix));
valIndexes = cell(1,length(FeatureMatrix));
res = cell(1,length(FeatureMatrix));

properties = GetProperties();
cellarray = cell(size(FeatureMatrix));
if(exist([GetPath(data) 'RegResults.mat']) > 0 && ...
        properties.rewrite == false)
    disp([id ' was already done']);
    npath = [GetPath(data) 'RegResults.mat'];
    for i=1:length(FeatureMatrix)
        content = FeatureMatrix{i};
        
        
        %Get arousal and valence Indexes
        arouIndexes{1,i} = content.arousalFeatures(3,:);
        valIndexes{1,i} = content.valenceFeatures(3,:);
    end
else
    for i=1:length(FeatureMatrix)
        content = FeatureMatrix{i};
        
        
        %Get arousal and valence Indexes
        arouIndexes{1,i} = content.arousalFeatures(3,:);
        valIndexes{1,i} = content.valenceFeatures(3,:);
        
        arousal = struct;
        valence = struct;
        
        
        %Get the data of feature selection
        AT = content.arousalCV{end,end}{1,1}{1,1};
        VT = content.valenceCV{end,end}{1,1}{1,1};
        
        %Restore the data
        [XA,YA] = restoration(AT);
        [XV,YV] = restoration(VT);
        
        %CV 5Fold
        cva = cvpartition(YA,'kfold',5);
        cvv = cvpartition(YV,'kfold',5);
        
        %ML Function
        %fun = @myRandomForestRegressionErr;
        fun = @myGaussianProcessRegression;
        
        %Define the parallel methods
        opts = statset('UseParallel',true);
        
        %Apply CV
        disp(['ID: ' id]);
        dt1 = datetime;
        disp(['Start Arousal CV of ' names{i}]);
        RA = crossval(fun,XA,YA,...
            'partition',cva,'Options',opts);
        dt2 = datetime;
        disp(['Start Valence CV of ' names{i}]);
        RV = crossval(fun,XV,YV,...
            'partition',cvv,'Options',opts);
        dt3 = datetime;
        
        arousal.trainTime = seconds(dt2-dt1);
        valence.trainTime = seconds(dt3-dt2);
        
        %Transofrm in table
        funVariables = {'Err' 'Pred' 'TrainY' 'TrainX' 'TestY' 'TestX' 'Fun'};
        RA = cell2table(RA,'VariableNames',funVariables);
        RV = cell2table(RV,'VariableNames',funVariables);
        
        %Get the info of the CV
        [arousal.pval,arousal.r,arousal.rmse,...
            arousal.mse,arousal.mae,...
            arousal.pred,arousal.testy, arousal.ML] = extractInfo(RA,cva);
        [valence.pval,valence.r,valence.rmse,...
            valence.mse,valence.mae,...
            valence.pred,valence.testy, valence.ML] = extractInfo(RV,cvv);
        
        %Prepare and Plot the data
        DrawPlotAndResults(arousal,rpath,id,'AROU','arousal',...
            ['Arousal-' names{i}],{[0.5,0,0],[0.9,0,0]})
        DrawPlotAndResults(valence,rpath,id,'VAL','valence',...
            ['Valence-' names{i}],{[0,0.5,0],[0,0.9,0]})
        
        
        result = content;
        result.arousalFin = arousal;
        result.valenceFin = valence;
        res{i} = result;
        
    end
    Results = cell2table(res,'VariableNames',names);
    npath = [GetPath(data) 'RegResults.mat'];
    save(npath, 'Results');
    
end
end

function [] = DrawPlotAndResults(mystruct,path,id,type,...
    ylab,plotName,colors)

properties = GetProperties();

%Plot the data
if(properties.savePlots)
    close all;
    h = figure('units','points','outerposition',[0 0 2560 1080]);
    if(properties.plotSilentMode)
        set(h, 'Visible', 'off');
    end
    plot(mystruct.testy,'color',colors{1},'LineWidth', 2);
    hold on
    plot(mystruct.pred,'o','color',colors{2});
    plot(mystruct.pred,'+','color',colors{2});
    title([id ' ' plotName]);
    xlabel('s');
    ylabel(ylab);
    
    print([path id '-' plotName 'points'],'-dpng');
    savefig([path id '-' plotName 'points']);
    
    plot(mystruct.pred,'color',colors{2});
    
    print([path id '-' plotName 'line'],'-dpng');
    savefig([path id '-' plotName 'line']);
    close all;
end

%Prepare TrainData.csv
tab = [mystruct.pred,mystruct.testy];
csvwrite([path id '-' plotName 'TrainData.csv'],tab);

%Prepare info.txt
mywar = 'MATLAB:DELETE:FileNotFound';
warning('off',mywar);
delete([path id '-' plotName 'Infos.txt']);
warning('on',mywar);

fileID = fopen([path id '-' plotName 'Infos.txt'],'w');
fprintf(fileID,'ID: %s\n',id);
fprintf(fileID,'Name %s\n\n',plotName);
fprintf(fileID,'R: %2.10f with p-val of corr %2.6f\n',...
    mystruct.r,mystruct.pval);
rsquare = mystruct.r^2;
fprintf(fileID,'R-square: %2.10f\n',rsquare);
fprintf(fileID,'RMSE: %2.10f\n',mystruct.rmse);
fprintf(fileID,'MSE: %2.10f\n',mystruct.mse);
fprintf(fileID,'MAE: %2.10f\n',mystruct.mae);
MLTab = mystruct.ML;
funs = MLTab.Fun;
err = MLTab.Err;
test = funs{1};
[~,C] = size(test.X);
[R,~] = size(MLTab.TestIdx{1});
fprintf(fileID,'N* of Features %i\n',C);
fprintf(fileID,'N* of instances %i\n',R);
fprintf(fileID,'\n\n\n Model Info:\nGaussian process regression (GPR) model\n');
fprintf(fileID,'PARAMETER:\n');
fprintf(fileID,'KernelFunction: %s\n',test.KernelFunction);
fprintf(fileID,'BasisFunction: %s\n',test.BasisFunction);

for i=1:length(funs)
    mymodel = funs{i};
    [R,~] = size(mystruct.ML.TestX{i});
    fprintf(fileID,'----------------------\n');
    fprintf(fileID,'Fold %i\nSigma: %2.10f\n',i,mymodel.Sigma);
    fprintf(fileID,'Beta: %2.10f\n',mymodel.Beta);
    fprintf(fileID,'N* of test %i\n',R);
    fprintf(fileID,'Fold RMSE %2.10f\n',err(i));
end
fclose(fileID);


end

function [pval,r,rmse,mymse,mymae,pred,testy,R] = extractInfo(R,cv)
y = R.TestY;
pr = R.Pred;


pred = zeros(length(test(cv,1)),1);
testy = zeros(length(test(cv,1)),1);
idx = cell(length(y),1);

for i = 1:length(y)
    idx{i} = test(cv,i);
    pred(idx{i},1) = pr{i};
    testy(idx{i},1) = y{i};
end

R.TestIdx = idx;
rmse = RMSE(pred,testy);
mymse = mse(pred,testy);
mymae = mae(pred,testy);
[r,pval] = corrcoef(pred,testy);
r = r(1,2);
pval = pval(1,2);


end

function [X,Y] = restoration(T)
ty = T.TestY;
tx = T.TestX;
idx = T.TestIdx;

[R,C] = size(tx{1});
X = zeros(length(idx{1}),C);
Y = zeros(length(idx{1}),1);

for i = 1:length(idx)
    X(idx{i},:) = tx{i};
    Y(idx{i},:) = ty{i};
end

end


function [ret] = myGaussianProcessRegression(XT,yT,xt,yt)
ops = struct('UseParallel',true);
md = fitrgp(XT,yT,'KernelFunction','ardrationalquadratic',...
    'OptimizeHyperparameters','auto');
%'HyperparameterOptimizationOptions',ops);
pr = predict(md,xt);
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,md};
end

function [ret] = myRandomForestRegressionErr(XT,yT,xt,yt)
ops = statset('UseParallel',true);
tb = TreeBagger(128,XT,yT,'Method','regression','Options',ops);
pr = predict(tb,xt);
err = sqrt(immse(yt,pr));
ret = {err,pr,yT,XT,yt,xt,tb};
end

function [xxx] = collapse(x)
xx = table();
for i = 1:length(x)
    myx = x{i};
    myx.Properties.RowNames = {};
    xx = [xx;myx];
end
names = xx.Properties.VariableNames;
p = table2array(xx);
p = sum(p,1);

xxx = array2table(p,'VariableNames',names);

end


function [] = plotBar(myindexes,subrootPath,type)


properties = GetProperties();

if(properties.savePlots)
    
    nsubrootPath = [subrootPath type '/'];
    mkdir(nsubrootPath);
    
    indNames = myindexes.Properties.VariableNames;
    rename = @(x) strrep(x,'collapse_','');
    indNames = cellfun(rename,indNames,'UniformOutput',false);
    
    for i = 1:width(myindexes)
        tab = myindexes{1,i};
        tabNames = tab.Properties.VariableNames;
        rename2 = @(x) strrep(x,'_','-');
        tabNames = cellfun(rename2,tabNames,'UniformOutput',false);
        myname = indNames{i};
        mydata = table2array(tab);
        
        %Get the data not empty
        toconsider = find(mydata ~= 0);
        tabNames = tabNames(toconsider);
        mydata = mydata(toconsider);
        
        %plot and save
        close all;
        f = figure('units','points','outerposition',[0 0 2560 1080]);
        if(properties.plotSilentMode)
            set(f, 'Visible', 'off');
        end
        bar(categorical(tabNames),mydata);
        ylabel('Selection Number');
        mytitle = [type '-' myname 'Bars'];
        title(mytitle);
        print([nsubrootPath mytitle],'-dpng');
        savefig([nsubrootPath mytitle]);
    end
end
end