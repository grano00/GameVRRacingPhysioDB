%#####################%
%    GENERAL NOTES    %
%######################

%START SETUP! 
clear all
close all

%The the function available in the subdir
addpath(genpath('Functions'));
addpath(genpath('Utility'));
addpath(genpath('Ledalab')); %GSR TOOLBOX


%directory = 'Examples';
directory = 'Real';
targetDir = 'MatlabDataset'; %directory where are stored the matlab files

tabOfFilDir = [targetDir '/' directory];
%Sample rate of physiological data
sr = 556;
%Frame per Second of labeling (Sample Rate)
fps = 60;

%read the list of the ids in the directory and return it in an array
dirs = dir(['Data/' directory '/']);
dirs = dirs(3:end,:); %delete dirs of root and superroot

dirExclude = {'15'};

%Remove the directories that must be excluded
dirs = ExcludeDirs(dirs,dirExclude);


%Display Parameters
properties = GetProperties();
fprintf([ '\n--------------------------\n'...
    'Parameters:\nRewrite = ' num2str(properties.rewrite) ...
    '\nSavePlots = ' num2str(properties.savePlots)...
    '\nHidePlots = ' num2str(properties.plotSilentMode)...
    '\nRewritePlots = ' num2str(properties.rewritePlot)...
    '\n--------------------------\n']);

%END SETUP!
try
    parpool;
catch ME
    disp(ME);
end

%load all the dataset and store them as matlab files
disp('Start to read the csv files');
tableOfFiles = loadcsvs(dirs,targetDir,directory,fps,sr);
savetable(tableOfFiles,tabOfFilDir);


%Split the games into in-game and non-ingame
disp('Split into ingame and not ingame');
tableOfFiles = ExtractInGame(tableOfFiles,14,15,sr,fps);
savetable(tableOfFiles,tabOfFilDir);

%filter and fix data of the labeling data
disp('Filt the labeling data');
tableOfFiles = LabelFilterAndScale(tableOfFiles,fps,sr);
savetable(tableOfFiles,tabOfFilDir);

%filter and fix data of physio
disp('Filt the physiological data');
tableOfFiles = PhysioFilter(tableOfFiles,fps,sr);
savetable(tableOfFiles,tabOfFilDir);

%extract the features at second precision
disp('Get the data features');
tableOfFiles = FeaturesExtraction(tableOfFiles,fps,sr);
savetable(tableOfFiles,tabOfFilDir);

%create the table for ML algos
disp('Create the tables for ML algos');
tableOfFiles = MatrixGenerator(tableOfFiles,sr);
savetable(tableOfFiles,tabOfFilDir);


%MACHINE LEARNING

%##############%
%  REGRESSION  %
%##############%

%Automaticaly Feature Selection & Prediction

warning('off','all');


%VR
disp('Applying Feature Selection Algorithm IN VR');
tableOfFiles = AutoRegFeatureSelection(tableOfFiles,'vr');
savetable(tableOfFiles,tabOfFilDir);

%NOVR
disp('Applying Feature Selection Algorithm IN NOVR');
tableOfFiles = AutoRegFeatureSelection(tableOfFiles,'novr');
savetable(tableOfFiles,tabOfFilDir);

%BY GAME
disp('Applying Feature Selection Algorithm BY GAME');
tableOfFiles = AutoRegFeatureSelection(tableOfFiles,'mergeGame');
savetable(tableOfFiles,tabOfFilDir);

warning('on','all');


disp('Using ML Algo for predict the data');
tableOfFiles = RegTrainAndPred(tableOfFiles,tabOfFilDir);
savetable(tableOfFiles,tabOfFilDir);

disp('Create boxplot and excel files');
tableOfFiles = produceResults(tableOfFiles,tabOfFilDir);
savetable(tableOfFiles,tabOfFilDir);

%%%%%%%%%%%%%%%%%%%%%%%%
%   MERGE REGRESSION   %
%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------
%Here will be merged the regression data and analyize them for general
%propouse and prediction
% -
% The results give the prediction of: 
%   ALL, G1, G2, VR, NOVR, 
%   G1VR, G2VR, G1NOVR, G2NOVR
%-----------------------

%Collapse the file and return a new table of file where in the row there
%are 9 line, one for each type

%Apply the same function of FeatureSelection

%Apply the same function of Regression and plot



%##################%
%  CLASSIFICATION  %
%##################%
