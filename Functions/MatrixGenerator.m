function [ newTableOfFile ] = MatrixGenerator( tableOfFiles,sr )
%MATRIXGENERATOR This function read the data and create a matrix for ML
%application.
%   It will store for each experiental session n different matrix with features
%   of physiological data(one for each game) and nx4 different matrix of
%   labeling and game event classification


vrFeatures = tableOfFiles.VRFeatures;
novrFeatures = tableOfFiles.NOVRFeatures;


%In the following function there is a loop where the data will be selected,
%filtered, and restored.
newVRFMatrix = FMGeneratorLoop(vrFeatures,sr);
newNOVRFMatrix = FMGeneratorLoop(novrFeatures,sr);


%create the new table
newTableOfFile = tableOfFiles;
newTableOfFile.VRfeaturesMatrix = newVRFMatrix;
newTableOfFile.NOVRfeaturesMatrix = newNOVRFMatrix;

end


%###########################################%
%# FUNCTIONS FOR FEATURES MATRIX GENERATOR #%
%###########################################%

%This function loop for each user and load the files
function [mytable] = FMGeneratorLoop(features,sr)

%Set the new empty path table
mytable = cell(size(features));
newMLM = cell(size(features));

for i = 1:length(features)
    data = features{i};
    load(data);
    
    matrixes = cell(1,length(GameFeatures));
    for p=1:length(GameFeatures)
        
        matrixes{p} = FMGenerator(GameFeatures{p},sr);
        
    end
    
    MLMatrixes = matrixes;
    
    newMLM{i} = [GetPath(data) 'MLMatrix.mat'];
    save(newMLM{i},'MLMatrixes');
end

%Return the table of the features file path
mytable = newMLM;

end


function [MATRIX] = FMGenerator(startingTables,sr)

zygo = getEmgMatrix(startingTables.EMG_Zygomatic,sr,'zygo_');
supercilii = getEmgMatrix(startingTables.EMG_Supercilii,sr,'supercilii_');
nasalis = getEmgMatrix(startingTables.EMG_Nasalis,sr,'nasalis_');
orbicularis = getEmgMatrix(startingTables.EMG_OrbicularisOculi,sr,'orbicularis_');
temporalis = getEmgMatrix(startingTables.EMG_Temporalis,sr,'temporalis_');
gsr = getGSRMatrix(startingTables.GSR,sr,'GSR_');
gsrTonic = getGSRMatrix(startingTables.GSRTonic,sr,'GSRTonic_');
gsrPhasic = getGSRMatrix(startingTables.GSRPhasic,sr,'GSRPhasic_');
respiration = getRespMatrix(startingTables.Respiration,sr,'RESP_');
hr = getHRMatrix(startingTables.HeartRate,sr,'HR_');


MATRIX = [zygo, supercilii, nasalis, orbicularis, temporalis, ...
    gsr, gsrTonic, gsrPhasic, respiration, hr];

end

%This function return a matrix starting from a table of EMG features
%It will also select the informative features that will store in the matrix

%THE STARTING FEATURES ARE:
%-----------------------------------------------------------------
% 'Data[TABLE OF SRx1 WITH THE DATA]','InterpData',...
% 'FFT'[TABLE OF FREQ WITH IN THE 1TH COL THERE ARE HZ],...
% 'PSD'[TABLE OF FREQ WITH IN THE 1TH COL THERE ARE HZ AND IN THE LAST
% COL THERE ARE THE VALUE IN 10LOG10 FORM],...
% 'AR[SEE FFT]','BandPower','Power',...
% 'Integrated', 'Mean','MeanABS','ModMeanABS1','ModMeanABS2',...
% 'MAVS','RMS','VAR','WL','ZC','SSC',...
% 'WAMP','SSI','FMD','FMN','MFMD','MFMN','FR'...
%-----------------------------------------------------------------
function M = getEmgMatrix(mtab,sr,name)
%Remove useless features
mtab.Data = [];
mtab.FFT = [];
mtab.PSD = [];
mtab.AR = [];

%Tranform in uniform matrix
mtab = tab2extendedTable(mtab);

%Change the variable names
mtab.Properties.VariableNames = strcat(name,mtab.Properties.VariableNames);


M = mtab;
end

%This function return a matrix starting from a table of HR features
%It will also select the informative features that will store in the matrix

%THE STARTING FEATURES ARE:
%-----------------------------------------------------------------
% 'hrMean','hrInterp'
%-----------------------------------------------------------------
function M = getHRMatrix(mtab,sr,name)
%Remove useless features
mtab.hrInterp = [];


%Change the variable names
mtab.Properties.VariableNames = strcat(name,mtab.Properties.VariableNames);

M = mtab;
end

%This function return a matrix starting from a table of RESPIRATION features
%It will also select the informative features that will store in the matrix

%THE STARTING FEATURES ARE:
%-----------------------------------------------------------------
% 'Data[TABLE OF SRx1 WITH THE DATA]','InterpData',...
% 'FFT'[TABLE OF FREQ WITH IN THE 1TH COL THERE ARE HZ],...
% 'PSD'[TABLE OF FREQ WITH IN THE 1TH COL THERE ARE HZ AND IN THE LAST
% COL THERE ARE THE VALUE IN 10LOG10 FORM],...
% 'AR[SEE FFT]','BandPower','Power',...
% 'Integrated', 'Mean','MeanABS','ModMeanABS1','ModMeanABS2',...
% 'MAVS','RMS','VAR','WL','ZC','SSC',...
% 'WAMP','SSI','FMD','FMN','MFMD','MFMN','FR'...
% 'RespirationRate{'respMean','respInterp'}'
%-----------------------------------------------------------------
function M = getRespMatrix(mtab,sr,name)
%Remove useless features
mtab.Data = [];
mtab.FFT = [];
mtab.PSD = [];
mtab.AR = [];
mtab.ZC = [];
mtab.SSC = [];
mtab.WAMP = [];
mtab.Mean = [];

%Select the mean insted than interpolation
resp = mtab.RespirationRate;
resp.respInterp = [];
mtab.RespirationRate = table2array(resp);


%Tranform in uniform matrix
mtab = tab2extendedTable(mtab);

%Change the variable names
mtab.Properties.VariableNames = strcat(name,mtab.Properties.VariableNames);

M = mtab;
end

%This function return a matrix starting from a table of GSR features
%It will also select the informative features that will store in the matrix

%THE STARTING FEATURES ARE:
%-----------------------------------------------------------------
% 'Data[TABLE OF SRx1 WITH THE DATA]','InterpData',...
% 'FFT'[TABLE OF FREQ WITH IN THE 1TH COL THERE ARE HZ],...
% 'PSD'[TABLE OF FREQ WITH IN THE 1TH COL THERE ARE HZ AND IN THE LAST
% COL THERE ARE THE VALUE IN 10LOG10 FORM],...
% 'AR[SEE FFT]','BandPower','Power',...
% 'Integrated', 'Mean','MeanABS','ModMeanABS1','ModMeanABS2',...
% 'MAVS','RMS','VAR','WL','ZC','SSC',...
% 'WAMP','SSI','FMD','FMN','MFMD','MFMN','FR'...
%-----------------------------------------------------------------
function M = getGSRMatrix(mtab,sr,name)
%Remove useless features
mtab.Data = [];
mtab.FFT = [];
mtab.PSD = [];
mtab.AR = [];
mtab.ZC = [];
mtab.SSC = [];
mtab.WAMP = [];


%Tranform in uniform matrix
mtab = tab2extendedTable(mtab);

%Change the variable names
mtab.Properties.VariableNames = strcat(name,mtab.Properties.VariableNames);

M = mtab;
end



%This function analyze the data of the table, check if there are variables
%with the matrix content and tranform it
function M = tab2extendedTable(tab)

indexToDelete = [];

for i = 1:width(tab)
    try
        
        arr = table2array(tab{1,i});
        %In the table, the first col are the frequencies meantime the
        %other col is the information
        
        %Get the variable names and replace the special characters
        varNames = strcat([tab(1,i).Properties.VariableNames],...
            cellstr(num2str(arr(:,1))));
        
        varNames = varNames';
        varNames = strrep(varNames,'.','p');
        varNames = strrep(varNames,' ','');
        
        %Generate the rows
        M = zeros(height(tab),length(varNames));
        for k = 1:height(tab)
            arr = table2array(tab{k,i});
            arr = arr';
            M(k,:) = arr(2,:);
        end
        
        %Create new table and merge
        mytab = cell2table(num2cell(M),'VariableNames',varNames);
        
        tab = [tab mytab];
        indexToDelete = [indexToDelete i];
    catch ME
        
    end
end


%Remove the variables
tab(:,indexToDelete) = [];

M = tab;

end
