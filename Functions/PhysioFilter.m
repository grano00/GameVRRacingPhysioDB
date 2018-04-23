function [ newTableOfFile ] = PhysioFilter( tableOfFiles,fps,sr )
%PHYSIOFILTER This function filter the physiological data
%   Detailed explanation goes here




vrInGame = tableOfFiles.VRInGame;
novrInGame = tableOfFiles.NOVRInGame;

%In the following function there is a loop where the data will be selected, filtered, and
%restored.
newInGVR = AllPhysioDataFilter(vrInGame,fps,sr);
newInGNoVR = AllPhysioDataFilter(novrInGame,fps,sr);

%create the new table
newTableOfFile = tableOfFiles;
newTableOfFile.VRFilteredPhysio = newInGVR;
newTableOfFile.NOVRFilteredPhysio = newInGNoVR;

end

function [newInGame] = AllPhysioDataFilter(inGame,fps,sr)
%This function loop for all files and call the filter function.

newInG = cell(size(inGame));
for i = 1:length(inGame)
    data = inGame{i};
    load(data);
    gp = gamePartsTable;
    
    for p=1:height(gp)
        for j=1:width(gp)
            %%Filter part here
            gp{p,j} = DataFilter(gp{p,j},fps,sr,GetPath(data),p,j);
        end
    end
    gamePartsFilt = gp;
    
    newInG{i} = [GetPath(data) 'physioFilt.mat'];
    save(newInG{i},'gamePartsFilt');
end
newInGame = newInG;
end

function [filtData] = DataFilter(gp,fps,sr,path,row,col)

    filtData = gp;
    
    for i=1:length(filtData)
        
       tempData = filtData{i};

       %THE DATA WAS ALREADY CENTERED AND SCALED INTO ExtractInGame
       %Function
       
       %Filt ECG
       ecgFiltLP = getMyFilter('bp',[8,30],[5,35],sr,60);
       tempData(:,2) = myfiltfilt(ecgFiltLP,tempData(:,2));
       
       %Filt EMG
       %see http://gpu.di.unimi.it/docs/iciap-virtual-emg.pdf

       emgFiltHP = getMyFilter('hp',20,19,sr,60);
       
       tempData(:,3:7) = produceEMGData(tempData(:,3:7),sr,emgFiltHP);
       
       %Filt GSR
       [tempData(:,end-3),tonic,phasic] = produceGSRData(tempData(:,end-3),...
           sr,path,...
           num2str(row),num2str(col),num2str(i));
       tempData = [tempData(:,1:end-3) tonic phasic tempData(:,end-2:end)];
       %add the new names at varName
       load([path 'startingData.mat']);
       varNames = [varNames(1:end-3) {'GSRTonic', 'GSRPhasic'} varNames(end-2:end)];
       save([path 'varNames.mat'],'varNames');
       
       %Filt Breath
       coeff = ones(1, sr)/sr;
       tempData(:,end-1) = myfiltfilt(coeff,1,tempData(:,end-1));
       
       properties = GetProperties();
       if(properties.savePlots && col == 1)
        indexes = [2:10,12:length(varNames)-1]; %It defines the interesting columns
        Plotta(tempData(:,indexes),varNames(:,indexes),sr,path,...
            [num2str(row) '-' num2str(col) '-' num2str(i)]);
       end
       
       %reassign the data to return
       filtData{i} = tempData;
    
    end
    
    
end

function [] = Plotta(mydata,names,sr,path,pendixFileName)

    properties = GetProperties();
    
    %find the correct number of rows/cols in the subplot
    p_rows = 0;
    p_cols = 0;
    for i = 2:length(names)
        if (mod(length(names),i) == 0)
            p_cols = i;
            p_rows = length(names) / i;
            break
        end
    end

    
    %Plot filtered data in time domain    
    if(exist([path 'plots\filteredTimeData' pendixFileName '.png'],'file') == 0 || ...
        properties.rewritePlot)
    
        h = figure('units','points','outerposition',[0 0 2560 1080]);
        if(properties.plotSilentMode)
            set(h, 'Visible', 'off');
        end
        PlotTimeFiltered(mydata,names,sr,p_rows,p_cols);
        print([path 'plots\filteredTimeData' pendixFileName],'-dpng');
        savefig([path 'plots\filteredTimeData' pendixFileName]);
        close all
    end
    
    %Plot filtered data in frequency domain
    if(exist([path 'plots\filteredFFTData' pendixFileName '.png'],'file') == 0 || ...
        properties.rewritePlot)
    
        t = figure('units','points','outerposition',[0 0 2560 1080]);
        if(properties.plotSilentMode)
            set(t, 'Visible', 'off');
        end
        PlotFFTFiltered(mydata,names,sr,p_rows,p_cols);
        print([path 'plots\filteredFFTData' pendixFileName],'-dpng');
        savefig([path 'plots\filteredFFTData' pendixFileName]);
        close all
    end

end

function [] = PlotTimeFiltered(mydata,names,sr,row,col)
    time = 1:length(mydata(:,1));
    time = time ./ sr;
    
    for i= 1:length(names)
        subplot(row,col,i);
        plot(time,mydata(:,i));
        title(names{i})
        xlabel('Sec.s');
    end
end

function [] = PlotFFTFiltered(mydata,names,sr,row,col)

    
    for i= 1:length(names)
        subplot(row,col,i);
        [y,f]= myfft(mydata(:,i),sr,-1);
        if(contains(lower(names{i}),'resp') == 1 || contains(lower(names{i}),'gsr') == 1)
           index = f<10;
           f = f(index);
           y = y(index);
        elseif (contains(lower(names{i}),'emg') == 1)
            index = f>18 | f<203;
            f = f(index);
            y = y(index);
            
        end
        plot(f,y);
        title(names{i})
        xlabel('Hz');
    end
end

function [emgs] = produceEMGData(emgs,sr,myFilter)
    [R,C] = size(emgs);
    for i=1:C
       emg = emgs(:,i);
       emg = myfiltfilt(myFilter,emg);
       %emg = abs(emg);
       %[up,down] = envelope(emg);
       %emgs(:,i) = up; 
       emgs(:,i) = emg;
    end
    
end

function [gsr,tonic,phasic] = produceGSRData(gsr,sr,path,row,col,subcol)
    
    save('testgsr1.mat')
    time = ((1:length(gsr)) / 556)';
    %necessary for ledalab
    events = (1:length(gsr))'*0;
    
    tempDir = dir;
    absPath = tempDir.folder;
    newPath = [absPath '/' path];
    
    properties = GetProperties();
    if( exist([newPath 'gsrTable' row '-' col '-' subcol '.mat'],'file') ~= 2 || properties.rewrite )
        tonic = gsr;
        phasic = gsr;

        if(strcmp(col,'1'))
            dlmwrite([newPath 'gsrTable' row '-' col '-' subcol '.txt'],[time gsr events],'delimiter','\t');

            Ledalab([newPath 'gsrTable' row '-' col '-' subcol '.txt'], ...
            'open', 'text',...
            'filter',[1 5],...
            'smooth',{'adapt',6},...
            'analyze','CDA', ...
            'optimize',4);

            load([newPath 'gsrTable' row '-' col '-' subcol '.mat'])

            cd(absPath);

            gsr = data.conductance;
            tonic = analysis.tonicData';
            phasic = analysis.phasicData';

        end
    else
        load([newPath 'gsrTable' row '-' col '-' subcol '.mat'])
        gsr = data.conductance;
        tonic = analysis.tonicData';
        phasic = analysis.phasicData';
    end
    
    save('testgsr2.mat')
end


