function [T] = loadcsvs( dirs, targetDir,directory,fps,sr )
%LOADCSVS This function load all the csv files and store them into a new
%directory (targetDir). It return the list of files where in the first col 
%there are the ID, in the second col are the vr data and in the tirth 
%are the novr data.


list = table(length(dirs),3);
ids = cell(length(dirs),1);
vr = cell(length(dirs),1);
novr = cell(length(dirs),1);

newdir = [targetDir '\' directory];%defines the new directory
mkdir(newdir); %create the dir

varNames = {'event','ECG','EMG_Zygomatic','EMG_Supercilii',...
            'EMG_Nasalis','EMG_OrbicularisOculi','EMG_Temporalis','GSR',...
            'LUX','Respiration','SampleRate'};

for i = 1:length(dirs)
    ids{i} = dirs(i).name;
    display(ids{i});
    properties = GetProperties();
    
    %VR
    if( exist([newdir '\' dirs(i).name '\vr\startingData.mat'],'file') ~= 2 || properties.rewrite )
        %read the physiological data
        physioData = csvread([dirs(i).folder '\' dirs(i).name '\vr\physio.csv']);
        %read the label value
        label = csvread([dirs(i).folder '\' dirs(i).name '\vr\labeling.csv']);
        
        %fix the sample rate indexes
        physioData(:,end) = physioData(:,end) + 1;

        %check if there are some data that should be recovered
        physioData = recoverData(physioData,sr,dirs(i).name,'vr');

        
        label(:,1) = mod((label(:,1) + 1),fps);


        %save them in matlab file in the new dir
        mkdir([newdir '\' dirs(i).name '\vr\'])
        save([newdir '\' dirs(i).name '\vr\startingData.mat'], 'physioData', 'label','varNames');
    end 
    
    %save the plots
    if(properties.savePlots)
      mkdir([newdir '\' dirs(i).name '\vr\plots']);
      load([newdir '\' dirs(i).name '\vr\startingData.mat']);
      PlotPrint(physioData,label,dirs(i).name,newdir,'vr\plots','rawData');
      st = find(physioData(:,1) == 14);
      PlotPrint(physioData((st(1)- 500):end,:), ...
        label,dirs(i).name,newdir,'vr\plots','inGame');
    end
    vr{i} = [newdir '\' dirs(i).name '\vr\startingData.mat'];

    
    %NOVR
    if( exist([newdir '\' dirs(i).name '\novr\startingData.mat'],'file') ~= 2 || properties.rewrite )
        %read the physiological data
        physioData = csvread([dirs(i).folder '\' dirs(i).name '\novr\physio.csv']);
        %read the label value
        label = csvread([dirs(i).folder '\' dirs(i).name '\novr\labeling.csv']);
        
        %fix the sample rate indexes
        physioData(:,end) = physioData(:,end) + 1;
        
        %check if there are some data that should be recovered
        physioData = recoverData(physioData,sr,dirs(i).name,'novr');
        
        label(:,1) = mod((label(:,1) + 1),fps);
        
        
        %save them in matlab file in the new dir
        mkdir([newdir '\' dirs(i).name '\novr\'])
        save([newdir '\' dirs(i).name '\novr\startingData.mat'], 'physioData', 'label','varNames');
    end 
    
    %save the plots
    if(properties.savePlots)
        mkdir([newdir '\' dirs(i).name '\novr\plots']);
        load([newdir '\' dirs(i).name '\novr\startingData.mat']);
        PlotPrint(physioData,label,dirs(i).name,newdir,'novr\plots','rawData');
        st = find(physioData(:,1) == 14);
        PlotPrint(physioData((st(1)- 500):end,:), ...
            label,dirs(i).name,newdir,'novr\plots','inGame');
    end
    novr{i} = [newdir '\' dirs(i).name '\novr\startingData.mat'];
    

end
    
T = table(ids,'VariableNames',{'ID'},'RowNames',ids);
T.rawData = table(vr,novr,'VariableNames',{'VR' 'NOVR'});
end


%It recovers the losted data
function [P] = recoverData(physio,sr,name,type)
    startingPoints = find(physio(:,end) == 1);
    endingPoints = startingPoints(2:end) - 1;
    endingPoints = [endingPoints; length(physio)];
    
    diffPoints = endingPoints - startingPoints;
    
    P = [];
    [R,C] = size(physio);
    
    for i = 1:length(diffPoints)
       if(diffPoints(i) ~= (sr-1))
           disp(['Error in ' name ' / ' type ' at sec ' num2str(i)]);
           subphysio = physio(startingPoints(i):endingPoints(i),:);
           newData = zeros(sr,C);
           for p = 2:(C-1)
                newData(:,p) = ResampleTo(subphysio(:,p),...
                    sr);
           end
           newData(:,end) = 1:sr;
           newData(:,1) = ones(sr,1)*16;
           if(length(find(subphysio(:,1) == 14))> 0)
               newData(1,1) = 14;
           elseif (length(find(subphysio(:,1) == 15))> 0)
               newData(1,1) = 15;
           end

           P = [P; newData];
       else
           P = [P; physio(startingPoints(i):endingPoints(i),:)];
       end
    end
    
    if(length(P) ~= R)
        disp(['User ' name ' has some error in ' type ' that were fixed']);
    end
end

function [] = PlotPrint(physio,label,name,newdir,type,filename)
    
    starting = find(physio(:,1)==14);
    ending = find(physio(:,1)==15);
    properties = GetProperties();
    if(exist([newdir '\' name '\' type '\' filename '.png'],'file') == 0 || ...
            properties.rewritePlot)
        physio(:,[1,9,11]) = [];

        titles = {'ECG','EMG1','EMG2','EMG3','EMG4','EMG5', ...
            'GSR','Resp','VAL','AUR'};
        %figure('units','normalized','outerposition',[0 0 1 1]);
        h = figure('units','points','outerposition',[0 0 2560 1080]);
        if(properties.plotSilentMode)
            set(h, 'Visible', 'off');
        end
        for i = 1:(length(titles)-2)
            subplot(5,2,i);
            plottaPhysio(starting,ending,physio(:,i),titles{i});
        end

        i = i + 1;
        subplot(5,2,i);
        plottaLabel(label(:,4)*100,label(:,2),titles{i});
        i = i + 1;
        subplot(5,2,i);
        plottaLabel(label(:,4)*100,label(:,3),titles{i});


        %save and close the plots!
        print([newdir '\' name '\' type '\' filename],'-dpng');
        savefig([newdir '\' name '\' type '\' filename]);
    end
    close all
end

function [] = plottaPhysio(starting,ending,data,mytitle)

    for i = 1:length(starting)
        rectangle('Position',[starting(i),min(data), ...
            ending(i) - starting(i),max(data) - min(data)], ...
            'FaceColor',[.8 .8 .8]);
        hold on;
    end
    plot(data);
    title(mytitle)
end

function [] = plottaLabel(ingame,data,mytitle)

    plot(ingame,'*');
    hold on;
    plot(data);
    title(mytitle)
end
