function [ newTableOfFiles ] = ExtractInGame( tableOfFiles,startingGame,endGame,sr,fps)
%EXTRACTINGAME Summary of this function goes here
%   Detailed explanation goes here+
    newTableOfFiles = tableOfFiles;
    data = tableOfFiles.rawData;

    [r c] = size(data);
    newdata = cell(r,c);
    
    for i = 1:r
        for j = 1:c
            %load all files
            load(char(data{i,j}));
           
            %Scale and center
            sac = @(x)(x-mean(x))./std(x); %Use the std for scale
            %sac = @(x)(x-mean(x))./max(abs(x-mean(x))); %Use the canonical normalization
            p = physioData(:,2:end-1);
            physioData(:,2:end-1) = sac(p);
            
            %Get the indexes of starting and ending game session
            startingGameList = find(physioData(:,1) == startingGame);
            stopGameList = find(physioData(:,1) == endGame);
            
            stopGameList = stopGameList - 1;
            
            %if there is no a equal ammount of starting and stopping value
            %means that there is an error..
            if(length(startingGameList) ~= length(stopGameList))
                error('there is a problem on the event list. They are not of the same length.');
            end
            
            % ## extract the different parts and store them in a table
            bef{1} = physioData(1:(startingGameList(1)-1),:);
            aft{1} = physioData((stopGameList(end)+1):end,:);
            
            ing = cell(length(startingGameList),1);
            outg = cell(length(startingGameList) - 1,1);
            
            for k = 1:length(ing)
                ing{k} = physioData(startingGameList(k):stopGameList(k),:);
            end
            
            for k = 1:length(outg)
                outg{k} = physioData(stopGameList(k)+1:startingGameList(k+1)-1,:);
            end
            
            %merge all
            gamePartsTable = table(ing',bef',outg',aft','VariableNames',{'InGame','BeforeGame','BetweenSessions','AfterGame'});
            
            % ##
            
            %Add the labeling at the game part
            
            %Extract the labeling
            labelData = cell(size(ing));
            start = 0;
            counter = 1;
            t = 1;
            %Define start and stop game
            while (t <= length(label(:,4)))
                if(label(t,4) == 1)
                   start = t;
                   while (label(t,4) == 1)
                       t = t + 1;
                   end
                   labelData{counter} = label(start:(t-1),1:3);
                   counter = counter + 1;

                end
                t = t + 1;
            end
            
            
            
            %Save the data
            tempPath = GetPath(char(data{i,j}));
            save([tempPath 'inGameData.mat'], 'gamePartsTable', 'labelData');
            newdata{i,j} = [tempPath 'inGameData.mat'];

            
            
        end
    end
    
    newdata = cell2table(newdata, 'VariableNames',{'VRInGame' 'NOVRInGame'});
    newTableOfFiles = [data,newdata];
    
    

end
