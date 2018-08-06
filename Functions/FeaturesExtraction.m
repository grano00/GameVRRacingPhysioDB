function [ newTableOfFile ] = FeaturesExtraction( tableOfFiles,fps,sr )
%FEATURESEXTRACTION This function extract the features for each data in
%table
%   Detailed explanation goes here


vrInGame = tableOfFiles.VRFilteredPhysio;
novrInGame = tableOfFiles.NOVRFilteredPhysio;

%In the following function there is a loop where the data will be selected, filtered, resampled, and
%restored.
newFInGVR = FExtractLoop(vrInGame,fps,sr);
newFInGNoVR = FExtractLoop(novrInGame,fps,sr);

%create the new table
newTableOfFile = tableOfFiles;
newTableOfFile.VRFeatures = newFInGVR;
newTableOfFile.NOVRFeatures = newFInGNoVR;

end

function [newInGame] = FExtractLoop(inGame,fps,sr)

%Set the new empty path table
newInG = cell(size(inGame));

for i = 1:length(inGame)
    data = inGame{i};
    load(data);
    load([GetPath(data) 'varNames.mat']);
    %gp will cointain the table with all physio data
    gp = gamePartsFilt;
    gfcell = cell(size(gamePartsFilt.InGame));
    

    
    
    for p=1:height(gp)

        gfcell(p,:) = FExtract(gp(p,:),varNames,fps,sr);
        
    end
    GameFeatures = gfcell;
    
    newInG{i} = [GetPath(data) 'Features.mat'];
    save(newInG{i},'GameFeatures');
end

%Return the table of the features file path
newInGame = newInG;

end

function [inGameCells] = FExtract(gp,varNames,fps,sr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The data in GP will be structured in one line.
% - The first column "InGame" will be splitted into n matrixes. Each matrix contains
% the data of the game session. (i.e. InGame(1) -> Proj. Cars, InGame(2) ->
% RedOut).
% - The second "BeforeGame" and last "AfterGame" columns contain only one
% matrix.
% - The tirth "BetweenSession" column cointains n-1 matrixes that
% rappresent the data between the game session.
%--------------------------------------------------------------------------
%%Here will be extracted the features at the precison of 1 sec (sr samples).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create the interesting array of cell
ingame = gp.InGame;
nogame = [gp.BeforeGame, gp.BetweenSessions, gp.AfterGame];
inGameCells = cell(1,length(ingame));

parfor i=1:length(ingame)

    %Get the cell of approximation window
    windowSize = 3; %It defines the multiplicator of sr for define the window size
    approxWindow = GetApproximationWindow(ingame{i},nogame{i},sr,windowSize);
    
    
    %Get the data features from approximation window data
    %The features are selected considering the paper:
    %------------------------------------------------------------
    %Stages for Developing Control Systems using EMG and EEG Signals: 
    %A survey - Ericka Janet Rechy, Ramirez, and Huosheng Hu
    %------------------------------------------------------------
    
    tableOfFeatures = generalExtractProcess(approxWindow,ingame{:,i},varNames,sr);
    
    
    
    %#####################%
    %#  OTHER FEATURES   #%
    %#####################%
    
    
    
    %Get the hr and rr variability
    adc = (3.3/(2^12))*(10^6);
    floatingWindow = 4;
    [hr, rr,sechr] =  getHeartRate(ingame{i}(:,2),nogame{i}(:,2),nogame{i+1}(:,2),...
        sr,adc, floatingWindow);
    
    tableOfFeatures.HeartRate = sechr;
    
    %Get the respiration features
    [respirationRate,secRespirationRate ]= getBreathFeatures(ingame{i}(:,end-1),nogame{i}(:,end-1),...
        nogame{i+1}(:,end-1), sr);
    
    tableOfFeatures.Respiration.RespirationRate = secRespirationRate;
    %####################%
    
    
    inGameCells{i} = tableOfFeatures;
    
end


end

function [windowCellArray] = GetApproximationWindow(ingame,before,...
    singleWindow,multiplier)

%this function returns a cell array of the approximation window.
%the data start from the end of nogame until the last of ingame.


%Get the ammount of data that must collected by before
ammountBefore = singleWindow * (multiplier);

%Create the new matrix of data
mydata = [before((end-(ammountBefore-1)):end,:) ; ingame];

%initialize the output
windowCellArray = cell(1,(length(ingame) / singleWindow));

for i = (singleWindow*multiplier):singleWindow:length(mydata)
    index = (i / (singleWindow)) - (multiplier-1); %It defines the cell index
    windowCellArray{index} = mydata(((i-(singleWindow*multiplier)+1)):i,:); %Split the data
end

end

function [hr,rr,sechr] = getHeartRate(mydata,before,after,...
    sr,adc,floatingWindow)


y = [before;mydata;after];


%Define the area where there are mydata data
startingMyData = length(before) + 1;
endingMyData = length(y) - length(after);



%Get the rr
[fiducial_points, BL]=embdac(y,sr,adc);

points = fiducial_points(:,1);
points(isnan(points)) = [];

difference = diff(points);


%Find the missing R data due of signal noise
rMiss = round(difference ./ mean(difference)) - 1;

newPoints = points;
for i=1:length(rMiss)
    if(rMiss(i) > 0)
        val = points(i+1) - points(i);
        temp = [];
        for j=1:(rMiss(i)+1)
            temp = [temp; (points(i)+(val*(j/(rMiss(i)+1))))];
        end
        temp = round(temp);
        prev = find(newPoints <= points(i));
        after = find(newPoints >= points(i+1));
        
        newPoints = [newPoints(prev); temp; newPoints(after)];
    elseif(rMiss(i) < 0)
        newPoints(newPoints == points(i)) = [];
    end
    
end

newPoints = unique(newPoints);

%If there are no enough data in after, it create a dummy point!
dummyPos = newPoints + (newPoints(end) - newPoints(end-1));
if(dummyPos > length(y))
   dummyPos = length(y) - 1;
end
newPoints = [newPoints; dummyPos];        


%Get HR from RR
k = diff(newPoints);
hr = k./sr;
hr = 60./hr;

%Approx the hr
coeff = ones(1, floatingWindow)/floatingWindow;
hr = filtfilt(coeff,1,hr);

%Select only mydata area

hr = ResampleTo(hr,length(newPoints(1):newPoints(end)),'i');

startingMyHr = startingMyData - newPoints(1);
endingMyHr = endingMyData - newPoints(1);

hr = hr(startingMyHr:endingMyHr);

hr = hr';

%Transform in second domain
%---
%Area Mean
%---
hrMean = zeros(length(hr)/sr,1);
idx = 1;
for k = 1:sr:length(hr)
    hrMean(idx) = mean(hr(k:(k+sr-1)));
    idx = idx + 1;
end

%---
%Interpolation
%---
hrInterp = ResampleTo(hr,length(hr)/sr,'i');

sechr = table(hrMean,hrInterp','VariableNames',{'hrMean' 'hrInterp'});


rr = newPoints;

rr(rr<startingMyData) = [];
rr(rr>endingMyData) = [];
rr = rr - startingMyData;

end

function [respirationRate,secRespirationRate] = getBreathFeatures(mydata,before,after,sr)

allRespiration = [before;mydata;after];

%Get the pos of all peaks (upper and lower)
[pks1, pos1] = findpeaks(allRespiration,'MinPeakDistance',sr/4);
[pks2, pos2] = findpeaks(-allRespiration,'MinPeakDistance',sr/4);
up = ones(length(pos1),1);
down = (ones(length(pos2),1))*-1;
pos = [[pos1;pos2],[up;down]];
[~,index] = sort(pos(:,1));
pos = pos(index,:);

rate = diff(pos(:,1));

resampledRate = ResampleTo(rate,length(allRespiration),'i');
resampledRate(1:length(before)) = [];
resampledRate(end-(length(after)-1):end) = [];

%return the respiration rate in sec.
respirationRate = (resampledRate/sr)';

%Transform in second domain
%---
%Area Mean
%---
respMean = zeros(length(respirationRate)/sr,1);
idx = 1;
for k = 1:sr:length(respirationRate)
    respMean(idx) = mean(respirationRate(k:(k+sr-1)));
    idx = idx + 1;
end

%---
%Interpolation
%---
respInterp = ResampleTo(respMean,length(respirationRate)/sr,'i');


secRespirationRate = table(respMean,respInterp','VariableNames',{'respMean' 'respInterp'});


end


function [mytable] = generalExtractProcess(appWindow,gameData,varNames,sr)
    %Create the empty table that will collect all the features
    %In each row will be stored the features of appWindow
    %In each col will be stored the information to the type of data
    
    %Remove event, LUX, and sample rate
    indexToRemove = [1,2,11,13];
    names = varNames;
    names(indexToRemove) = [];
    
    %Generate the emtpy table and set the names of the variable into the table
    mytable =  cell2table(cell(0,length(names)), ...
        'VariableNames', names);
    
    %Allocate a variable with length appwindow - 2 because the first and
    %the last element are not drawed for ingame area
    featuresContent = cell(length(appWindow)-2,length(names));
    
    %Exctract the features from all columns
    %NB: the first and the last element of appWindow are outsite of ingame
    %part. It will be used only for extract information.
    for i = 2:length(appWindow)
        mydata = appWindow{i};
        before = appWindow{i-1};
        
        %Remove the un-used indexes
        mydata(:,indexToRemove) = [];
        before(:,indexToRemove) = [];
        
        %Get the iterpolated of raw data
        interpData = ResampleTo(gameData,length(appWindow)-1);
        interpData(:,indexToRemove) = [];
        
        %get the table of features for each type
        for c = 1:length(names)
            
            %Check if it is an emg for additional features extraction
            if(contains(lower(names{c}),'emg') == 1)
                featuresContent{i-1,c} = getFeatures(mydata(:,c),...
                    before(:,(c)),interpData(i-1,c),sr,'emg');
            elseif(contains(lower(names{c}),'resp') == 1)
                featuresContent{i-1,c} = getFeatures(mydata(:,c),...
                    before(:,(c)),interpData(i-1),sr,'resp');
            elseif(contains(lower(names{c}),'gsr') == 1)
                featuresContent{i-1,c} = getFeatures(mydata(:,c),...
                    before(:,(c)),interpData(i-1),sr,'gsr');
            else
                featuresContent{i-1,c} = getFeatures(mydata(:,c),...
                    before(:,(c)),interpData(i-1),sr,'nada');
            end
        end
    end
    
    
    %Assign the data to the table
    mytable = [mytable;featuresContent];
end


function [featuresTable] = getFeatures(mydata,before,uniqueData,sr,signaltype)
        
        idx = 1;% it is the index of the current features

        featuresNames = {'Data','InterpData',...
            'FFT','PSD','AR','BandPower','Power',...
            'Integrated', 'Mean','MeanABS','PMean',...
            'ModMeanABS1','ModMeanABS2',...
            'MAVS','RMS','VAR','WL','ZC','SSC',...
            'WAMP','SSI','FMD','FMN','MFMD','MFMN','FR'...
            };
        featuresCell = cell(1,length(featuresNames));
        
        %Get the data
        featuresCell{idx} = mydata(end-sr+1:end);
        idx = idx + 1;
        
        %Get the interpolated data
        featuresCell{idx} = uniqueData;
        idx = idx + 1;
        
        
        %Get the FFT
        [y,f] = myfft(mydata,sr,-1);
        featuresCell{idx} = [f',y];
        idx = idx + 1;
        
        %Get the PSD
        [y,f] = mypsd(mydata,sr);
        psd = [f',y];
        featuresCell{idx} = psd;
        idx = idx + 1;
        
        %Get PSD in AR using Yule-Walker method
        %- NB : Check for automatically order selection 
        %-> https://it.mathworks.com/help/signal/ug/ar-order-selection-with-partial-autocorrelation-sequence.html
        order = 10;
        dimension = sr * 2;
        [y,f] = pyulear(mydata,order,dimension,sr);
        ar = [f,y];
        featuresCell{idx} = ar;
        idx = idx + 1;
        
        %Get the BandPower of the signal. It is the average power in a
        %freq. range. 

        [bands,steps,mymod,lfreq,hfreq] = getBandParameters(sr,signaltype);
        
        mpow = zeros((((bands - mymod) / steps)+1),2);
        k = 1;

        for i = 0:steps:(bands-steps-mymod)
            mpow(k,:) = [i, bandpower(mydata,sr,[i (i+steps)])];
            k = k + 1;
        end
        if(mymod > 0)
            mpow(k,:) = [i+steps, bandpower(mydata,sr,[i+steps bands])];
        else
            mpow(end,:) = [];
        end
        featuresCell{idx} = mpow;
        idx = idx + 1;
        
        %Get the total power of the signal
        featuresCell{idx} = bandpower(mydata);
        idx = idx + 1;
        
        absData = abs(mydata);
        absBefore = abs(before);
        
        %Get the Integral (area)
        featuresCell{idx} = sum(absData);
        idx = idx + 1;
        
        %Get the mean
        featuresCell{idx} = mean(mydata);
        idx = idx + 1;
        
        %Get the mean of the absoulte values
        featuresCell{idx} = mean(absData);
        idx = idx + 1;
        
        %PMean, precise mean. It consider the mean of only the last sr
        %value of the signal
        featuresCell{idx} = mean(mydata((end-sr+1):end));
        idx = idx + 1;
        
        %%%% THIS CODE WAS TAKEN BY THE CITED PAPER.
        %%%% ANYWAY, IT CONSIDER THE CENTER OF THE SIGNAL.
        %%%% THE MOST IMPORTANT INFORMATION, IN MY STUDY
        %%%% ARE AT THE END OF THE SIGNAL..
%        
%         %Modfied means ABS. At the value will be added the weight 1 if the
%         %elements i (0,..,i,...,N) is 0.25N <= i <= 0.75N, 
%         %0.5 otherwise
%         myabsData1 = absData;
%         under = round(0.25 * length(mydata));
%         top = round(0.75 * length(mydata));
%         myabsData1(1:under) = myabsData1(1:under) .* 0.5;  
%         myabsData1(top:end) = myabsData1(top:end) .* 0.5;
%         featuresCell{idx} = sum(myabsData1) / length(mydata);
%         idx = idx + 1;
%
%         %Modfied means ABS. At the value will be added the weight 1 if the
%         %elements i (0,..,i,...,N) is 0.25N <= i <= 0.75N, 
%         %4i/N for i<0.25N, and 4(i-N)/N i>0.75N
%         myabsData2 = absData;
%         under = round(0.25 * length(mydata));
%         top = round(0.75 * length(mydata));
%         myabsData2(1:under) = myabsData2(1:under) .* ...
%            (4*(1:under)/length(mydata))';  
%         myabsData2(top:end) = myabsData2(top:end) .* ...
%            (4*((top:length(mydata))-length(mydata))/length(mydata))';
%         featuresCell{idx} = sum(myabsData2) / length(mydata);
%         idx = idx + 1;  

        %%%%END
      
        %The following code define as most important information the end of
        %the signal, similar to the guidelines described by the paper
        
        %Modfied means ABS. At the value will be added the weight 1 if the
        %elements i (0,..,i,...,N) is 0.75N >= i, 
        %0.5 otherwise
        myabsData1 = absData;
        top = round(0.75 * length(mydata));
        myabsData1(1:top) = myabsData1(1:top) .* 0.5;
        featuresCell{idx} = sum(myabsData1) / length(mydata);
        idx = idx + 1;
        
        %Modfied means ABS. At the value will be added the weight 1 if the
        %elements i (0,..,i,...,N) is 0.75N >= i, 
        %and an incremental weight from 0 to 1 otherwhise.
        myabsData2 = absData;
        approx = ((1:top) ./ top)';
        myabsData2(1:top) = myabsData2(1:top) .* approx;
        featuresCell{idx} = sum(myabsData2) / length(mydata);
        idx = idx + 1;
        
        %Mean Absolute Value Slope (MAVS). The difference between mean
        %absolute values of adjacent segment k and k-1
        featuresCell{idx} = mean(absData) -mean(absBefore);
        idx = idx + 1;
        
        %Root Mean Square (RMS). It is modeled as amplitude modulated
        %Gaussian random process whose RMS is related to the constant force
        %and non fatiguing contraction
        featuresCell{idx} = rms(mydata);
        idx = idx + 1;
        
        %Variance (VAR)
        featuresCell{idx} = var(mydata);
        idx = idx + 1;
        
        %Waveform Length (WL) is the cumulative length of the waveform over
        %the segment. (summarize in one parameter amplitude,freq, and
        %duration)
        featuresCell{idx} = sum(diff(mydata));
        idx = idx + 1;
        
        %Set threshold 
        [zctr, ssctr] = setThreshold(signaltype);
        
        %Zero Crossing (ZC) is the number of times that the waveform
        %crosses zero (number of time that it change the sign).
        %Here was added a threshold for avoid noise (multiplecrossing).
        zc = 0;
        zcarray = [];
        for i = 1:length(mydata)-1
           if(...
               (...
                   (mydata(i) > 0 && mydata(i+1) < 0) || ...
                   (mydata(i) < 0 && mydata(i+1 > 0 ))...
               ) && ...
                   (abs(mydata(i) - mydata(i+1)) >= zctr)...
              )
              zc = zc + 1; 
              zcarray = [zcarray i];
           end
        end
        featuresCell{idx} = zc;
        idx = idx + 1;
        
        %Slope Sign Changes (SSC) is the number of times the slope
        %changes sign. 
        sc = 0;
        scarray = [];
        for i = 2:length(mydata)-1
           if(...
               (...
                   (mydata(i) > mydata(i-1) && mydata(i+1) < mydata(i)) || ...
                   (mydata(i) < mydata(i-1) && mydata(i+1) > mydata(i))...
               ) && (...
                      (abs(mydata(i) - mydata(i+1)) >= ssctr) ||...
                      (abs(mydata(i) - mydata(i-1)) >= ssctr) ...
                     )...
              )
              sc = sc + 1; 
              scarray = [scarray i];
           end
        end
        featuresCell{idx} = sc;
        idx = idx + 1;
        
        %Willson Amplitude (WAMP) calculates the number of times that the
        %absolute value of the difference between signal amplitude of two
        %consecutive samples exceeds a predetermined threshold
        tr = 0.5;
        w = zeros(1,length(mydata)-1);
        for i = 1:(length(mydata)-1)
            if(abs(mydata(i)-mydata(i+1)) > tr)
                w(i) = 1;
            end
        end
        featuresCell{idx} = sum(w);
        idx = idx + 1;
        
        %Simple Square Integral (SSI) defines the energy of the signal
        ssi = abs(mydata.^2);
        featuresCell{idx} = sum(ssi);
        idx = idx + 1;
        
        
        %Frequency Median (FMD) splits the power spectrum density into two
        %equal parts
        [y,f] = mypsd(mydata,sr);
        featuresCell{idx} = (sum(y)/2);
        idx = idx + 1;
        
        %Frequency Mean (FMN) (see Oskoei, Mohammadreza Asghari, and 
        %Huosheng Hu. "GA-based feature subset selection for myoelectric
        %classification.") 
        dem = sum(y);
        num = sum(f(:).*y(:));
        featuresCell{idx} = num/dem;
        idx = idx + 1;
        
        %Modified Frequency Median (MFMD) is the frequency at which the
        %spectrum is divided into two regions with equal amplitude
        %----------------
        %SEE Phinyomark, Angkoon, Chusak Limsakul, and Pornchai
        %Phukpattaranont. "A novel feature extraction for robust EMG pattern recognition"
        %----------------
        L = length(mydata);
        NFFT = 2^nextpow2(L);
        mfmd = (1/2) * sum(abs(fft(mydata,NFFT)./L));
        featuresCell{idx} = mfmd;
        idx = idx + 1;
        
        
        %Modified Frequency Mean (MFMN) was proposed for Phinyomark (see
        %MFMD) an it is the average of the frequency.
        myf = sr/2*linspace(0,1,NFFT);
        Fy = abs(fft(mydata,NFFT)/L);
        featuresCell{idx} = sum(Fy.*myf') / sum(Fy);
        idx = idx + 1;
        
        %Frequency Ratio (FR) distinguish the difference between
        %contraction and relaxation of a muscle in frequency domain, by
        %appling fft 
        highfreq = bandpower(mydata,sr,hfreq);
        lowfreq = bandpower(mydata,sr,lfreq);
        featuresCell{idx} = abs(lowfreq)/abs(highfreq);
        %idx = idx + 1;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     IF IT IS AN EMG      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(exist('signaltype') > 0 && contains(lower(signaltype),'emg') > 0)
            emgCell = cell(1,1);
            emgNames = {'bho'};
            
            
            %featuresCell = [featuresCell,emgCell];
            %featuresNames = [featuresNames,emgNames];
        end
        %%%%%%%%%%%
        % END EMG %
        %%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     IF IT IS AN GSR      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(exist('signaltype') > 0 && contains(lower(signaltype),'gsr') > 0)
            gsrNames = {'MIN','MAX','NPeaks','PeaksMean'};
            gsrCell = cell(1,length(gsrNames));
            
            %Get the max and min value
            gsrCell{1} = min(mydata);
            gsrCell{2} = max(mydata);
            
            %Get the number of peaks
            [A,B] = findpeaks(medfilt1(mydata,sr),'MinPeakDistance',sr);
            gsrCell{3} = length(B);
            
            %Get the mean amp of the peaks
            gsrCell{4} = mean(mydata(B));
            
            
            featuresCell = [featuresCell,gsrCell];
            featuresNames = [featuresNames,gsrNames];
        end
        %%%%%%%%%%%
        % END EMG %
        %%%%%%%%%%%
        
        featuresTable =  cell2table(featuresCell, ...
            'VariableNames', featuresNames);
        

end


%TODO Look for features of GSR



function [zctr, ssctr] = setThreshold(signaltype)
    if(exist('signaltype') > 0 && contains(lower(signaltype),'emg') > 0)
        zctr = 0.5;
        ssctr = 0.5;

    elseif (exist('signaltype') > 0 && contains(lower(signaltype),'resp') > 0)
        zctr = 0.05;
        ssctr = 0.005;
    elseif (exist('signaltype') > 0 && contains(lower(signaltype),'gsr') > 0)
        zctr = 0.001;
        ssctr = 0.001;
    else
        zctr = 0.5;
        ssctr = 0.5;
    end
end        
        

function [bands,step,mymod,lfreq,hfreq] = getBandParameters(sr,signaltype)

if(exist('signaltype') > 0 && contains(lower(signaltype),'emg') > 0)
    bands = round(sr/2);
    step = 5;
    mymod = mod(bands,step);
    lfreq = [20,100];
    hfreq = [100,bands]; 
elseif (exist('signaltype') > 0 && contains(lower(signaltype),'resp') > 0)
    bands = 1;
    step = 0.05;
    mymod = 0;
    lfreq = [0.01, 0.2];
    hfreq = [0.2, 1];
elseif (exist('signaltype') > 0 && contains(lower(signaltype),'gsr') > 0)
    bands = 2;
    step = 0.05;
    mymod = 0;
    lfreq = [0.01, 0.2];
    hfreq = [0.2, 2];
else
    bands = round(sr/2);
    step = 5;
    mymod = mod(bands,step);
        lfreq = [20,100];
    hfreq = [100,bands];
end

end
