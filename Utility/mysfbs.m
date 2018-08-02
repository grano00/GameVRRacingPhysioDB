function [ indexes_f, hist_f, err_f, selected_f, cvinfo_f] = mysffs(X, Y, ...
    fun, log, niteration,cv,varNames, randomSeed)
%SFFS Definition of a personal way of SFFS. 
%   Parameters are:

%X: Data to elaborate
%Y: Target variable
%fun: AI Method. It must return a cell array with in the 1th col the
%discrimation info
%log: is an integer, it can be 0 for don't display anything, 1 for log, 2
%for display dots..
%niteration: Number of iteration of the algo. With a greater ammout of
%iteration the algorithm should perform better
%cv: a cvpartition
%varNames: the variable names that return fun
%randomSeed: a starting seed (rng)


%Definition of the new variable names.
%On they will be added the information of indexes used for CV
varNames = [varNames 'TestIdx'];

%Set new rng and generate a new seed
rseed = randi([1 500]);
selected = 0;
if(exist('randomSeed') > 0 && floor(randomSeed) == randomSeed)
    rseed = randomSeed;
end
rng(rseed);

%Set empty final variables. The err_f is equal to inf in order to replace
%it at the first iteration
err_f = inf;
cvinfo_f = {};

%Loop for each iteration
for interation = 1:niteration
    %Set new Seed
    rseed = randi([1 500]);
    rng(rseed);
    
    %Define the parallel methods
    opts = statset('UseParallel',true);
    
    displ(log,['Start interation ' num2str(interation) ' with seed ' num2str(rseed)]);
    
    %Setup starting variables
    [R, C] = size(X);
    indexes = zeros(1,C);
    hist = indexes;
    err = 0;
    k = 1;
    maxErr = inf;
    cvinfo = {};
    cvinfoID = 1;
    
    %Loop until find a minimun local
    myexit = true;
    while myexit
        
        tempErr = zeros(1,C);
        ti = indexes;
        
        %Forward Selection
        displ(log,['Selected ' num2str(sum(ti)) '/' num2str(length(ti))]);
        tempCVInfo = cell(1,C);
        parfor i = 1:C
            %If the feature are not already selected
            if ti(i) == 0
                %add it
                lti = logical(ti);
                lti(i) = true;
                rng(rseed);
                %Cross validation
                redo = true;
                while redo
                    try
                        
                        mycvinfo = crossval(fun,X(:,lti),Y,...
                            'partition',cv,'Options',opts);
                        redo = false;
                    catch ME
                        disp('FIND EXCEPTION:');
                        disp(ME);
                        disp('REDO THE CV');
                    end
                end
                
                %Save the data
                cvidx = getIdx(cv);
                tempCVInfo{i} = {cv2tab(mycvinfo,cvidx,varNames)};
                myerr = cell2deepmat(mycvinfo(:,1));
                tempErr(i) = mean((myerr ./ (cv.TestSize')));
                
            else
                tempErr(i) = maxErr;
            end
        end
        
        %Select the feature that better improve prediction
        [A,I] = min(tempErr);
        [maxErr,~] = max(tempErr);
        displ(log,['MIN = ' num2str(A) ' in POS ' num2str(I)]);
        
        err(k) = A;
        indexes(I) = 1;
        
        
        cvinfo=[cvinfo; {{cvinfoID},tempCVInfo(I)}];
        cvinfoID = cvinfoID + 1;
        
        
        %Check if it was already appeared or the misclassification is perfect
        %or it is a lot of time (10) that was not find any better error
        [~,index] = ismember(hist,indexes,'rows');
        
        hist(k,:) = indexes;
        [~,s] = min(err);
        if(sum(index) > 0 || A == 0 || s < (length(err) - 10))
            displ(log,'find in hist, exit');
            [nerr,selected] = min(err);
            indexes = hist(selected,:);
            mycvinfoindex = cell2deepmat(cvinfo(:,1));
            tCVInfo= cvinfo(mycvinfoindex==selected,:);
            if(nerr ~= err(k))
                err(k+1) = nerr;
                hist(k+1,:) = indexes;
                cvinfo=[cvinfo; tCVInfo];
                selected = k+1;
            end
            %If it is true, exit from the loop
            disp('d');
            fprintf('\n\n');
            break;
        end
        
        %Check if cvinfo is greater than 15 and delete the previous
        %infos (for save the memory).
        mybuffer = 12;
        if(length(cvinfo) > mybuffer)
            cvinfo(1:(end-mybuffer),:) = [];
        end
        
        [HR HC] = size(hist);
        %It check if are setted at least 3 variables for start the backward
        %part
        if(HR > 2)
            tempErr = zeros(1,C);
            ti = indexes;
            %Backward
            
            displ(log,['Selected ' num2str(sum(ti)) '/' num2str(length(ti))]);
            parfor i = 1:C
                if ti(i) == 1
                    lti = logical(ti);
                    lti(i) = false;
                    rng(rseed);
                    redo = true;
                    while redo
                        try
                            mycvinfo = crossval(fun,X(:,lti),Y,...
                                'partition',cv,'Options',opts);
                            redo = false;
                        catch ME
                            disp('FIND EXCEPTION: ');
                            disp(ME);
                            disp('REDO THE CV')
                        end
                    end
                    cvidx = getIdx(cv);
                    tempCVInfo{i} = {cv2tab(mycvinfo,cvidx,varNames)};
                    myerr = cell2deepmat(mycvinfo(:,1));
                    tempErr(i) = mean((myerr ./ (cv.TestSize')));
                    
                else
                    tempErr(i) = maxErr;
                end
            end
            
            [A,I] = min(tempErr);
            displ(log,['BACK MIN ' num2str(A) ' is in pos ' num2str(I)]);
            
            [maxErr,~] = max(tempErr);
            if(A < err(k))
                displ(log,'less so delete');
                k = k + 1;
                indexes(I) = 0;
                hist(k,:) = indexes;
                err(k) = A;
                cvinfo=[cvinfo;  {{cvinfoID},tempCVInfo(I)}];
            else
                
            end
        end
        
        k = k + 1;
    end
    
    %Save info of the better iteration that must return
    if(err(end) < err_f(end))
        displ(log,['iteration ' num2str(interation) ...
            'has a lower error: SELECTED']);
        indexes_f= indexes;
        hist_f = hist;
        err_f = err;
        selected_f = selected;
        cvinfo_f = cvinfo;
    end
end
end

%Get the list of indexes starting from one col of cv cell
function cellcv = getIdx(cv)
n = cv.NumTestSets;
cellcv = cell(n,1);
for i = 1:n
    myidx = test(cv,i);
    cellcv{i} = myidx;
end
end

%Convert cross validation into a table
function [T] = cv2tab(mycell,cv,myVarNames)
[R,C] = size(mycell);
if(R==1 && C ==1), mycell = mycell{1}; end

T = cell2table([mycell cv],'VariableNames',myVarNames);
end

%Show the info. -Debug Version-
function [] = displ(log,text)

switch log
    case 1
        fprintf([text ' .-. \n']);
    case 2
        fprintf('.')
end
end