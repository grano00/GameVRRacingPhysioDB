function [  ] = produceResults( tableOfFiles,mainDir )
%PRODUCERESULTS Summary of this function goes here
%   Detailed explanation goes here

%Create matrixes
names = {'NOVRGame1','NOVRGame2','NOVR',...
    'VRGame1','VRGame2','VR',...
    'Game1','Game2'};
vrFiles = tableOfFiles.VRbyUserResults;
novrFiles = tableOfFiles.NOVRbyUserResults;
gameFiles = tableOfFiles.GAMESbyUserResults;

A_rmse = zeros(length(vrFiles),length(names));
A_mse = A_rmse;
A_mae = A_rmse;
A_r = A_rmse;
A_r2 = A_rmse;

V_rmse = zeros(length(vrFiles),length(names));
V_mse = V_rmse;
V_mae = V_rmse;
V_r = V_rmse;
V_r2 = V_rmse;

for i = 1:length(novrFiles)
   load(novrFiles{i});
   for j = 1:3 
       mydata = Results{1,j};
       arousalData = mydata.arousalFin;
       valenceData = mydata.valenceFin;
       
       %Assing arousal data
       A_rmse(i,j) = arousalData.rmse;
       A_mse(i,j) = arousalData.mse;
       A_mae(i,j) = arousalData.mae;
       A_r(i,j) = arousalData.r;
       A_r2(i,j) = ((arousalData.r)^2);
       
       %Assing valence data
       V_rmse(i,j) = valenceData.rmse;
       V_mse(i,j) = valenceData.mse;
       V_mae(i,j) = valenceData.mae;
       V_r(i,j) = valenceData.r;
       V_r2(i,j) = ((valenceData.r)^2);
   end
end
for i = 1:length(vrFiles)
   load(vrFiles{i});
   for j = 4:6 
       mydata = Results{1,j};
       arousalData = mydata.arousalFin;
       valenceData = mydata.valenceFin;
       
       %Assing arousal data
       A_rmse(i,j) = arousalData.rmse;
       A_mse(i,j) = arousalData.mse;
       A_mae(i,j) = arousalData.mae;
       A_r(i,j) = arousalData.r;
       A_r2(i,j) = ((arousalData.r)^2);
       
       %Assing valence data
       V_rmse(i,j) = valenceData.rmse;
       V_mse(i,j) = valenceData.mse;
       V_mae(i,j) = valenceData.mae;
       V_r(i,j) = valenceData.r;
       V_r2(i,j) = ((valenceData.r)^2);
   end
end
for i = 1:length(gameFiles)
   load(gameFiles{i});
   for j = 7:8 
       mydata = Results{1,j};
       arousalData = mydata.arousalFin;
       valenceData = mydata.valenceFin;
       
       %Assing arousal data
       A_rmse(i,j) = arousalData.rmse;
       A_mse(i,j) = arousalData.mse;
       A_mae(i,j) = arousalData.mae;
       A_r(i,j) = arousalData.r;
       A_r2(i,j) = ((arousalData.r)^2);
       
       %Assing valence data
       V_rmse(i,j) = valenceData.rmse;
       V_mse(i,j) = valenceData.mse;
       V_mae(i,j) = valenceData.mae;
       V_r(i,j) = valenceData.r;
       V_r2(i,j) = ((valenceData.r)^2);
   end
end


arousal = struct;
arousal.rmse = array2table(A_rmse,'VariableNames',names);
arousal.mse = array2table(A_mse,'VariableNames',names);
arousal.mae = array2table(A_mae,'VariableNames',names);
arousal.r = array2table(A_r,'VariableNames',names);
arousal.r2 = array2table(A_r2,'VariableNames',names);

valence = struct;
valence.rmse = array2table(V_rmse,'VariableNames',names);
valence.mse = array2table(V_mse,'VariableNames',names);
valence.mae = array2table(V_mae,'VariableNames',names);
valence.r = array2table(V_r,'VariableNames',names);
valence.r2 = array2table(V_r2,'VariableNames',names);

npath = [mainDir '/IndividualResults/'];
mkdir(npath);
save([npath 'PerfIndexes.mat'],'arousal','valence');

%Write Excel
fname = [npath 'arousalPerformanceIndexes.xlsx'];
writetable(arousal.rmse,fname,'Sheet','rmse');
writetable(arousal.mse,fname,'Sheet','mse');
writetable(arousal.mae,fname,'Sheet','mae');
writetable(arousal.r,fname,'Sheet','r');
writetable(arousal.r2,fname,'Sheet','rsquare');

fname = [npath 'valencePerformanceIndexes.xlsx'];
writetable(valence.rmse,fname,'Sheet','rmse');
writetable(valence.mse,fname,'Sheet','mse');
writetable(valence.mae,fname,'Sheet','mae');
writetable(valence.r,fname,'Sheet','r');
writetable(valence.r2,fname,'Sheet','rsquare');

%Create boxplots
properties = GetProperties();
h = figure('units','points','outerposition',[0 0 2560 1080]);
if(properties.plotSilentMode)
  set(h, 'Visible', 'off');
end
bplotPath = [npath 'boxPlot/'];
mkdir(bplotPath);

%Arousal
%-------
%RMSE
boxplot(A_rmse,names);
title('Arousal RMSE BoxPlots');
xlabel('Experiments');
ylabel('RMSE');
print([bplotPath 'arousal-rmse-boxplot'],'-dpng');
savefig([bplotPath 'arousal-rmse-boxplot']);
%MSE
boxplot(A_mse,names);
title('Arousal MSE BoxPlots');
xlabel('Experiments');
ylabel('MSE');
print([bplotPath 'arousal-mse-boxplot'],'-dpng');
savefig([bplotPath 'arousal-mse-boxplot']);
%MAE
boxplot(A_mae,names);
title('Arousal MAE BoxPlots');
xlabel('Experiments');
ylabel('MAE');
print([bplotPath 'arousal-mae-boxplot'],'-dpng');
savefig([bplotPath 'arousal-mae-boxplot']);
%R
boxplot(A_r,names);
title('Arousal R BoxPlots');
xlabel('Experiments');
ylabel('R');
print([bplotPath 'arousal-r-boxplot'],'-dpng');
savefig([bplotPath 'arousal-r-boxplot']);
%R^2
boxplot(A_r2,names);
title('Arousal R-Squared BoxPlots');
xlabel('Experiments');
ylabel('R-Squared');
print([bplotPath 'arousal-r2-boxplot'],'-dpng');
savefig([bplotPath 'arousal-r2-boxplot']);


%Valence
%-------
%RMSE
boxplot(V_rmse,names);
title('Valence RMSE BoxPlots');
xlabel('Experiments');
ylabel('RMSE');
print([bplotPath 'valence-rmse-boxplot'],'-dpng');
savefig([bplotPath 'valence-rmse-boxplot']);
%MSE
boxplot(V_mse,names);
title('Valence MSE BoxPlots');
xlabel('Experiments');
ylabel('MSE');
print([bplotPath 'valence-mse-boxplot'],'-dpng');
savefig([bplotPath 'valence-mse-boxplot']);
%MAE
boxplot(V_mae,names);
title('Valence MAE BoxPlots');
xlabel('Experiments');
ylabel('MAE');
print([bplotPath 'valence-mae-boxplot'],'-dpng');
savefig([bplotPath 'valence-mae-boxplot']);
%R
boxplot(V_r,names);
title('Valence R BoxPlots');
xlabel('Experiments');
ylabel('R');
print([bplotPath 'valence-r-boxplot'],'-dpng');
savefig([bplotPath 'valence-r-boxplot']);
%R^2
boxplot(V_r2,names);
title('Valence R-Squared BoxPlots');
xlabel('Experiments');
ylabel('R-Squared');
print([bplotPath 'valence-r2-boxplot'],'-dpng');
savefig([bplotPath 'valence-r2-boxplot']);

close all
end

