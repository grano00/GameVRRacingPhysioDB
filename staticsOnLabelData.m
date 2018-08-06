vrLabel = tableOfFiles.VRInGame;
novrLabel = tableOfFiles.NOVRInGame;

VRAURG1 = [];
VRAURG2 = [];
VRVALG1 = [];
VRVALG2 = [];

VRAUR = [];
VRVAL = [];

NOVRAURG1 = [];
NOVRAURG2 = [];
NOVRVALG1 = [];
NOVRVALG2 = [];

NOVRAUR = [];
NOVRVAL = [];


AURG1 = [];
AURG2 = [];
VALG1 = [];
VALG2 = [];


    cs = @(x)((x-50)./50);
    
PAUR = [];

PVAL = [];

for i = 1:length(vrLabel)
    load(vrLabel{i});
    disp(['work on ' vrLabel{i}]);
    labelData{1} = cs(labelData{1});
    labelData{2} = cs(labelData{2});
    vrAurG1 = labelData{1}(:,2);
    vrAurG2 = labelData{2}(:,2);
    vrValG1 = labelData{1}(:,3);
    vrValG2 = labelData{2}(:,3);
    vrAur = [vrAurG1; vrAurG2];
    vrVal = [vrValG1; vrValG2];
    
    load(novrLabel{i});
    labelData{1} = cs(labelData{1});
    labelData{2} = cs(labelData{2});
    novrAurG1 = labelData{1}(:,2);
    novrAurG2 = labelData{2}(:,2);
    novrValG1 = labelData{1}(:,3);
    novrValG2 = labelData{2}(:,3);
    novrAur = [novrAurG1; novrAurG2];
    novrVal = [novrValG1; novrValG2];
    
    AurG1 = [vrAurG1; novrAurG1];
    AurG2 = [vrAurG2; novrAurG2];
    
    ValG1 = [vrValG1; novrValG1];
    ValG2 = [vrValG2; novrValG2];
    
    VRAURG1 = [VRAURG1;vrAurG1];
    VRAURG2 = [VRAURG2;vrAurG2];
    VRVALG1 = [VRVALG1;vrValG1];
    VRVALG2 = [VRVALG2;vrValG2];
    
    VRAUR = [VRAUR; vrAur];
    VRVAL = [VRVAL; vrVal];
    
    NOVRAURG1 = [NOVRAURG1;novrAurG1];
    NOVRAURG2 = [NOVRAURG2;novrAurG2];
    NOVRVALG1 = [NOVRVALG1;novrValG1];
    NOVRVALG2 = [NOVRVALG2;novrValG2];
    
    NOVRAUR = [NOVRAUR; novrAur];
    NOVRVAL = [NOVRVAL; novrVal];

    AURG1 = [AURG1;AurG1];
    AURG2 = [AURG2;AurG2];
    
    VALG1 = [VALG1;ValG1];
    VALG2 = [VALG2;ValG2];
    
    PAUR = [PAUR; mean([AURG1;AURG2])];
    PVAL = [PVAL; mean([VALG1;VALG2])];
    
end

clearvars -except AU* VA* VR* NOVR* PAUR PVAL tableOfFiles

AUR = [AURG1;AURG2];
VAL = [VALG1;VALG2];
    
DATA = cell2table(cell(0,18), 'VariableNames', {'AUR','VAL','AURG1', 'AURG2', 'VALG1', 'VALG2', 'VRAUR', 'NOVRAUR', 'VRVAL', ...
    'NOVRVAL', 'VRAURG1','VRAURG2', 'NOVRAURG1','NOVRAURG2','VRVALG1','VRVALG2','NOVRVALG1','NOVRVALG2'});

dev = @(x)(sum((x - mean(x)).^2));
sdr = @(x)(std(x)/abs(mean(x))); %standard deviation relative


dataNames = {'MEAN','STD','VAR','DEV','MIN','MAX','SDR','MAD'};
DATA = [DATA;{mean(AUR),mean(VAL), mean(AURG1), mean(AURG2), mean(VALG1), mean(VALG2), mean(VRAUR), mean(NOVRAUR), mean(VRVAL), mean(NOVRVAL), mean(VRAURG1), mean(VRAURG2), mean(NOVRAURG1), mean(NOVRAURG2), mean(VRVALG1), mean(VRVALG2), mean(NOVRVALG1), mean(NOVRVALG2)}];
DATA = [DATA;{std(AUR),std(VAL),std(AURG1), std(AURG2), std(VALG1), std(VALG2), std(VRAUR), std(NOVRAUR), std(VRVAL), std(NOVRVAL), std(VRAURG1), std(VRAURG2), std(NOVRAURG1), std(NOVRAURG2), std(VRVALG1), std(VRVALG2), std(NOVRVALG1), std(NOVRVALG2)}];
DATA = [DATA;{var(AUR),var(VAL),var(AURG1), var(AURG2), var(VALG1), var(VALG2), var(VRAUR), var(NOVRAUR), var(VRVAL), var(NOVRVAL), var(VRAURG1), var(VRAURG2), var(NOVRAURG1), var(NOVRAURG2), var(VRVALG1), var(VRVALG2), var(NOVRVALG1), var(NOVRVALG2)}];
DATA = [DATA;{dev(AUR),dev(VAL),dev(AURG1), dev(AURG2), dev(VALG1), dev(VALG2), dev(VRAUR), dev(NOVRAUR), dev(VRVAL), dev(NOVRVAL), dev(VRAURG1), dev(VRAURG2), dev(NOVRAURG1), dev(NOVRAURG2), dev(VRVALG1), dev(VRVALG2), dev(NOVRVALG1), dev(NOVRVALG2)}];
DATA = [DATA;{min(AUR),min(VAL),min(AURG1), min(AURG2), min(VALG1), min(VALG2), min(VRAUR), min(NOVRAUR), min(VRVAL), min(NOVRVAL), min(VRAURG1), min(VRAURG2), min(NOVRAURG1), min(NOVRAURG2), min(VRVALG1), min(VRVALG2), min(NOVRVALG1), min(NOVRVALG2)}];
DATA = [DATA;{max(AUR),max(VAL),max(AURG1), max(AURG2), max(VALG1), max(VALG2), max(VRAUR), max(NOVRAUR), max(VRVAL), max(NOVRVAL), max(VRAURG1), max(VRAURG2), max(NOVRAURG1), max(NOVRAURG2), max(VRVALG1), max(VRVALG2), max(NOVRVALG1), max(NOVRVALG2)}];
DATA = [DATA;{sdr(AUR),sdr(VAL),sdr(AURG1), sdr(AURG2), sdr(VALG1), sdr(VALG2), sdr(VRAUR), sdr(NOVRAUR), sdr(VRVAL), sdr(NOVRVAL), sdr(VRAURG1), sdr(VRAURG2), sdr(NOVRAURG1), sdr(NOVRAURG2), sdr(VRVALG1), sdr(VRVALG2), sdr(NOVRVALG1), sdr(NOVRVALG2)}];
DATA = [DATA;{mad(AUR),mad(VAL),mad(AURG1), mad(AURG2), mad(VALG1), mad(VALG2), mad(VRAUR), mad(NOVRAUR), mad(VRVAL), mad(NOVRVAL), mad(VRAURG1), mad(VRAURG2), mad(NOVRAURG1), mad(NOVRAURG2), mad(VRVALG1), mad(VRVALG2), mad(NOVRVALG1), mad(NOVRVALG2)}];
DATA.Properties.RowNames = dataNames;
DATA.Indexes = dataNames';
DATA = [DATA.Indexes, DATA(:,1:end-1)];
writetable(DATA,'labelAnalysis.xls','Sheet',1);

PDATA = array2table([PAUR,PVAL], 'VariableNames', {'AUR','VAL'});

writetable(PDATA,'labelAnalysis.xls','Sheet',2);