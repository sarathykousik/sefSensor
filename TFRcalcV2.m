%% TFRcalc
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc;    clear all;     close all;

ft_defaults;
addpath('J:\MEG_Research\SEF\SEFcleanData')

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFcleanData';
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFcleanTFR';
mkdir(proc.save_folder)

%% Files

[temp.num,temp.txt,temp.raw] = xlsread('J:\MEG_Research\SEF\SEF-TSSS\subjectDataList.xls');
temp.num_cell = num2cell(temp.num);

subjName = {temp.raw{1,6:10}};
condName = {temp.raw{4:10,1}};

%% Import, epoch, filter, artefact data

for subjLoop = 1:length(subjName)

    %Choose filenames
    [row, col]=find(strcmp(temp.raw, subjName{subjLoop})==1);
    disp('==========================================')
    disp('==========================================')
    disp(['   *******      ',subjName{subjLoop}, '    *********'])
    disp('==========================================')
    disp('==========================================')
    mkdir([proc.save_folder,'\',subjName{subjLoop}]);

    fileList = {temp.raw{4:10,col}}; 
    cd([proc.data_folder,'\', subjName{subjLoop}])
    
    for condLoop  = 1:length(fileList)
        disp('#######################################')
        disp(['   ****    ',num2str(condName{condLoop}),'     ****    '])
        disp('#######################################')

        dataFile = [proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'-SaveData.mat'];
        cleanData = load(dataFile);
        [a b c]     =   fileparts(dataFile);


        cfg              = [];
        cfg.paramter     = 'trial';
        cfg.keeptrials   = 'yes';
        cfg.output       = 'pow';
        cfg.channel      = 'MEG';
        cfg.method       = 'mtmconvol';
        cfg.pad          = 2;
        cfg.taper        = 'hanning';
        cfg.foi          = 1:5:100;                          
        cfg.t_ftimwin    = 2./cfg.foi;   
        cfg.toi          = -0.5:0.01:0.5;                 
        TFRhann          = ft_freqanalysis(cfg, cleanData.saveData.dataArtRej);

        TFRsave.file     = dataFile;
        TFRsave.subject  = cleanData.saveData.subject;
        TFRsave.condName = cleanData.saveData.condName;
        TFRsave.TFR      = TFRhann;

        save ([proc.save_folder,'\', subjName{subjLoop},'\',b, '-TFR'], 'TFRsave', '-v7.3');
        disp('*********************************')
        disp(['Saved data to:   ', [proc.save_folder,'\', subjName{subjLoop},'\',b, '-TFR']]);
        disp('*********************************')

%        Check for files
%                 disp(subjName{subjLoop})
%                 disp(fileList{condLoop})
%                 disp(num2str(condName{condLoop}))
%                 exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'-SaveData.mat'],'file')
    end
end




