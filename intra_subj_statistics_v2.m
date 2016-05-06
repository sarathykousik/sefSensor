clc;    clear all;  close all;

%%
data_folder = 'J:\data-tsss-sef\SEF-TSSS';
_folder                 = 'J:\results_FT\avgs';
cd(save_folder)
filename                    = []; 
filename                    = dir('*.mat');

[num,txt,raw] = xlsread('J:\data-tsss-sef\SEF-TSSS\subjectDataList.xls');
num_cell = num2cell(num);

subjName = {raw{1,6:10}};
condName = {raw{4:10,1}};


%% Initial processing
cfg = [];
for subjLoop = 1:length(subjName)
    %Choose filenames
    [row, col]=find(strcmp(raw, subjName{subjLoop})==1);
    disp('==========================================')
    disp('==========================================')
    disp(['   *******      ',subjName{subjLoop}, '    *********'])
    disp('==========================================')
    disp('==========================================')
    fileList = {raw{4:10,col}};
    
        for condLoop  = 1:fileList
            disp('#######################################')
            disp(['   ****    ',num2str(condName{condLoop}),'     ****    '])
            disp('#######################################')
            % Load data
            workData{loop}     = load(filename(loop).name);
            workData{loop}.filename = filename(loop).name;

            % Combineplanar
            cmbData_tlck{loop}          = ...
                            ft_combineplanar(cfg, workData{loop}.saveData.tlck_avg);
            cmbData_tlck{loop}.filename = filename(loop).name;
            cmbData_TFR{loop}           = ...
                            ft_combineplanar(cfg, workData{loop}.saveData.TFR_hann);
            cmbData_TFR{loop}.filename  = filename(loop).name;


            % Re-arrange data

        end
end

%% Timelock
% Statistics
cfg = [];
cfg.channel     = 'MEG';
cfg.latency     = [0.3 0.7];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT'
cfg.alpha       = 0.05;
cfg.correctm    = 'yes';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
 
Nsub = 2;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,cmbData_tlck{1:2},cmbData_tlck{2:3});
% Plots
% make the plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'neuromag306cmb.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, cmbData_tlck{1})
title('Nonparametric: significant without multiple comparison correction')

%% TFR
% Statistics

% Plots
 
