clc;    clear all;     close all;

ft_defaults;

%%
proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFVisClean\';
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFGammaInduced';
mkdir(proc.save_folder)

load('J:\MEG_Research\SEF\SEFepresults\allTlck.mat');
allTlck=allTlck{1:10};
%
cd(proc.data_folder)
filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
    
%     for cond = 1:col   %col
    [a b c]     =   fileparts(char(file_sub(subj, cond)));
    disp(['   **** visClean   ',b,'     ****    '])
    disp(['   #### visClean   ',tlck{subj,cond}.name,'     ####'])
    load(char(file_sub(subj, cond)));
    visCleanGrad = [];  ind_dat = [];

    cfg= []; cfg.channel = 'MEGGRAD';
    visCleanGrad = ft_selectdata(cfg,visClean);
    
    tlckGrad     = ft_timelockanalysis([],visCleanGrad);

    for trialLoop = 1: length(visCleanGrad.trial)
        ind_dat{1,trialLoop} = ...
            visCleanGrad.trial{1,trialLoop} - tlck{subj,cond}.avg;
    end

    visCleanInd = visCleanGrad;
    visCleanInd.trial = ind_dat;
    visCleanInd.filename =char(file_sub(subj, cond));
        
%         save (['J:\MEG_Research\SEF\SEFCleanInd\',b, '-visCleanInd'], 'visCleanInd', '-v7.3');
        
end
% end