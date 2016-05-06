%% Calculate Induced, evoked gamma
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc; clear all; close all;

%%
filePath = 'J:\MEG_Research\SEF\SEFVisClean';  
save_folder = 'J:\MEG_Research\SEF\SEFVisClean\SEFIndGamma';

cd(filePath)
filenames      = dir('*.mat');



for loop=1:length(filenames)
    [a b c]     =   fileparts(filenames(loop).name);    
    disp(['#########  ',b])
    load(filenames(loop).name)
   % Calc ERF
    cfg = []; cfg.channel = 'MEGGRAD'; %cfg.trials = 'all';
    tlck = ft_timelockanalysis(cfg, visClean)
   
   % Calc Ind data (visClean-ERF)
    for trialLoop = 1: length(visClean.trial)
        ind_dat{1,trialLoop} = ...
            visClean.trial{1,trialLoop} - tlck.avg;
    end

   
   % Calc Ind gamma
   indGamma = calcGamma(ind_dat);
   
   % Calc Evoked gamma
   evokedGamma = calcGamma(tlck);
   
   save ([save_folder, '/',b, '-indGamma'], 'indGamma');
   save ([save_folder, '/',b, '-evokedGamma'], 'evokedGamma');
end

% calc gamma
function gammaStruct = calcGamma(cleanData)
    
    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'MEGGRAD';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
    cfg.foi          = 75;                          
    cfg.t_ftimwin    = 5./cfg.foi;   
    cfg.toi          = -0.500:0.010:0.500;           
    cfg.tapsmofrq    = 20;
    gammaStruct        = ft_freqanalysis(cfg, cleanData);


return




