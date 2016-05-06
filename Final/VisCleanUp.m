
clear all; clc; close all;
 
%%
filenames           = dir('*.mat');
% 

for loop = 1:length(filenames)
    load(filenames(loop).name)
    cfg                 = [];
    visClean            = ft_rejectvisual(cfg,saveData.dataArtRej);

    [a b c]     =   fileparts(filenames(loop).name);
    save (['J:\MEG_Research\SEF\SEFVisClean','\',b, '-visClean'], 'visClean', '-v7.3');
    saveData = [];
end

%%
for loop = 1:length(filenames)
    disp('********************************')
    load(filenames(loop).name)
    [a b c]     =   fileparts(filenames(loop).name);

    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'MEG';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 1;
    cfg.taper        = 'hanning';
    cfg.foi          = [40:90];                          
    cfg.t_ftimwin    = 4./cfg.foi;   
    cfg.toi          = -0.5:0.01:0.5;           
    cfg.tapsmofrq    = cfg.foi*0.1;
    TFR              = ft_freqanalysis(cfg, visClean);

    save (['J:\MEG_Research\SEF\SEFcleanTFR','\',b, '-TFR'], 'TFR', '-v7.3');
    TFR      = [];
    saveData = [];
    disp('###############################')
end 

%% Gamma mtp: Fc = 65 Hz, BW: 50 Hz

cfg              = [];
cfg.paramter     = 'trial';
cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
% cfg.pad          = 1;
% cfg.taper        = 'hanning';
cfg.foi          = 65;                          
cfg.t_ftimwin    = 4./cfg.foi;   
cfg.toi          = -0.5:0.01:0.5;           
cfg.tapsmofrq    = 25;
% cfg.keeptapers   = 'yes';
gammaFreq        = ft_freqanalysis(cfg, visClean);

cfg = [];
cfg.baseline                = [-0.090 -0.01];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg, gammaFreq);

cfg                         = [];
gammaFreq_cmb                = ft_combineplanar(cfg, gammaFreq_bsl);

cfg                         = [];
% cfg.zlim                    = 'maxmin';
cfg.xlim                    = [0.02 0.35]; 
cfg.ylim                    = 'maxmin';    
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.axes                    = 'no';
figure,ft_multiplotER(cfg, gammaFreq_cmb);


























