clear all; clc; close all;

filenames = dir('*.mat');
%%
% for i=1:length(filenames)
%     disp(i)
%     load(filenames(i).name);
   
    cfg = [];
    cfg.baseline                = [-0.1 -0.02];
    cfg.baselinetype            = 'relative';
    cfg.parameter               = 'powspctrm';
    TFR_bsl                     = ft_freqbaseline(cfg, TFRwave);

    cfg                         = [];
    TFR_bsl                     = ft_combineplanar(cfg, TFR_bsl);

    cfg                         = [];
%     cfg.channel                 = TFR_bsl.label{16};
    cfg.zlim                    = [0 3.5];
    cfg.xlim                    = [0.02 0.08]; 
    cfg.ylim                    = [60 100];    
    cfg.maskstyle               = 'saturation';	
%     cfg.colormap   = 'gray';
    cfg.shading                 = 'interp';
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    figure,ft_topoplotTFR(cfg, TFR_bsl);
    
    
    cfg                         = [];
    cfg.zlim                    = 'maxmin';
    cfg.xlim                    = [0 0.1]; 
    cfg.ylim                    = [10 100];    
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    figure,ft_multiplotTFR(cfg, TFR_bsl);
    
%     TFRsave = [];
    
% end

%%

cfg              = [];
cfg.paramter     = 'trial';
cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.pad          = 1;
cfg.taper        = 'hanning';
cfg.foi          = [1:2:30];                          
cfg.t_ftimwin    = 2./cfg.foi;   
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.01;
cfg.toi          = -0.5:0.01:0.5;                 
TFRhann_1_30          = ft_freqanalysis(cfg, saveData.dataArtRej);

%%
cfg              = [];
cfg.paramter     = 'trial';
% cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.pad          = 1;
cfg.taper        = 'hanning';
cfg.foi          = [1:100];                          
cfg.t_ftimwin    = 4./cfg.foi;   
cfg.toi          = -0.5:0.005:0.5;           
cfg.tapsmofrq    = cfg.foi*0.1;
TFRhann_1_100   = ft_freqanalysis(cfg, dummy);

cfg = [];
cfg.baseline                = [-0.090 -0.01];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';
TFRhann_bsl = ft_freqbaseline(cfg, TFRhann_1_100);

cfg                         = [];
TFRhann_bsl                 = ft_combineplanar(cfg, TFRhann_bsl);

cfg                         = [];
cfg.zlim                    = [-15e-24 15e-24];
cfg.xlim                    = [-0.1 0.4]; 
cfg.ylim                    = [1 100];    
cfg.maskstyle               = 'saturation';	
cfg.parameter               = 'powspctrm';
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, TFRhann_bsl);

%% Evoked power
cfg = [];
avgSEF = ft_timelockanalysis(cfg, saveData.dataArtRej);

cfg              = [];
cfg.paramter     = 'trial';
% cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.pad          = 1;
cfg.taper        = 'hanning';
cfg.foi          = [1:100];                          
cfg.t_ftimwin    = 4./cfg.foi;   
cfg.toi          = -0.5:0.005:0.5;           
cfg.tapsmofrq    = cfg.foi*0.1;
evokedPow   = ft_freqanalysis(cfg, avgSEF);

cfg = [];
cfg.baseline                = [-0.090 -0.01];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';
TFRhann_bsl = ft_freqbaseline(cfg, TFRhann_31_100);

cfg                         = [];
TFRhann_bsl                 = ft_combineplanar(cfg, TFRhann_bsl);

cfg                         = [];
cfg.zlim                    = [-5e-24 5e-24];
cfg.xlim                    = [-0.1 0.4]; 
cfg.ylim                    = [1 100];    
cfg.maskstyle               = 'saturation';	
cfg.parameter               = 'powspctrm';
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, TFRhann_bsl);



