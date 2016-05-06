
clear all; clc; close all;
 
%%
cd('J:\MEG_Research\SEF\SEF-TSSS\Patients\cpyTrans')
filenames           = dir('*.mat');
 
for loop = 1:length(filenames)
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['######## ', b])
    load(filenames(loop).name)
    
    
    cfg                 = [];
    cfg.channel         = 'MEGGRAD';
    visClean            = ft_rejectvisual(cfg, proc.data_epoched);

%     save (['J:\MEG_Research\SEF\SEFVisClean\',b, '-visClean'], 'visClean', '-v7.3');
%     saveData = [];
end


%%

tlck_cmb_new = ft_combineplanar([],ft_timelockanalysis([],proc.data_epoched))
figure
cfg                     = [];
% cfg.parameter           = 'powspctrm';
% cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
cfg.hlim = [0 0.2];
% cfg.channel = {'all','-MEG1442+1443'};
ft_multiplotER(cfg, tlck_cmb,tlck_cmb_new)


%% TFR

cfg                     = [];
cfg.output              = 'pow';
cfg.channel             = 'MEGGRAD';
cfg.method              = 'wavelet';
cfg.pad                 = 2;
cfg.foi                 = 10:2:100;                          
cfg.width               = 5; 
cfg.toi                 = -0.2:0.010:0.2;            
TFRwave_vis             = ft_freqanalysis(cfg, visClean);

cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relchange';
cfg.parameter               = 'powspctrm';
TFRwave_bsl               = ft_freqbaseline(cfg,TFRwave);


cfg                     = [];
cfg.parameter           = 'powspctrm';
% cfg.channel = {'all','-MEG0432+0433', '-MEG0442+0443', '-MEG1812+1813'}
% cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
% cfg.zlim = [0.9 3];
cfg.xlim = [0 0.2];
% cfg.channel = {'all','-MEG1442+1443'};
% cfg.graphcolor = 'rbgkcm';
% cfg.zlim=[2 3];
ft_multiplotTFR(cfg, ft_combineplanar([],TFRwave_bsl))


%% Gamma
cfg              = [];
cfg.paramter     = 'trial';
cfg.taper        = 'dpss';
cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel      = 'MEGGRAD';
cfg.method       = 'mtmconvol';
cfg.pad          = 2;
cfg.foi          = 75;                          
cfg.t_ftimwin    = 7./cfg.foi;   
cfg.toi          = -0.200:0.020:0.200;           
cfg.tapsmofrq    = 20;
gammaFreq        = ft_freqanalysis(cfg, saveData.dataArtRej);

%
cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relchange';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);

cfg             = [];
cfg.jackknife   = 'yes';
gammaDesc       = ft_freqdescriptives(cfg, gammaFreq_bsl);

cfg                         = [];
gammaCmb                    = ft_combineplanar(cfg, gammaFreq_bsl);

cfg                     = [];
cfg.parameter           = 'powspctrm';
% cfg.channel = {'all','-MEG0432+0433', '-MEG0442+0443', '-MEG1812+1813'}
% cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
% cfg.zlim = [0.9 3];
% cfg.xlim = [0 0.1];
% cfg.channel = {'all','-MEG1442+1443'};
% cfg.graphcolor = 'rbgkcm';
% cfg.zlim=[2 3];
ft_multiplotER(cfg, gammaCmb)

%%

cfg_l             = [];
cfg_l.method      = 'mtmfft';
cfg_l.foi         = [1:200];
cfg_l.taper       = 'hanning';
cfg_l.channel     = 'MEGGRAD'
pow  = ft_freqanalysis(cfg_l, saveData.dataArtRej);

% figure,
% semilogy(pow.freq, mean(pow.powspctrm)), hold on

cfg                     = [];
cfg.parameter           = 'powspctrm';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
ft_multiplotER(cfg, ft_combineplanar([], pow))





