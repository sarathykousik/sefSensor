%% TFRcalc
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\ArtRejtSSS\SEFdata';
proc.save_folder                 = 'J:\MEG_Research\ArtRejtSSS\SEFdata\results';
mkdir(proc.save_folder)

%%
cd(proc.data_folder)
filenames      = dir('*.mat');

for loop  = 2:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
   %
    cfg              = [];
%     cfg.paramter     = 'trial';
% cfg.trials = [1:10];
%     cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'MEGGRAD';
%     cfg.channelcmb = {visClean.label{31} 'MEG',  visClean.label{32} 'MEG'};
    cfg.method       = 'wavelet';
    cfg.pad          = 2;
    cfg.foi          = 10:2:100;                          
    cfg.width        = 5; 
    cfg.toi          = -0.12:0.020:0.2;            
    TFRwave           = ft_freqanalysis(cfg, saveData.dataArtRej);
%     save ([proc.save_folder,'\',b, '-TFRwave'], 'TFRwave', '-v7.3');

    cfg = [];
    cfg.baseline                = [-0.1 -0.02];
    cfg.baselinetype            = 'relative';
    cfg.parameter               = 'powspctrm';
    TFRwave_bsl                 = ft_freqbaseline(cfg, TFRwave);

    cfg                         = [];
    TFRwave_bsl_cmb             = ft_combineplanar(cfg, TFRwave_bsl);
    save ([proc.save_folder,'\',b, '-TFRwave_bsl'], 'TFRwave_bsl', '-v7.3');
       

end

%%
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
TFR_logpow    = ft_math(cfg, TFRwave);

%%

    cfg                         = [];
    cfg.zlim                    = 'maxmin';
    cfg.xlim                    = [0 0.1]; 
    cfg.ylim                    = [1 100];
    cfg.maskstyle               = 'saturation';	
        cfg.shading                 = 'interp';

    cfg.channel                 = {'all','-MEG1322+1323', '-MEG1512+1513',...
        '-MEG1522+1523', '-MEG1912+1913'};
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.shading                 = 'interp';
    figure,
    ft_multiplotTFR(cfg, ft_combineplanar([],TFRwave_bsl));


%%

chan_left =  {'MEG0112+0113', 'MEG0122+0123', 'MEG0132+0133',...
    'MEG0142+0143', 'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233',...
    'MEG0242+0243', 'MEG0312+0313', 'MEG0322+0323', 'MEG0332+0333', ...
    'MEG0342+0343', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433',...
    'MEG0442+0443', 'MEG0512+0513', 'MEG0522+0523', 'MEG0532+0533', ...
    'MEG0542+0543', 'MEG0612+0613', 'MEG0622+0623', 'MEG0632+0633', ...
    'MEG0642+0643', 'MEG0712+0713', 'MEG0742+0743', 'MEG0812+0813', ...
    'MEG0822+0823', 'MEG1012+1013', 'MEG1512+1513', 'MEG1522+1523',...
    'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623',...
    'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723',...
    'MEG1732+1733', 'MEG1742+1743', 'MEG1812+1813', 'MEG1822+1823', ...
    'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923',...
    'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2042+2043',...
    'MEG2112+2113', 'MEG2122+2123', 'MEG2142+2143'};
%%

chan_left = {'MEG0112+0113', 'MEG0122+0123', 'MEG0132+0133', 'MEG0142+0143',...
    'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243',...
    'MEG0322+0323', 'MEG0332+0333', 'MEG0342+0343', 'MEG0412+0413',...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0612+0613',...
    'MEG0622+0623', 'MEG0632+0633', 'MEG0642+0643', 'MEG0712+0713',...
    'MEG0742+0743', 'MEG1012+1013', 'MEG1512+1513', 'MEG1522+1523',...
    'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623',...
    'MEG1632+1633', 'MEG1642+1643', 'MEG1722+1723', 'MEG1812+1813',...
    'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913',...
    'MEG2012+2013', 'MEG2042+2043'};

cfg = [];
cfg.baseline                = [-0.15 -0.10];
cfg.baselinetype            = 'relchange';
cfg.parameter               = 'powspctrm';
TFRwave_bsl               = ft_freqbaseline(cfg, gammaFreq);

% cfg           = [];
% cfg.parameter = 'powspctrm';
% cfg.operation = 'log10';
% TFR_logpow    = ft_math(cfg, gammaFreq);

cfg = [];
cfg.latency = [0 0.2];
% cfg.channel = chan_left;
gamma_baseline = ft_selectdata(cfg, ft_combineplanar([],TFRwave_bsl));
gamma_baseline.powspctrm = zeros(size(gamma_baseline.powspctrm));

cfg = [];
cfg.latency = [0 0.2];
% cfg.channel = chan_left;
gamma_activation = ft_selectdata(cfg, ft_combineplanar([],TFRwave_bsl));
gamma_baseline.time = gamma_activation.time;

% load('j:\MEG_Research\SEF\neighbours.mat')

cfg = [];
% cfg.latency           = [0.02 0.08];
cfg.statistic           = 'ft_statfun_actvsblT';
% cfg.channel             = chan_left;
% cfg.method              = 'analytic';
cfg.correctm            = 'bonferonni';
cfg.method            = 'montecarlo';
% cfg.correctm          = 'cluster';
cfg.clusteralpha        = 0.01;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 1;
cfg.tail                = 2;
cfg.clustertail         = 2;
cfg.alpha               = 0.01;
cfg.numrandomization    = 500;
cfg_neighb.method       = 'distance';
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb, gamma_activation);

    ntrials = size(gamma_activation.powspctrm,1);
    design  = zeros(2,2*ntrials);
    design(1,1:ntrials) = 1;
    design(1,ntrials+1:2*ntrials) = 2;
    design(2,1:ntrials) = [1:ntrials];
    design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat_log     = ft_freqstatistics(cfg,  gamma_activation, gamma_baseline);

cfg                     = [];
cfg.parameter           = 'stat';
cfg.maskparameter       = 'mask';
cfg.maskstyle           = 'box';
cfg.interplimits        = 'electrodes';
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
figure, ft_multiplotER(cfg, stat_log);


%%

cfg = [];
cfg.baseline                = [-0.15 -0.10];
cfg.baselinetype            = 'relchange';
cfg.parameter               = 'powspctrm';
TFRwave_bsl               = ft_freqbaseline(cfg, gammaFreq);

% cfg           = [];
% cfg.parameter = 'powspctrm';
% cfg.operation = 'log10';
% TFR_logpow    = ft_math(cfg, gammaFreq);

cfg = [];
cfg.latency = [0 0.2];
% cfg.channel = chan_left;
gamma_activation = ft_selectdata(cfg, ft_combineplanar([],TFRwave_bsl));
% gamma_baseline.time = gamma_activation.time;

%%

cfg=[];
% cfg.latency=[0.02 0.2];
cfg.method='stats';
cfg.statistic='ttest_samples_vs_const';
% cfg.tail = 'right';
cfg.constantvalue = 0;
cfg.alpha = 0.01;
% cfg.method    = 'analytic';
% cfg.correctm  = 'bonferonni';

    ntrials = size(SEF_ITC_cond{1,1}.powspctrm,1);
    design  = zeros(1,1*ntrials);
    design(1,1:ntrials) = 1;
%     design(1,ntrials+1:2*ntrials) = 2;
%     design(2,1:ntrials) = [1:ntrials];
%     design(2,ntrials+1:2*ntrials) = [1:ntrials];
cfg.design  = design;
cfg.ivar     = 1;
cfg.uvar     = 2;

[stat] = ft_freqstatistics(cfg, SEF_ITC_cond{1,1});

% figure(1),imagesc(squeeze(stat.mask)), colorbar, title('mask')
% figure(2),imagesc(squeeze(stat.prob)), colorbar, title('prob')
cfg                     = [];
cfg.parameter           = 'stat';
cfg.maskparameter       = 'mask';
cfg.maskstyle           = 'box';
cfg.interplimits        = 'electrodes';
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
figure(3),ft_multiplotER(cfg, stat);

