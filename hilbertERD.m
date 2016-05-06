cfg = []; 
cfg.hilbert = 'abs';
cfg.bpfilter  = 'yes';
cfg.bpfreq   = [40 90];    
cfg.bpfiltdir = 'twopass';
cfg.bpfilttype = 'but';
cfg.bpinstabilityfix = 'split';
cfg.bpfiltord     = 50;
% cfg.trials = 1:100;
hilbData = ft_preprocessing(cfg, visClean);

%
cfg                     = [];
hilbComb                = ft_combineplanar(cfg, hilbData);

cfg = [];
tlckHilb = ft_timelockanalysis(cfg, hilbComb);

cfg = [];
cfg.baseline = [-0.09 -0.01];
cfg.baselinetype = 'absolute';
% cfg.parameter = 'trial';
hilbBase = ft_timelockbaseline(cfg, tlckHilb);

cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.xlim = [-0.1 0.1];
% cfg.ylim = [0.2e-11 1.19e-11];
ft_multiplotER(cfg,hilbBase)


% cfg = [];
% % cfg.parameter = 'avg';
% cfg.layout    = 'neuromag306cmb.lay';
% cfg.xlim = [0.02 0.06];
% ft_topoplotER(cfg,hilbBase)