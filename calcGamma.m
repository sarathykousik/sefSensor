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