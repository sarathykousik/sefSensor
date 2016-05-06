%%
% Sim 10, 15, 45, 70, 90 Hz =>  append   => TFR with same settings as dpss
loop =1;

for freqloop = [91:3:120]
    disp(['################## ', num2str(loop)])
    cfg = [];
    cfg.method  = 'superimposed';
    cfg.fsample = 1000;
    cfg.numtrl  = 50;
    cfg.trllen  = 0.2;

    cfg.s1.freq  = [freqloop]; 
    cfg.s1.ampl = 1;
    cfg.s2.phase = 0;

    % cfg.s2.freq = [60];
    % cfg.s2.ampl = 1;
    % cfg.s2.phase = 0;

    % cfg.s3.freq = [80];
    % cfg.s3.ampl = 1;
    % cfg.s3.phase = 0;

    data_1 = ft_freqsimulation(cfg);

    cfg = [];
    cfg.method  = 'superimposed';
    cfg.fsample = 1000;
    cfg.numtrl  = 50;
    cfg.trllen  = 0.2;

    cfg.s1.freq  = [6]; 
    cfg.s1.ampl = 1;
    cfg.s1.phase = 0;

    cfg.s2.freq = [12];
    cfg.s2.ampl = 1;
    cfg.s2.phase = 0;

    cfg.s2.freq = [25];
    cfg.s2.ampl = 1;
    cfg.s2.phase = 0;

    data_2 = ft_freqsimulation(cfg);

    % data = ft_appenddata([], data_1, data_2);
    data = [];
    for tr  = 1:numel(data_2.trial)
        data.trial{tr} = [data_2.trial{tr} data_2.trial{tr} ...
            data_1.trial{tr} data_2.trial{tr} data_2.trial{tr}];
        data.time{tr} = 1e-3:1e-3:1.1;
        data.fsample = 1000;
        data.label = data_1.label;
        data.cfg = data_1.cfg;
        data.dimord = 'rpt_time';
    end

%
% figure(1)
% plot(data.time{1}, data.trial{1}(1,:))

   cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'mix';
    cfg.method       = 'mtmconvol';
%     cfg.pad          = 2;
%     cfg.foi          = 21;                          
    cfg.foi          = 65;                          
    cfg.t_ftimwin    = 5./cfg.foi;   
    cfg.toi          = 0:0.010:1.200;           
    cfg.tapsmofrq    = 25;
%     cfg.tapsmofrq    = 8;
    gammaFreq        = ft_freqanalysis(cfg, data);
    gammaDesc{loop} = ft_freqdescriptives([],gammaFreq);
    loop=loop+1;
end

freq_str = num2str(91:3:120);
leg_freq = strsplit(freq_str);

    cfg                         = [];
    cfg.zlim                    = 'maxmin';
    cfg.channel = 'mix'
    cfg.ylim                    = [0.0 0.2]; 
%     cfg.xlim                    = [0.12 0.5]
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.axes                    = 'no';
    cfg.graphcolor              = 'brgkycmbrgkycm';
    ft_singleplotER(cfg, gammaDesc{1:10});
    legend(leg_freq{1:10})
    
title('Simulation(Noise:6,12,25) - dpss,Fc=65Hz BW=25Hz')
%     figure,plot(gammaFreq.freq,squeeze(mean(gammaFreq.powspctrm(:,1,:))))
    