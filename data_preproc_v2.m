% clc;    close all;    clear all;

ft_defaults;
% Get data blanked
proc.data_folder                 = ...
    'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO' ;
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -2;
proc.post_stim_time              = 7;
proc.save_folder                 = ...
    'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO';

%% Declarations

proc.blanked_folder         = proc.save_folder;                                  %'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO\blanked';
cd(proc.blanked_folder)
filename                    = []; 
filename                    = dir('*.fif');
loop                        = 2;

% Import data, Define trials, Cut trials
proc.data_import            = [];
cfg                         = [];
cfg.dataset                 = filename(loop).name;
proc.data_import            = ft_preprocessing(cfg);

proc.data_trial = [];
cfg = [];
cfg.trialdef.prestim        = 0.15;
cfg.trialdef.poststim       = 0.35;
cfg.trialdef.eventtype      = proc.stim_chan ; 
cfg.trialdef.eventvalue     = 5
cfg.dataset                 = filename(loop).name;
proc.data_trial             = ft_definetrial(cfg);

proc.blank_seq          = ones(size(proc.data_import.trial{1},2),1);
for ind = proc.pre_stim_time:proc.post_stim_time
    proc.blank_seq(proc.data_trial.trl(:,1)+ind) = 0;
end

proc.MEG_channel = find(strncmp(proc.data_import.label, 'MEG',3)==1);                    % Pick only MEG channels
proc.data_import.trial{1}(proc.MEG_channel,:) = ...
        bsxfun(@times, proc.data_import.trial{1}(proc.MEG_channel,:), proc.blank_seq');  % Multiply blanking sequence with MEG channels

% Filtering, Baseline correction and extraction of MEG channels only
proc.preproc_data_MEG       = [];
cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfreq                  =  90;
cfg.hpfilter                = 'yes';
cfg.hpfreq                  =  1;
cfg.dftfilter               = 'yes';
cfg.channel                 = 'MEG'; 
cfg.detrend                 = 'yes';
% cfg.hilbert                 = 'angle';
proc.preproc_data_MEG       = ft_preprocessing(cfg, proc.data_import);

% proc.data_trial = [];
% cfg = [];
% cfg.trialdef.prestim        = 0.1;
% cfg.trialdef.poststim       = 0.35;
% cfg.trialdef.eventtype      = 'STI001'; 
% cfg.trialdef.eventvalue     = 5;
% cfg.dataset                 = filename(loop).name;
% proc.data_trial             = ft_definetrial(cfg);

proc.data_epoched           = [];
proc.data_epoched           = ft_redefinetrial(proc.data_trial,proc.preproc_data_MEG);

%% Timelock-Averaging

cfg                         = [];
cfg.removemean              = 'yes';
cfg.trials                  = 5:length(proc.data_epoched.trial)-5;   
cfg.keeptrials              = 'yes';
tlk_2                       = ft_timelockanalysis(cfg, proc.data_epoched);

% cfg                         = [];
% cfg.parameter               = 'trial';
% cfg.baseline                = [-0.090 -0.005];
% proc.timelocked_data_bsl    = ft_timelockbaseline(cfg, proc.timelocked_data);

% 
% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.trials                  = 10:length(proc.data_epoched.trial)-10;   
% cfg.channel                 =  'MEG***1'; 
% proc.timelocked_data_1      = ft_timelockanalysis(cfg, proc.data_epoched);
% 
% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.trials                  = 10:length(proc.data_epoched.trial)-10;   
% cfg.channel                 =  'MEG***2'; 
% proc.timelocked_data_2      = ft_timelockanalysis(cfg, proc.data_epoched);
% 
% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.trials                  = 10:length(proc.data_epoched.trial)-10;   
% cfg.channel                 =  'MEG***3'; 
% proc.timelocked_data_3      = ft_timelockanalysis(cfg, proc.data_epoched);

cfg                         = [];
cfg.parameter               = 'avg';    
cmb_tlkdata                 = ft_combineplanar(cfg, tlk_2);

cfg                         = [];
cfg.showlabels              = 'no';
cfg.channels                = 'MEG***1';
cfg.parameter               = 'avg';    
cfg.layout                  = 'neuromag306all.lay';
cfg.xlim                    = [0 0.1]; 
cfg.ylim                    = 'maxmin';
ft_multiplotER(cfg, cmb_tlkdata); 

%% Moving average

cfg                         = [];
% cfg.demean                  = 'yes';
% cfg.baseline                = [-0.05 -0.005];
cmb_tlkdata                 = ft_combineplanar(cfg, proc.timelocked_data);

cfg                         = [];
cfg.showlabels              = 'no';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.ylim                    = 'maxmin'; %[-30e-13 100e-13];
% cfg.xlim                    = [-0.005 0.06];
ft_multiplotER(cfg, cmb_tlkdata); 

%% Moving average

cfg                         = [];
cfg.removemean              = 'yes';


tlck = [];
trl_loop = 1;

for loop = 10:15:length(proc.data_epoched.trial)-35
    disp(['**********  ', num2str(loop)]);
    cfg.trials                 = loop:loop+30;
    tlck.mov_ave(trl_loop)     = ft_timelockanalysis(cfg, proc.data_epoched);
    trl_loop = trl_loop+1;
end

% cfg                         = [];
% cfg.showlabels              = 'no';
% cfg.layout                  = 'neuromag306all.lay';
% cfg.channel                 =  'MEG***1'; 
% 
% cfg.ylim                    = [-3e-13 3e-13];
% % cfg.ylim                    = [-80e-13 80e-13];
% % cfg.ylim                    = [-80e-13 80e-13];
% cfg.xlim                    = [-0.005 0.06];
% ft_multiplotER(cfg, tlck.mov_ave(1)); 

%% Extract channels of interest

MEG_grad = {'MEG0413', 'MEG0443', 'MEG1813'};
MEG_mag  = {'MEG0221', 'MEG0231', 'MEG0411', 'MEG0441'};

for loop = 1:length(tlck.mov_ave)

    maps.avg_mag(loop,:) = ...
        mean(tlck.mov_ave(1,loop).avg(match_str(tlck.mov_ave(1).label, MEG_mag),:),1);
    
    maps.avg_grad(loop,:) = ...
        mean(tlck.mov_ave(1,loop).avg(match_str(tlck.mov_ave(1).label, MEG_grad),:),1);
    
end

figure, imagesc(tlck.mov_ave(1).time, 1:length(tlck.mov_ave) ,maps.avg_grad)
xlabel('Time in seconds'), ylabel('Epoch windows')
axis xy, title('Average of grad data')

figure, imagesc(tlck.mov_ave(1).time, 1:length(tlck.mov_ave) ,maps.avg_mag)
xlabel('Time in seconds'), ylabel('Epoch windows')
axis xy, title('Average of mag data')

% Peak amplitude, latency  

N20_time = find(tlck.mov_ave(1).time>0.021 & tlck.mov_ave(1).time<0.030);
N30_time = find(tlck.mov_ave(1).time>0.028 & tlck.mov_ave(1).time<0.036);

for loop = 1:size(maps.avg_grad,1)

    [N20_amp(loop) N20_lat(loop)] = findpeaks(maps.avg_grad(loop,N20_time));
    [N30_amp(loop) N30_lat(loop)] = findpeaks(-maps.avg_grad(loop,N30_time));
    
end

N20_lat = N20_lat + min(N20_time);
N30_lat = N30_lat + min(N30_time);

figure,
subplot 221, plot(tlck.mov_ave(1).time(N20_lat), '-o'),axis([1 30 0.02 0.03])
title('N20 Latency'),  ylabel('Latency in seconds')
subplot 222, plot(tlck.mov_ave(1).time(N30_lat), '-o'),axis([1 30 0.028 0.036])
title('N30 Latency'), 
subplot 223, plot(N20_amp, '-o')
title('N20 Amplitude'),  xlabel('Epochs'), ylabel('Amplitude a.u.')
subplot 224, plot(N30_amp, '-o')
title('N30 Amplitude'), xlabel('Epochs')


%%
TFRhann = [];
TFRhann_bsl = [];
TFRhann_abs = [];
cfg              = [];
cfg.paramter     = 'trial';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.pad          = 0.5;
cfg.taper        = 'hanning';
cfg.foi          = 2.5:5:90;                          
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.020;   
cfg.toi          = -0.1:0.020:0.35;                 
TFRhann          = ft_freqanalysis(cfg, proc.data_epoched);

% TFRhann.powspctrm(isnan(TFRhann.powspctrm)) = 0;
cfg.baseline     = [-0.090 -0.01];
cfg.baselinetype = 'absolute';
cfg.parameter    = 'powspctrm';
TFRhann_bsl = ft_freqbaseline(cfg, TFRhann);
% TFRhann_bsl.powspctrm(isnan(TFRhann_bsl.powspctrm)) = 0;

cfg                         = [];
TFRhann_abs                 = ft_combineplanar(cfg, TFRhann_bsl);
TFRhann_abs.powspctrm(isnan(TFRhann_abs.powspctrm)) = 0;

cfg                         = [];
cfg.zlim                    = 'maxmin';
cfg.xlim                    = [0.01 0.34]; 
cfg.ylim                    = [1 90];    
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306all.lay';
ft_multiplotTFR(cfg, TFRhann_abs);

%%

cfg              = [];
cfg.output       = 'pow';
cfg.paramter     = 'trial';
cfg.pad          = 1;
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = [1:12,14:2:30,35:5:100];                          
cfg.t_ftimwin    = 2./cfg.foi;  % 7 cycles per time window
cfg.toi          = -0.1:0.01:0.35;                 
TFRhann7         = ft_freqanalysis(cfg,proc.data_epoched);

cfg.baseline     = [-0.090 -0.001];
cfg.baselinetype = 'absolute';
cfg.parameter    = 'powspctrm';
TFRhann_bsl_7    = ft_freqbaseline(cfg, TFRhann7);

cfg              = [];
TFRhann_abs_7    = ft_combineplanar(cfg, TFRhann_bsl_7);

cfg                         = [];
cfg.colorbar                = 'yes';
cfg.maskstyle               = 'saturation';	
cfg.layout                  = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, TFRhann_abs_7);

%% Morlet wavelet

cfg = [];
cfg.channel    = 'MEG';	                
cfg.method     = 'wavelet';                
cfg.width      = 2; 
cfg.gwidth     = 1;
cfg.output     = 'pow';	
cfg.foi        = [1:12,14:2:30,35:5:100];
cfg.toi        = -0.15:0.015:0.35;		              
TFRwave        = ft_freqanalysis(cfg, proc.data_epoched);

cfg = [];
cfg.baseline     = [-0.090 -0.01];
cfg.baselinetype = 'absolute';
cfg.parameter    = 'powspctrm';
TFRwave_bsl = ft_freqbaseline(cfg, TFRwave);

cfg                         = [];
TFRwave_cmb                 = ft_combineplanar(cfg, TFRwave_bsl);

cfg                         = [];
cfg.colorbar                = 'yes';
cfg.maskstyle               = 'saturation';	
cfg.layout                  = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, TFRwave_cmb);


%% ERD/ERS

% Alpha
cfg                         = [];
cfg.dftfilter               = 'yes';
cfg.channel                 =  'MEG'; 
cfg.detrend                 = 'yes';
cfg.bpfilter                = 'yes';
cfg.bpfreq                  = [8 12];
cfg.detrend                 = 'yes';
proc.alpha.MEG              = ft_preprocessing(cfg, proc.data_epoched);
% proc.alpha.MEG.trial{1}     = proc.alpha.MEG.trial{1}.^2;
proc.alpha.MEG_epoched      = ft_redefinetrial(proc.data_trial,proc.alpha.MEG);

cfg                         = [];
cfg.removemean              = 'yes';
cfg.keeptrials              = 'yes';
proc.alpha.tl               = ft_timelockanalysis(cfg, proc.alpha.MEG_epoched);

cfg                         = [];
cfg.parameter               = 'trial';
cfg.baseline                = [-0.090 -0.005];
proc.alpha.tl_bsl           = ft_timelockbaseline(cfg, proc.alpha.tl);
proc.alpha.MEG              = [];
proc.alpha.MEG_epoched      = [];


% Beta
cfg                         = [];
cfg.dftfilter               = 'yes';
cfg.channel                 =  'MEG'; 
cfg.detrend                 = 'yes';
cfg.bpfilter                = 'yes';
cfg.bpfreq                  = [13 25];
cfg.detrend                 = 'yes';
proc.beta.MEG               = ft_preprocessing(cfg, proc.data_import);
% proc.beta.MEG.trial{1}      = proc.beta.MEG.trial{1}.^2;
proc.beta.MEG_epoched       = ft_redefinetrial(proc.data_trial,proc.beta.MEG);

cfg                         = [];
cfg.removemean              = 'yes';
cfg.keeptrials              = 'yes';
proc.beta.tl                = ft_timelockanalysis(cfg, proc.beta.MEG_epoched);

cfg                         = [];
cfg.parameter               = 'trial';
cfg.baseline                = [-0.090 -0.005];
proc.beta.tl_bsl            = ft_timelockbaseline(cfg, proc.beta.tl);
proc.beta.MEG               = [];
proc.beta.MEG_epoched       = [];

% Low gamma
cfg                         = [];
cfg.dftfilter               = 'yes';
cfg.channel                 =  'MEG'; 
cfg.detrend                 = 'yes';
cfg.bpfilter                = 'yes';
cfg.bpfreq                  = [26 40];
cfg.detrend                 = 'yes';
proc.LGamma.MEG             = ft_preprocessing(cfg, proc.data_import);
% proc.LGamma.MEG.trial{1}      = proc.LGamma.MEG.trial{1}.^2;
proc.LGamma.MEG_epoched     = ft_redefinetrial(proc.data_trial,proc.LGamma.MEG);

cfg                         = [];
cfg.removemean              = 'yes';
cfg.keeptrials              = 'yes';
proc.LGamma.tl              = ft_timelockanalysis(cfg, proc.LGamma.MEG_epoched);

cfg                         = [];
cfg.parameter               = 'trial';
cfg.baseline                = [-0.090 -0.005];
proc.LGamma.tl_bsl          = ft_timelockbaseline(cfg, proc.LGamma.tl);
proc.LGamma.MEG             = [];
proc.LGamma.MEG_epoched     = [];

% High gamma
cfg                         = [];
cfg.dftfilter               = 'yes';
cfg.channel                 =  'MEG'; 
cfg.detrend                 = 'yes';
cfg.bpfilter                = 'yes';
cfg.bpfreq                  = [60 90];
cfg.detrend                 = 'yes';
proc.HGamma.MEG              = ft_preprocessing(cfg, proc.data_import);
% proc.HGamma.MEG.trial{1}      = proc.HGamma.MEG.trial{1}.^2;
proc.HGamma.MEG_epoched      = ft_redefinetrial(proc.data_trial,proc.HGamma.MEG);

cfg                         = [];
cfg.removemean              = 'yes';
cfg.keeptrials              = 'yes';
proc.HGamma.tl              = ft_timelockanalysis(cfg, proc.HGamma.MEG_epoched);

cfg                         = [];
cfg.parameter               = 'trial';
cfg.baseline                = [-0.090 -0.005];
proc.HGamma.tl_bsl          = ft_timelockbaseline(cfg, proc.HGamma.tl);
proc.HGamma.MEG             = [];
proc.HGamma.MEG_epoched     = [];

%%

% cfg                         = [];
% plot_1                      = ft_combineplanar(cfg, ERS_Alpha_wave_bsl);
% 
% cfg                         = [];
% cfg.showlabels              = 'no';
% cfg.parameter               = 'avg';
% cfg.layout                  = 'neuromag306cmb.lay';
% cfg.ylim                    = 'maxmin';
% cfg.xlim                    = [-0.09 0.4];
% ft_multiplotER(cfg, plot_1); 

%%

cfg              = [];
cfg.paramter     = 'trial';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.pad          = 1;
cfg.taper        = 'hanning';
cfg.foi          = 60:1:90;                          
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.020;   
cfg.toi          = -0.1:0.020:0.35;                 
ERS_Alpha_wave   = ft_freqanalysis(cfg, proc.HGamma.MEG_epoched);

cfg = [];
cfg.baseline     = [-0.090 -0.01];
cfg.baselinetype = 'absolute';
cfg.parameter    = 'powspctrm';
ERS_Alpha_wave_bsl = ft_freqbaseline(cfg, ERS_Alpha_wave);

cfg                         = [];
TFRwave_cmb                 = ft_combineplanar(cfg, ERS_Alpha_wave_bsl);

cfg                         = [];
cfg.colorbar                = 'yes';
cfg.maskstyle               = 'saturation';	
cfg.layout                  = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, TFRwave_cmb);

TFRwave_cmb.powspctrm = squeeze(mean(TFRwave_cmb.powspctrm,2));

cfg                         = [];
cfg.showlabels              = 'no';
% cfg.parameter               = 'avg';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.ylim                    = 'maxmin';
cfg.xlim                    = [-0.09 0.4];
ft_multiplotER(cfg, TFRwave_cmb); 
title('High gamma')

%%

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
% cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, proc.data_epoched);

%% Statistics
cfg =[];
% tlk_1_avg = ft_timelockanalysis(cfg, tlk_1)
freq_cmb  = ft_combineplanar(cfg, freq);

cfg                         = [];
cfg.showlabels              = 'no';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.ylim                    = 'maxmin'; %[-30e-13 100e-13];
% cfg.xlim                    = [-0.005 0.06];
ft_multiplotER(cfg, freq_cmb); 


cfg                         = [];
tlk_1_cmb                   = ft_combineplanar(cfg, tlk_1);

cfg                         = [];
tlk_2_cmb                   = ft_combineplanar(cfg, tlk_2);

%%
cfg = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.feedback    = 'yes'; % show a neighbour plot 
cfg.template    = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\neighbours\neuromag306cmb_neighb.mat';
neighbours      = ft_prepare_neighbours(cfg, tlk_1_cmb); % define neighbouring channels


%%
cfg = [];
cfg.channel     = 'MEG';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [0.015 0.060];
cfg.avgovertime = 'yes';
cfg.parameter   = 'trial';
cfg.tail = 0;
cfg.method      = 'montecarlo';
cfg.statistic   = 'indepsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
 
design = zeros(1,size(tlk_1_cmb.trial,1) + size(tlk_2_cmb.trial,1));
design(1,1:size(tlk_1_cmb.trial,1)) = 1;
design(1,(size(tlk_1_cmb.trial,1)+1):(size(tlk_1_cmb.trial,1) + size(tlk_2_cmb.trial,1)))= 2;

cfg.design = design;             % design matrix
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 stat = [];
stat = ft_timelockstatistics(cfg,tlk_1_cmb,tlk_2_cmb)
 
% make a plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'neuromag306cmb.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, tlk_1_avg_cmb)
title('Nonparametric: significant with cluster multiple comparison correction')


cfg = [];
cfg.highlightsymbolseries = ['*','*','.','.','.'];
cfg.layout    = 'neuromag306cmb.lay';
cfg.contournum = 0;
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.parameter='stat';
cfg.zlim = [-5 5];
ft_clusterplot(cfg,stat);
