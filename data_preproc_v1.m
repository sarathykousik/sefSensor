clc;    close all;    clear all;

% Get data blanked
data_folder                 = ...
    'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO' ;
stim_chan                   = 'STI001';
pre_stim_time               = -2;
post_stim_time              = 7;
save_folder                 = ...
    'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO_blanked';
% blank_ft2fif(data_folder, stim_chan, pre_stim_time, post_stim_time, save_folder);

%% Declarations

blanked_folder              = save_folder;                                  %'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO\blanked';
cd(blanked_folder)
filename                    = []; 
filename                    = dir('*.fif');
loop                        = 2;

% Import data, Define trials, Cut trials
data_import                 = [];
cfg                         = [];
cfg.dataset                 = filename(loop).name;
data_import                 = ft_preprocessing(cfg);

% Filtering, Baseline correction and extraction of MEG channels only
preproc_data_MEG            = [];
cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfreq                  =  90;
cfg.hpfilter                = 'yes';
cfg.hpfreq                  =  1;
cfg.dftfilter               = 'yes';
cfg.channel                 =  'MEG*'; 
% cfg.detrend               =  'yes'; 
preproc_data_MEG            = ft_preprocessing(cfg, data_import);

data_trial = [];
cfg = [];
cfg.trialdef.prestim        = 0.05;
cfg.trialdef.poststim       = 0.35;
cfg.trialdef.eventtype      = 'STI001'; 
cfg.trialdef.eventvalue     = 5;
cfg.dataset                 = filename(loop).name;
data_trial                  = ft_definetrial(cfg);

data_epoched                = [];
data_epoched                = ft_redefinetrial(data_trial,preproc_data_MEG);

%%
% Artefact rejection
% cfg                         = [];
% cfg.artfctdef.feedback      = 'yes';
% clean_data                  = ft_rejectartifact(cfg, data_epoched);

% ICA

%%

% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.channel                 = {'MEG***3','MEG***2'};
% 
% for loop = 1:length(data_epoched.trial)-60
%     
%     loop
%     cfg.trials               = [loop:loop+50];
%     mov_avg(loop)            = ft_timelockanalysis(cfg, data_epoched);
%     
% end

%%
% Timelock-Averaging

cfg                         = [];
cfg.removemean              = 'yes';
% cfg.channel                 = {'MEG***3','MEG***2'};
% cfg.trials                  = [6:length(data_epoched.trial)];
timelocked_data             = ft_timelockanalysis(cfg,data_epoched );

cfg = [];
avgComb = ft_combineplanar(cfg,timelocked_data);

% Plotting Averaged ERP
cfg                         = [];
cfg.showlabels              = 'no'; 
% cfg.fontsize                = 6; 
% cfg.showlabels = 'yes';
% cfg.channel                 = {'MEG***3','MEG***2'};
cfg.layout                  = 'neuromag306all.lay';
cfg.ylim                    = [-5e-12 5e-12];
ft_multiplotER(cfg, avgCmb); 

%%


cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:90  % analysis 2 to 30 Hz in steps of 2 Hz 
% cfg.tapsmofreq   = 2;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.002   % length of time window = 0.5 sec
cfg.toi          = -0.05:0.002:0.3;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
TFRhann = ft_freqanalysis(cfg, timelocked_data);

%%
cfg = [];
cfg.baseline     = [-0.05 -0.005]; 
cfg.baselinetype = 'absolute'; 
cfg.clim         = [min(min(min(TFRhann.powspctrm)))*.85 max(max(max(TFRhann.powspctrm)))*1.25];	        
% cfg.showlabels   = 'yes';	
cfg.layout                  = 'neuromag306all.lay';
figure 
ft_multiplotTFR(cfg, TFRhann);

cfg                         = [];
cfg.showlabels              = 'no'; 
% cfg.fontsize                = 6; 
% cfg.showlabels = 'yes';
% cfg.channel                 = {'MEG***3','MEG***2'};
cfg.layout                  = 'neuromag306all.lay';
% cfg.ylim                    = [-5e-12 5e-12];
ft_multiplotER(cfg, TFRhann); 
