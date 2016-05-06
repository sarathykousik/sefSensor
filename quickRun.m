% function output = fif_file_preproc(file_name,subj, condN ,proc)
clc;    clear all;     close all;
ft_defaults;
% addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')
no = 8

proc                             = [];
proc.data_folder                 = 'J:\MEGData';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              =  10;
% proc.save_folder                 = 'J:\MEGData';

cd(proc.data_folder)
filenames      = dir('*.fif');

%%
cfg                         = [];
cfg.dataset                 = filenames(no).name;
proc.data_import            = ft_preprocessing(cfg);

% Find epochs
cfg = [];
cfg.trialdef.prestim        = 0.5;
cfg.trialdef.poststim       = 0.5;
cfg.trialdef.eventtype      = proc.stim_chan ; 
cfg.trialdef.eventvalue     = 5;
cfg.dataset                 = filenames(no).name;%filenames(condLoop).name;
proc.data_trial             = ft_definetrial(cfg);

% Blanking
proc.blank_seq              = ones(size(proc.data_import.trial{1},2),1);
for ind = proc.pre_stim_time:proc.post_stim_time
    proc.blank_seq(proc.data_trial.trl(:,1)+(cfg.trialdef.prestim*1000)+ind) = 0;
end

proc.MEG_channel            = find(strncmp(proc.data_import.label, 'MEG',3)==1);                    % Pick only MEG channels
proc.data_import.trial{1}(proc.MEG_channel,:) = ...
        bsxfun(@times, proc.data_import.trial{1}(proc.MEG_channel,:), proc.blank_seq');  % Multiply blanking sequence with MEG channels

% Filtering, Baseline correction and extraction of MEG channels only
cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfreq                  =  100;
cfg.hpfilter                = 'yes';
cfg.hpfreq                  =  2;
cfg.dftfilter               = 'yes';
cfg.channel                 =  {'MEG'}; 
cfg.detrend                 = 'yes';
proc.preproc_data_MEG       = ft_preprocessing(cfg, proc.data_import);

proc.data_epoched           = ft_redefinetrial(proc.data_trial,proc.preproc_data_MEG);

%
cfgtlck = [];
% cfgtlck.channel = 'MEGGRAD';
cfgtlck.removemean = 'yes';
proc.tlck = ft_timelockanalysis(cfgtlck,proc.data_epoched);

proc.tlckCmb = ft_combineplanar([],proc.tlck);

% Use ft_math to diff datasets
%%
cfg                     = [];
cfg.parameter           = 'avg';
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.interplimits   = 'electrodes';
cfg.maskstyle           = 'saturation';	 
% Delete artefacted channels
cfg.channel = {'all', '-MEG0242+0243', '-MEG1612+1613','-MEG0232+0233' ,'-MEG0212+0213'};

figure(1);
cfg.xlim                = [0.02 0.03];
ft_topoplotER(cfg, proc.tlckCmb)


figure(2);
cfg.xlim                = [0.03 0.045];
ft_topoplotER(cfg, proc.tlckCmb)

figure(3);
cfg.xlim                = [0.05 0.065];
ft_topoplotER(cfg, proc.tlckCmb)

