clc;    close all;    clear all;

ft_defaults;
% Get data blanked
proc.data_folder                 = ...
    'J:\data-tsss-sef\SEF-TSSS' ;
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              = 10;
proc.save_folder                 = ...
    'J:\data-tsss-sef\SEF-TSSS';
% blank_ft2fif(data_folder, stim_chan, pre_stim_time, post_stim_time, save_folder);

%% Declarations

proc.blanked_folder         = proc.save_folder;                                  %'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO\blanked';
cd(proc.blanked_folder)
cd('010_XKG')

filename                    = []; 
filename                    = dir('*.fif');
loop                        = 1;

% Import data, Define trials, Cut trials
proc.data_import            = [];
cfg                         = [];
cfg.dataset                 = filename(loop).name;
proc.data_import            = ft_preprocessing(cfg);

proc.data_trial = [];
cfg = [];
cfg.trialdef.prestim        = 0.1;
cfg.trialdef.poststim       = 0.35;
cfg.trialdef.eventtype      = 'STI001'; 
cfg.trialdef.eventvalue     = 5;
cfg.dataset                 = filename(loop).name;
proc.data_trial             = ft_definetrial(cfg);

proc.blank_seq          = ones(size(proc.data_import.trial{1},2),1);
for ind = proc.pre_stim_time:proc.post_stim_time
    proc.blank_seq(proc.data_trial.trl(:,1)+100+ind) = 0;
end

proc.MEG_channel = find(strncmp(proc.data_import.label, 'MEG',3)==1);                    % Pick only MEG channels
proc.data_import.trial{1}(proc.MEG_channel,:) = ...
        bsxfun(@times, proc.data_import.trial{1}(proc.MEG_channel,:), proc.blank_seq');  % Multiply blanking sequence with MEG channels

% Filtering, Baseline correction and extraction of MEG channels only
proc.preproc_data_MEG       = [];
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

%% Jump
cfg                    = [];
cfg.trl                = proc.data_trial.trl;
cfg.continuous         = 'yes';
 
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel    = 'MEG';
cfg.artfctdef.zvalue.cutoff     = 30;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
 
% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';
 
% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';
[cfg, artifact_jump] = ft_artifact_zvalue(cfg, proc.data_import);

%% EOG
cfg                              = [];
cfg.trl                          = proc.data_trial.trl;
cfg.continuous                   = 'yes'; 

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = 'EOG*';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.bpfreq      = [1 15];
cfg.artfctdef.zvalue.bpfiltord   = 4;
cfg.artfctdef.zvalue.hilbert     = 'yes';

% feedback
cfg.artfctdef.zvalue.interactive = 'yes';
[cfg,artifact_EOG] = ft_artifact_zvalue(cfg, proc.data_import);

%% EMG

% muscle
cfg            = [];
cfg.trl        = proc.data_trial.trl;
cfg.continuous = 'yes'; 

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = 'EMG*';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.fltpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter    = 'yes';
cfg.artfctdef.zvalue.bpfreq      = [110 140];
cfg.artfctdef.zvalue.bpfiltord   = 9;
cfg.artfctdef.zvalue.bpfilttype  = 'but';
cfg.artfctdef.zvalue.hilbert     = 'yes';
cfg.artfctdef.zvalue.boxcar      = 0.2;

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_muscle] = ft_artifact_zvalue(cfg, proc.data_import);
  
  
%% Artifact reject

cfg                             =   []; 
cfg.artfctdef.reject            = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
cfg.artfctdef.eog.artifact      = artifact_EOG; % 
cfg.artfctdef.jump.artifact     = artifact_jump;
cfg.artfctdef.muscle.artifact   = artifact_muscle;
cfg.artfctdef.crittoilim        = [0.015 0.15];
% cfg.artfctdef.minaccepttim      = 0.1;
data_no_artifacts               = ft_rejectartifact(cfg, proc.data_epoched);

%% 
cfg             = [];
cfg.channel     = 'MEG*';
data            = ft_preprocessing(cfg, proc.preproc_data_MEG);

cfg             = [];
cfg.resamplefs  = 150;
cfg.detrend     = 'no';
data            = ft_resampledata(cfg, data);

cfg             = [];
cfg.runica.pca  = 50;
cfg.method      = 'runica';
comp            = ft_componentanalysis(cfg, data);

cfg             = [];
cfg.component   = [1:50];       % specify the component(s) that should be plotted
cfg.layout      = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
cfg.channel     = {'MEG***1'};%, 'MEG***3'};
cfg.comment     = 'no';
ft_topoplotIC(cfg, comp)

cfg             = [];
cfg.channel     = [6,14]; % components to be plotted
cfg.viewmode    = 'component';

cfg.layout      = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
ft_databrowser(cfg, comp)

%% Timelock-Averaging

cfg                         = [];
cfg.removemean              = 'yes';
proc.timelocked_data        = ft_timelockanalysis(cfg, proc.data_epoched);

cfg                         = [];
cmb_tlkdata                 = ft_combineplanar(cfg, proc.timelocked_data);

cfg                         = [];
cfg.showlabels              = 'no';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.ylim                    = 'maxmin'; 
ft_multiplotER(cfg, cmb_tlkdata); 

cfg = [];
cfg.channel     = 'MEG' ; %{'MEG1432', 'MEG1433', 'MEG1332', 'MEG1333', 'MEG1322'};
figure,ft_databrowser(cfg, proc_noBlank.data_import)