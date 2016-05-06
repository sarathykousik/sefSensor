%% Data preprocessing

cfg             = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      =  90;
cfg.hpfilter    = 'yes';
cfg.hpfreq      =  1;
cfg.dftfilter   = 'yes';
cfg.channel     =  'MEG*'; 
preproc_data    = ft_preprocessing(cfg, data);

%% Epoching

cfg3                         = [];
% cfg3.vartrllength            = 1;
trl = [data_trial.trl(:,1)'-50 ;  data_trial.trl(:,1)'+250; ...
                zeros(1,size(data_trial.trl(:,1)),1)];
cfg3.trl = trl';

data_epoched = [];
data_epoched = ft_redefinetrial(cfg3,preproc_data);

for loop = 1:length(data_epoched.time)
    data_epoched.time{loop} = data_epoched.time{loop}-.05;
end

%% Timelock

cfg = [];
avg_trials = ft_timelockanalysis(cfg, data_epoched);

%%
cfg = [];
cfg.showlabels = 'no'; 
cfg.fontsize = 6; 
% cfg.showlabels = 'yes';
cfg.channel = {'MEG***1'};
cfg.layout = 'neuromag306all.lay';
cfg.ylim = [-5e-13 5e-13];
ft_multiplotER(cfg, avg_trials); 



%% Head model

mri = ft_read_mri('C:\Program Files\MATLAB\R2012b\toolbox_add_on\bauer_m.mri');

cfg          = [];
% cfg.coordsys = 'fif'; % the MRI is expressed in the CTF coordinate system, see below
segmentedmri = ft_volumesegment(cfg, mri);

segmentedmri.transform = mri.transform;
segmentedmri.anatomy   = mri.anatomy;

% cfg.method = 'interactive';
% % cfg.parameter = 'anatomy';
% [realign_mri] = ft_volumerealign(cfg, mri)


cfg              = [];
cfg.funparameter = 'gray';
ft_sourceplot(cfg,segmentedmri);

cfg        = [];
cfg.method = 'singleshell';
hdm        = ft_prepare_headmodel(cfg, segmentedmri);
save prep_hdm.hdm hdm 
%%
templatedir  = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel';
template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm')); 

cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template.sourcemodel;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri;
sourcemodel        = ft_prepare_sourcemodel(cfg);

%%
cfg = [];
cfg.latency = [0.025 0.035];  % specify latency window around M50 peak
cfg.numdipoles = 1;
cfg.headmodel = hdm;
cfg.feedback = 'textbar';
cfg.grid.resolution = 2;
dipM50 = ft_dipolefitting(cfg, avg_trials);
cfg.latency = [0.100 0.120]; % specify latency window around M100 peak
dipM100 = ft_dipolefitting(cfg, avg_trials);
