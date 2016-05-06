clc; close all; clear all

%% Load registered Neuromag MRI
cd('J:\MEG_Research\SEF-SourceAnalysis\015_IHM_MRI\slices')
mri_fif = ft_read_mri('subjectno_0015_85643521_01.fif');
% ft_determine_coordsys(mri_fif, 'interactive', 'yes');
% mri_fif_cm     = ft_convert_units(mri_fif, 'cm');
mri_fif_mm     = ft_convert_units(mri_fif, 'mm'); 

%% Reslice, segment

cfg            = [];
cfg.dim        = [256 256 256];
mrirs          = ft_volumereslice(cfg, mri_fif_mm);
mrirs_cm        = ft_convert_units(mrirs, 'cm');

%% Re-align to Talairach coordinates
cfg        = [];
cfg.parameter = 'anatomy';

cfg.coordsys  = 'neuromag';
% cfg.fiducial.rpa  = [86.7   0    0]; 
% cfg.fiducial.nas  = [0    115.3  0]; 
% cfg.fiducial.lpa  = [-73.5  0    0]; 
mri_neuromag    = ft_volumerealign(cfg, mrirs);
% mrirs_cm        = ft_convert_units(mri_neuromag, 'cm');

mri_neuromag.transform

%%
cfg           = [];
cfg.coordsys  = 'neuromag';
cfg.output    = {'skullstrip' 'brain'};
seg           = ft_volumesegment(cfg, mri_neuromag);
seg_cm        = ft_convert_units(seg, 'cm');

% cfg =[];
% cfg.funparameter = 'brain';
% ft_sourceplot(cfg,seg);

%%
cfg           = [];
cfg.method    = 'singleshell';
vol       = ft_prepare_headmodel(cfg,seg_cm);

%%

% vol_orig = ft_transform_geometry((mri_neuromag_2.transform)*inv(mrirs_cm.transform), vol);
% vol_cm_T        = ft_convert_units(vol_T, 'm');

figure;hold on;
ft_plot_vol(vol, 'facecolor', 'none');
% ft_plot_mesh(sourcespace_T, 'edgecolor', 'none'); camlight 
ft_plot_sens(grad_tsss3, 'style', '*b')

% vol_T         = vol;
% vol_T.bnd       = ft_transform_geometry(T, vol_T.bnd);



%% Save Tal-MRI files
cfg             = [];
cfg.filename    = '015_IHM';
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mrirs);

cfg.filename    = '015_IHM_masked';
ft_volumewrite(cfg, seg);

%%

%  Freesurfer pipeline
%  MNE pipeline




%% Import fif file
file_name                   = dir('*.fif')
asfasdfasdfasdfasdf
cfg                         = [];
cfg.dataset                 = file_name(1).name;
data_import                 = ft_preprocessing(cfg);
grad_tsss3                  = data_import.grad;

%% orig
cfg           = [];
cfg.coordsys  = 'neuromag';
cfg.output    = {'skullstrip' 'brain'};
segorig           = ft_volumesegment(cfg, mrirs_cm);
segorig_cm        = ft_convert_units(seg, 'cm');

cfg           = [];
cfg.method    = 'singleshell';
volorig       = ft_prepare_headmodel(cfg,segorig_cm);

figure;hold on;
ft_plot_vol(vol_orig, 'facecolor', 'none');
ft_plot_sens(grad_tsss3, 'style', '*b')

%% Questions
% Is transformations right?
% Either transform both or don't transform both
% Answer: Don't transform anything. Sensor positions align only when there is no transformation