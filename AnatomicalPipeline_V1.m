clc; close all; clear all

%% Load registered Neuromag MRI
% cd('J:\MEG_Research\SEF-SourceAnalysis\015_IHM_MRI\slices')
mri_fif = ft_read_mri('015_IHM.nii');
% ft_determine_coordsys(mri_fif, 'interactive', 'yes');
% mri_fif_cm     = ft_convert_units(mri_fif, 'cm');
mri_fif_cm     = ft_convert_units(mri_fif, 'cm'); 

%% Reslice, segment

cfg            = [];
cfg.dim        = [256 256 256];
mrirs          = ft_volumereslice(cfg, mri_fif_mm);
% mrirs_cm        = ft_convert_units(mrirs, 'cm');

%% Re-align to Talairach coordinates
cfg        = [];
cfg.parameter = 'anatomy';

cfg.coordsys  = 'neuromag';
% cfg.fiducial.rpa  = [86.7   0    0]; 
% cfg.fiducial.nas  = [0    115.3  0]; 
% cfg.fiducial.lpa  = [-73.5  0    0]; 
mri_neuromag    = ft_volumerealign(cfg, mrirs);
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



%% Save Tal-MRI files
cfg             = [];
cfg.filename    = '015_IHM';
cfg.filetype    = 'mgz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mrirs_cm);

cfg.filename    = '014_NJO_masked';
ft_volumewrite(cfg, seg);

%%

%  Freesurfer pipeline
%  MNE pipeline


%% Load cortical reconstruction
bnd = ft_read_headshape('015_IHM-oct-6-src.fif', 'format', 'mne_source');
bnd_cm = ft_convert_units(bnd, 'cm');
figure;ft_plot_mesh(bnd);

%% Load Tal MRI and register back to MNI
% mri_nom           = ft_read_mri('j:\MEG_Research\SEF\SEF-TSSS\014_NJO.mgh');
% mri_nom.coordsys  = 'tal';
cfg               = [];
cfg.fiducial.rpa  = [86.7   0    0]; 
cfg.fiducial.nas  = [0    115.3  0]; 
cfg.fiducial.lpa  = [-73.5  0    0]; 
mri_nom_neuromag  = ft_volumerealign(cfg, seg);

mri_nom_neuromag.coordsys = 'neuromag';
mri_nom_cm        = ft_convert_units(mri_nom_neuromag, 'cm');

%% Segment original Neuromag MRI - create headmodel from it
%  Load sourcespace

% mri_nom_neuromag = ft_convert_units(mri_nom_neuromag, 'cm');
% T   = mri_nom_cm.transform * inv(mri_nom_cm.transformorig);

sourcespace = ft_convert_units(bnd, 'cm');
sourcespace_T = ft_transform_geometry(inv(T), sourcespace);

% cfg           = [];
% cfg.output    = {'brain'};
% seg           = ft_volumesegment(cfg, mri_nom_neuromag);
% seg_cm        = ft_convert_units(seg, 'cm');

%% seg on tal
cfg           = [];
cfg.coordsys  = 'neuromag';
cfg.output    = {'skullstrip' 'brain'};
seg_tal       = ft_volumesegment(cfg, mri_tal);

%%
cfg           = [];
cfg.method    = 'singleshell';
vol       = ft_prepare_headmodel(cfg,seg_cm);

%%

vol_T = ft_transform_geometry((mri_neuromag.transform)*inv(mri_tal.transformorig), vol);
% vol_cm_T        = ft_convert_units(vol_T, 'm');

figure;hold on;
ft_plot_vol(vol, 'facecolor', 'none');
% ft_plot_mesh(sourcespace_T, 'edgecolor', 'none'); camlight 
ft_plot_sens(grad_tsss3, 'style', '*b')

% vol_T         = vol;
% vol_T.bnd       = ft_transform_geometry(T, vol_T.bnd);


%% Import fif file
file_name                   = dir('*.fif')

cfg                         = [];
cfg.dataset                 = file_name(1).name;
data_import                 = ft_preprocessing(cfg);
grad_tsss3                  = data_import.grad;

shape   = ft_read_headshape(file_name(1).name,'unit','cm');
grad    = ft_read_sens(file_name(1).name,'senstype','meg');

%%
figure;
hold on
ft_plot_headshape(shape);
% ft_plot_sens(grad, 'style', '*b');
ft_plot_mesh(bnd_cm);


%% Questions
% Is transformations right?
% Either transform both or don't transform both
% Answer: Don't transform anything. Sensor positions align only when there is no transformation