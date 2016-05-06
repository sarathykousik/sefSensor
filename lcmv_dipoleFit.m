%%
% ft_read_mri -> ft_volumesegment -> ft_prepare_headmodel ...
%                 -> ft_prepare_leadfield ->  ft_source_analysis

%%
% Read MRI
cd('J:\data-tsss-sef\SEF-TSSS\014_NJO\sets')
mri_fif = ft_read_mri(...
    '_stergaard_niels_j_rgen_4061_01-erikj-140128.fif');
mri_mm     = ft_convert_units(mri_fif, 'mm');


cfg             = [];
cfg.dim = [256 256 256];
mri_resl = ft_volumereslice(cfg, template);

% cfg                     = [];
% cfg.coordsys            = 'neuromag';    
% % cfg.fiducial.nas        = [0        105.4    0];
% % cfg.fiducial.lpa        = [-73.6     0       0];
% % cfg.fiducial.rpa        = [85.2      0       0];
% realign_mri             = ft_volumerealign(cfg, mri_resl);

cfg                     = [];
[segmentedmri]          = ft_volumesegment(cfg, mri_resl);
segmentedmri.transform  = mri_resl.transform;
segmentedmri.anatomy    = mri_resl.anatomy;

segmentedmri_cm         = ft_convert_units(segmentedmri, 'cm');

cfg                     = [];
cfg.funparameter        = 'gray';
ft_sourceplot(cfg,segmentedmri);

%% Prepare headmodel / volume conductor model
cfg                 = [];
cfg.method          = 'singleshell';
cfg.grad            = proc.timelocked_data.grad;
hdm                 = ft_prepare_headmodel(cfg, segmentedmri); % vol
hdm_cm              = ft_convert_units(hdm, 'cm');
% grad_mm = ft_convert_units(proc.data_import.grad, 'mm');

%%  Source model /grid
% Warp tp MNI space for inter-subject comparison

templatedir         = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel';
template            = load(fullfile(templatedir, 'standard_sourcemodel3d8mm')); 

% inverse-warp the subject specific grid to the template grid 
cfg                = [];
cfg.vol            = hdm;
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template.sourcemodel;
cfg.grid.nonlinear = 'no'; % use non-linear normalization
% cfg.sourceunits    = 'mm';
sourcemodel        = ft_prepare_sourcemodel(cfg);

hdm_cm             = ft_convert_units(hdm, 'cm');
sourcemodel_cm     = ft_convert_units(sourcemodel, 'cm');

figure; hold on     % plot all objects in one figure
ft_plot_vol(hdm_cm, 'edgecolor', 'none')
alpha 0.4           % make the surface transparent
ft_plot_mesh(sourcemodel_cm.pos(sourcemodel_cm.inside,:)); hold on
ft_plot_sens(proc.timelocked_data.grad, 'style', '*b');

%%
cfg             = [];
cfg.grid        = sourcemodel_cm;
cfg.vol         = hdm_cm;
cfg.channel     = {'MEG'};
cfg.grad        = saveData.dataArtRej.grad;
sourcemodel_lf  = ft_prepare_leadfield(cfg, proc.timelocked_data.grad);

%%

cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = [0.015 0.15];
cfg.removemean       = 'no';
tlckavgpst           = ft_timelockanalysis(cfg, proc.timelocked_data);
cfg.covariancewindow = [-0.09 -0.01];
tlckavgpre           = ft_timelockanalysis(cfg, proc.timelocked_data);

cfg        = [];
cfg.method = 'mne';
cfg.grid   = sourcemodel_lf;
cfg.grad   = proc.timelocked_data.grad;
cfg.vol    = hdm_cm;
cfg.mne.lambda = '5%';
sourcepst  = ft_sourceanalysis(cfg, tlckavgpst);
sourcepre  = ft_sourceanalysis(cfg, tlckavgpre);

sourcepst.avg.nai = sourcepst.avg.pow./sourcepre.avg.pow;
sourcepst.pos = template.sourcemodel.pos;
sourcepst.dim = template.sourcemodel.dim;

cfg              = [];
cfg.parameter    = 'avg.pow';
cfg.interpmethod = 'nearest';
source_diff_int  = ft_sourceinterpolate(cfg, sourcepst, segmentedmri_cm);

cfg              = [];
cfg.funparameter = 'avg.nai';
cfg.method       = 'ortho';
cfg.location     = 'max';
ft_sourceplot(cfg, sourcepst);

%%
cfg        = [];                                           
cfg.toilim = [-0.09 -0.01];                       
data_bsl   = ft_redefinetrial(cfg, proc.data_epoched);
     
cfg.toilim = [0.005 0.3];                       
data_exp   = ft_redefinetrial(cfg, proc.data_epoched);

cfg      = [];
data_cmb = ft_appenddata(cfg, data_bsl, data_exp);
% give a number to each trial: 0 = baseline, 1 = experimental condition
data_cmb.trialinfo = [zeros(length(data_bsl.trial), 1); ones(length(data_exp.trial), 1)];

cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 15;
cfg.foi        = 55;
freq_cmb       = ft_freqanalysis(cfg, data_cmb);

cfg                = [];
cfg.trials         = freq_cmb.trialinfo == 0;
freq_bsl           = ft_selectdata(cfg, freq_cmb);
% remember the number of tapers per trial
freq_bsl.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
freq_bsl.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);

cfg.trials         = freq_cmb.trialinfo == 1;
freq_exp           = ft_selectdata(cfg, freq_cmb);
% remember the number of tapers per trial
freq_exp.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
freq_exp.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);

%%
cfg             = [];
cfg.grid        = sourcemodel;
cfg.vol         = hdm;
cfg.channel     = {'MEG'};
cfg.grad        = freq_cmb.grad;
sourcemodel_lf  = ft_prepare_leadfield(cfg, freq_cmb);

cfg                   = [];
cfg.frequency         = freq_cmb.freq;
cfg.grad              = freq_cmb.grad;
cfg.method            = 'dics';
cfg.keeptrials        = 'yes';
cfg.grid              = sourcemodel_lf;
cfg.vol               = hdm;
cfg.keeptrials        = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.fixedori     = 'yes';
cfg.dics.realfilter   = 'yes';
source                = ft_sourceanalysis(cfg, freq_cmb);

% beam pre- and poststim by using the common filter
cfg.grid.filter  = source.avg.filter;
source_bsl       = ft_sourceanalysis(cfg, freq_bsl);
source_exp       = ft_sourceanalysis(cfg, freq_exp);

source_diff = source_exp;
source_diff.avg.pow = (source_exp.avg.pow ./ source_bsl.avg.pow) - 1;

source_diff.pos = template.sourcemodel.pos;
source_diff.dim = template.sourcemodel.dim;

% note that the exact directory is user- and platform-specific
templatefile = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\external\spm8\templates\T1.nii';
template_mri = ft_read_mri(templatefile);

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'avg.pow';
cfg.interpmethod = 'nearest';
cfg.coordsys     = 'mni';
source_diff_int  = ft_sourceinterpolate(cfg, source_diff, template_mri);

%%
cfg               = [];
cfg.method        = 'slice';
cfg.coordsys      = 'mni';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg,source_diff_int);








