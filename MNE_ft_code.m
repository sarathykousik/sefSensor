cd('J:\MEG_Research\SEF\SEF-TSSS\014_NJO\sets')
mri_fif = ft_read_mri('_stergaard_niels_j_rgen_4061_01-erikj-140128.fif');
mri_mm     = ft_convert_units(mri_fif, 'mm');

% cfg=[];
% cfg.method = 'interactive';
% mri_ral = ft_volumerealign(cfg, mri_1);

cfg            = [];
cfg.dim        = [256 256 256];
mrirs          = ft_volumereslice(cfg, mri_mm);
mrirs_cm       = ft_convert_units(mrirs, 'cm');

cfg           = [];
cfg.coordsys  = 'neuromag';
cfg.output    = {'skullstrip' 'brain'};
seg           = ft_volumesegment(cfg, mrirs);

cfg              = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg,seg_cm);

cfg           = [];
cfg.method    = 'singleshell';
vol           = ft_prepare_headmodel(cfg,seg);
vol_cm        = ft_convert_units(vol, 'cm');

figure; hold on     % plot all objects in one figure
ft_plot_vol(vol_cm, 'edgecolor', 'none')
alpha 0.4           % make the surface transparent
ft_plot_sens(proc.dataArtRej.grad, 'style', '*b');

ft_plot_mesh(leadfield.pos(grid_gamma_pow,:));


%% saving for MNE
% seg.transform = mri_tal.transform;

cfg             = [];
cfg.filename    = '014_NJO';
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mrirs);

cfg.filename    = '014_NJO_masked';
ft_volumewrite(cfg, seg);

%%
mri = ft_read_mri('J:\MNE_files\orig.nii');
cfg = [];
% cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);

%%

bnd = ft_read_headshape('j:/014_NJO/bem/014_NJO-oct-6-src.fif', 'format', 'mne_source');
figure;ft_plot_mesh(bnd);

%%

mri_nom = ft_read_mri('j:/014_NJO/mri/014_NJO.nii');
mri_nom_cm = ft_convert_units(mri_nom, 'cm');

mri_brain = ft_read_mri('014_NJO_masked.nii');

%%
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';
mri_nom_realign = ft_volumerealign(cfg, mri_nom_cm);

% mri_nom_ctf = ft_convert_units(mri_nom_ctf, 'cm');
T   = mri_nom_realign.transform*inv(mri_nom_realign.transformorig);
bnd  = ft_read_headshape('j:/014_NJO/bem/014_NJO-oct-6-src.fif', 'format', 'mne_source');
sourcespace = ft_convert_units(bnd, 'cm');
% T= eye(4,4);
sourcespace_tr = ft_transform_geometry(T, sourcespace);

%%
% 
cfg           = [];
cfg.coordsys  = 'spm'; 
cfg.output    = {'brain'};
cfg.unit      = mri_nom_realign.unit;
seg           = ft_volumesegment(cfg, mri_nom_realign);
seg_cm        = ft_convert_units(seg,'cm');

cfg           = [];
cfg.method    = 'singleshell';
vol           = ft_prepare_headmodel(cfg,seg);
vol.bnd       = ft_transform_geometry(T, vol.bnd);
vol_cm        = ft_convert_units(vol,'cm');

figure;hold on;
ft_plot_vol(vol);alpha 0.5;
% cfg = [], cfg.parameter = 'brain', ft_sourceplot(cfg,seg_cm)
ft_plot_mesh(sourcespace_tr, 'edgecolor', 'none'); camlight 
ft_plot_sens(ft_transform_geometry(T,proc.timelocked_data.grad), 'style', '*b');

new_grad = ft_transform_geometry(T,proc.timelocked_data.grad);

%%
cfg             = [];
cfg.grad        = ft_transform_geometry(T,proc.timelocked_data.grad);                      % sensor positions
cfg.grid.pos    = sourcespace_tr.pnt;              % source points
cfg.grid.inside = 1:size(sourcespace_tr.pnt,1); % all source points are inside of the brain
cfg.vol         = vol;                               % volume conduction model
leadfield       = ft_prepare_leadfield(cfg);

%%
% cfg = []
% cfg.channel = 'MEG';
% proc.timelocked_data = ft_timelockanalysis(cfg, proc.timelocked_data);
cfg                  = [];
cfg.covariance       = 'yes';
cfg.channel          = 'MEG';
cfg.removemean       = 'no';
cfg.covariancewindow = [-0.09 -0.01];
covariance           = ft_timelockanalysis(cfg, proc.timelocked_data);

%%
cfg                     = [];
cfg.method              = 'mne';
cfg.grid                = leadfield;
cfg.vol                 = vol;
% cfg.mne.prewhiten       = 'yes';
cfg.mne.lambda          = 2;
% cfg.mne.normalize       = 'yes';
% cfg.mne.normalizeparam  = 1e4;
% cfg.mne.noisecov        = covariance.cov;
source_30               = ft_sourceanalysis(cfg,proc.timelocked_data);

toi                     =  find(proc.timelocked_data.time>0.02 & proc.timelocked_data.time<0.030);
bnd.pnt                 = sourcespace.pnt;
bnd.tri                 = sourcespace.tri;
m                       = mean(source_30.avg.pow(:,toi),2); 
ft_plot_mesh(bnd,'vertexcolor', m, 'facealpha', 0.4,'facecolor', 'skin', 'edgecolor', 'none');

%%

cfg        = [];
cfg.method = 'lcmv';
cfg.grid   = leadfield;
cfg.vol    = vol_cm;
cfg.lcmv.lambda = '5%';
cfg.lcmv.keepfilter = 'yes';
cfg.frequency = 15;
sourcepst  = ft_sourceanalysis(cfg, workData{1,1}.saveData.TFR_hann);
% sourcepre  = ft_sourceanalysis(cfg, tlckavgpre);

% sourcepst.avg.nai = sourcepst.avg.pow./sourcepre.avg.pow;

cfg              = [];
% cfg.funparameter = 'nai';
cfg.method       = 'ortho';
cfg.location     = 'max';
figure;
ft_sourceplot(cfg, sourcepst);


%%
cfg = [];
cfg.projectmom = 'yes';
sdFC = ft_sourcedescriptives(cfg,source_30);
sdFC.tri = sourcespace.tri;

cfg = [];
cfg.mask = 'avg.pow';
ft_sourcemovie(cfg,sdFC);

%%
templatefile = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\external\spm8\templates/T1.nii';
template_mri = ft_read_mri(templatefile);

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'avg.pow';
cfg.interpmethod = 'nearest';
cfg.coordsys     = 'mni';
source_diff_int  = ft_sourceinterpolate(cfg, source_30, template_mri);

cfg               = [];
cfg.method        = 'ortho';
cfg.coordsys      = 'mni';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0.0 1.2]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg,source_diff_int);
