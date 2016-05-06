cd('J:\MEG_Research\SEF\SEF-TSSS\014_NJO\slices')
mri_fif = ft_read_mri('_stergaard_niels_j_rgen_4061_01-erikj-140128.fif');
mri_mm     = ft_convert_units(mri_fif, 'mm');

cfg            = [];
cfg.dim        = [256 256 256];
mrirs          = ft_volumereslice(cfg, mri_mm);

cfg=[];
cfg.method = 'interactive';
cfg.coordsys = mrirs.coordsys;
mri_ral = ft_volumerealign(cfg, mrirs);


cfg           = [];
cfg.coordsys  = 'neuromag';
cfg.output    = {'skullstrip' 'brain'};
seg           = ft_volumesegment(cfg, mri_ral);

% cfg = [];
% cfg.template = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\external\spm8\templates\T1.nii';
% cfg.parameter = 'anatomy';
% cfg.coordsys = 'neuromag';
% [mri] = ft_volumenormalise(cfg, vol_cm)
% 

% seg.transform = mri_ral.transform;
% 
% cfg              = [];
% cfg.funparameter = 'brain';
% ft_sourceplot(cfg,seg);

% saving for MNE

cfg             = [];
cfg.filename    = '014_NJO';
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_ral);

cfg.filename    = '014_NJO_masked';
ft_volumewrite(cfg, seg);

%
mri = ft_read_mri('J:\MNE_files\014_NJO_masked.nii');
cfg = [];
% cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);

%%

bnd = ft_read_headshape('j:/MNE_files/014_NJO-oct-6-src.fif', 'format', 'mne_source');
figure;ft_plot_mesh(bnd);

high_tess_l = ft_read_headshape('J:\old MNE files\014_NJO\surf\lh.white');
high_tess_r = ft_read_headshape('J:\old MNE files\014_NJO\surf\rh.white');
figure;ft_plot_mesh(high_tess_r, high_tess_l);

%%

mri_nom = ft_read_mri('j:/MNE_files/014_NJO.nii');
mri_nom_cm = ft_convert_units(mri_nom, 'cm');

mri_filled = ft_read_mri('J:\old MNE files\014_NJO\mri\filled.nii');

%%
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';
mri_nom_realign = ft_volumerealign(cfg, mri_nom_cm);

%% 

T   = mri_nom_realign.transform*inv(mri_nom_realign.transformorig);
bnd  = ft_read_headshape('j:/MNE_files/014_NJO-oct-6-src.fif', 'format', 'mne_source');
sourcespace = ft_convert_units(bnd, 'cm');
% T= eye(4,4);
sourcespace_tr = ft_transform_geometry(T, sourcespace);



%% 
cfg           = [];
cfg.coordsys  = 'spm'; 
cfg.output    = {'brain'};
cfg.unit      = mri_nom_realign.unit;
seg           = ft_volumesegment(cfg, mri_nom_realign);
seg_cm        = ft_convert_units(seg,'cm');

cfg           = [];
cfg.method    = 'singleshell';
% cfg.grad      = proc.data_epoched.grad;
vol           = ft_prepare_headmodel(cfg,seg_cm);

vol.bnd       = ft_transform_geometry(T, vol.bnd);
vol_cm        = ft_convert_units(vol,'cm');

figure;hold on;
ft_plot_vol(hdm);alpha 0.5;
ft_plot_mesh(grid.pos(grid.inside,:), 'edgecolor', 'none'); camlight 
ft_plot_sens(grad, 'style', '*b');

new_grad = ft_transform_geometry(T,proc.timelocked_data.grad);

%% warping to mni space

templatedir  = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel';
template     = ...
    'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel\standard_sourcemodel3d10mm';
template_cortex =  ft_read_headshape(fullfile(templatedir, 'cortex_8196.surf.gii')); 
ft_plot_mesh(template_cortex, 'facecolor', 'brain', 'edgecolor', 'none')

cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.vol            = vol_cm;
sourcemodel        = ft_prepare_sourcemodel(cfg);  

%%
cfg             = [];
cfg.grad        = proc.timelocked_data.grad;                      % sensor positions
% cfg.grid.pos    = sourcespace.pnt;              % source points
% cfg.grid.inside = 1:size(sourcespace.pnt,1); % all source points are inside of the brain
cfg.grid.inside = 1:size(grid.inside,2);
cfg.grid.pos    = grid.pos(grid.inside,:);
cfg.vol         = hdm;                               % volume conduction model
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
cfg = [];
cfg.channel = 'MEG';
proc.timelocked_data = ft_preprocessing(cfg, proc.timelocked_data);

cfg = [];

avg_tlck = ft_timelockanalysis(cfg, proc.timelocked_data);

%%
cfg                     = [];
cfg.method              = 'mne';
% cfg.grid                = grid;
cfg.grid.leadfield      = leadfield.leadfield;
cfg.grid.inside = 1:size(grid.inside,2);
cfg.grid.pos    = grid.pos(grid.inside,:);
cfg.vol                 = hdm;
cfg.keepfilter = 'yes';
% cfg.mne.prewhiten       = 'yes';
cfg.mne.lambda          = 2;
% cfg.mne.normalize       = 'yes';
% cfg.mne.normalizeparam  = 1e4;
% cfg.mne.noisecov        = covariance.cov;
source_11               = ft_sourceanalysis(cfg,workData{1,1}.saveData.tlck_avg);
source_41               = ft_sourceanalysis(cfg,workData{4,1}.saveData.tlck_avg);
source_71               = ft_sourceanalysis(cfg,workData{7,1}.saveData.tlck_avg);

toi                     =  find(proc.timelocked_data.time>0.02 & proc.timelocked_data.time<0.03);
% grid.pnt                 = grid.pos(grid.inside,:);
% grid.tri                 = grid.tri;

bnd.pnt = grid.pnt

[pnt, tri] = icosahedron162;
xyz.pnt = pnt;
xyz.tri = tri;


m_11                       = mean(source_11.avg.pow(:,toi),2); 
m_41                       = mean(source_41.avg.pow(:,toi),2); 
m_71                       = mean(source_71.avg.pow(:,toi),2); 
sourcemodel.pnt = sourcemodel.pos;

figure, ft_plot_mesh(hdm.bnd,'vertexcolor', [(1:8196)./8196]');

figure, ft_plot_mesh(sourcespace,'vertexcolor', m_41, 'facealpha', 0.8,...
    'facecolor', 'skin', 'edgecolor', 'none');

figure, ft_plot_mesh(sourcespace,'vertexcolor', m_71, 'facealpha', 0.8,...
    'facecolor', 'skin', 'edgecolor', 'none');

%
cfg = [];
cfg.projectmom = 'yes';
sdFC = ft_sourcedescriptives(cfg,source_11);
% sdFC.tri = sourcespace.tri;

cfg = [];
cfg.mask = 'avg.pow';
figure, ft_sourcemovie(cfg,sdFC);

%% DICS

cfg                   = [];
cfg.frequency         = 79;
cfg.grad              = workData{1,1}.saveData.TFR_hann.grad;
cfg.method            = 'dics';
cfg.keeptrials        = 'yes';
cfg.grid              = leadfield;
cfg.latency           = 0.023; 
cfg.vol               = vol_cm;
% cfg.keeptrials        = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.fixedori     = 'yes';
cfg.dics.realfilter   = 'yes';
source_act                = ft_sourceanalysis(cfg,saveData.TFR_hann);
cfg.latency           = -0.02; 
source_bsl            = ft_sourceanalysis(cfg,saveData.TFR_hann);

figure, ft_plot_mesh(sourcespace,'vertexcolor', (1-(source_act.avg.pow./source_bsl.avg.pow)), 'facealpha',1,...
    'facecolor', 'skin', 'edgecolor', 'none');

%%

cfg        = [];
cfg.method = 'sam';
cfg.grid   = leadfield;
% cfg.lcmv.keepfilter = 'yes';
cfg.vol               = vol_cm;
cfg.lambda = '5%';
sourcepst  = ft_sourceanalysis(cfg, workData{1,1}.saveData.tlck_avg);

%%
cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.vol            = vol_cm;
sourcemodel        = ft_prepare_sourcemodel(cfg);  

cfg          = [];
cfg.method   = 'singleshell';
template_vol = ft_prepare_headmodel(cfg, template_seg);
template_vol = ft_convert_units(template_vol, 'cm'); % Convert the vol to cm, since the 

%%

template_atlas = ...
    ft_read_atlas('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
templatedir  = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel';
template_sourcemodel = ft_read_headshape(fullfile(templatedir, 'cortex_8196.surf.gii')); 


cfg = []; 
cfg.interpmethod = 'nearest'; 
cfg.parameter = 'tissue'; 
sourcemodel2 = ft_sourceinterpolate(cfg,template_atlas,sourcespace);

figure, ft_plot_mesh(sourcemodel2, 'facealpha',1,'facecolor', 'skin');
figure, ft_plot_mesh(sourcespace,  'facealpha',1,'facecolor', 'skin');


%%

template = ft_read_mri('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip/external/spm8/templates/T1.nii');
template.coordsys = 'spm'; % so that FieldTrip knows how to interpret the coordinate system
 
cfg          = [];
template_seg = ft_volumesegment(cfg, template);
 
cfg          = [];
cfg.method   = 'singleshell';
template_vol = ft_prepare_headmodel(cfg, template_seg);
template_vol = ft_convert_units(template_vol, 'cm'); % Convert the vol to cm, since the grid will also be expressed in cm
 
% construct the dipole grid in the template brain coordinates
% the source units are in cm
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
load('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel/standard_sourcemodel3d10mm.mat')

cfg = [];
cfg.grid.xgrid  = sourcemodel.xgrid;
cfg.grid.ygrid  = sourcemodel.ygrid;
cfg.grid.zgrid  = sourcemodel.zgrid;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.vol        = template_vol;
template_grid  = ft_prepare_sourcemodel(cfg);
 
% make a figure with the template head model and dipole grid
figure
hold on
ft_plot_vol(template_cortex, 'facecolor', 'cortex');alpha 0.5; 
ft_plot_mesh(template_cortex.pos(template_cortex.inside,:), 'vertexcolor', 'r');
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:))

%%
cd('J:\MEG_Research\SEF\SEF-TSSS\014_NJO\sets')
mri_fif = ft_read_mri('_stergaard_niels_j_rgen_4061_01-erikj-140128.fif');
mri_mm     = ft_convert_units(mri_fif, 'mm');

cfg            = [];
cfg.dim        = [256 256 256];
mrirs          = ft_volumereslice(cfg, mri_mm);

cfg=[];
cfg.method = 'interactive';
cfg.coordsys = mrirs.coordsys;
mri_ral = ft_volumerealign(cfg, mrirs);

cfg = [];
% cfg.downsample = 2;
% cfg.coordsys = 'ctf';
seg = ft_volumesegment(cfg, mri_ral);
seg_cm     = ft_convert_units(seg, 'cm');
seg_cm.anatomy  =  mri_ral.anatomy;

cfg = [];
cfg.method = 'singleshell';
% cfg.grad = grad;
hdm = ft_prepare_headmodel(cfg, sourcespace);

cfg = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
% cfg.spheremesh = length(hdm.bnd.pnt);
% cfg.vol = hdm;
cfg.mri            = seg_cm;
grid               = ft_prepare_sourcemodel(cfg);

cfg = [];
% cfg.grid.warpmni   = 'yes';
% cfg.grid.template  = template_grid;
% cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.spheremesh = length(hdm.bnd.pnt);
cfg.vol = hdm;
% cfg.mri            = seg_cm;
grid_vol               = ft_prepare_sourcemodel(cfg);


cfg = [];
cfg.vol = hdm_mri;
cfg.spheremesh = length(hdm_mri.bnd.pnt);
grid_mri               = ft_prepare_sourcemodel(cfg);

figure;
ft_plot_vol(hdm); hold on
ft_plot_mesh(grid_vol.pos(grid_vol.inside,:));
hold on, ft_plot_mesh(grid.pos(grid.inside,: ), 'vertexcolor', 'r' )

%%
% atlas = ft_read_atlas('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\atlas\spm_anatomy\AllAreas_v18_MPM')
atlas = ft_read_atlas('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
atlas_cm = ft_convert_units(atlas, 'cm');
% load the template sourcemodel with the resolution you need (i.e. the resolution you used in your beamformer grid)
% load('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel/standard_sourcemodel3d10mm.mat')
grid_mm     = ft_convert_units(grid, 'mm');
template_cortex_cm = ft_convert_units(template_cortex, 'cm');
% And call ft_sourceinterpolate: 
cfg = []; 
% cfg.interpmethod = 'spline'; 
cfg.parameter = 'pnt'; 
sourcemodel2 = ft_sourceinterpolate(cfg, template_cortex_cm, sourcespace); 
sourcemodel2_cm = ft_convert_units(sourcemodel2, 'cm');
figure,
ft_plot_mesh(grid.pos(grid.inside,:),'vertexcolor', 'r', 'edgealpha', 0.2), hold on
ft_plot_vol(hdm,'facealpha',0.2,'vertexcolor', 'none')

cfg = [];
cfg.atlas = atlas_cm;
cfg.roi='Frontal_Sup_L';
cfg.inputcoord = 'mni';
mas = ft_volumelookup(cfg,grid);

cfg= [];
cfg.atlas = atlas_cm;
cfg.roi='Frontal_Sup_L';
cfg.coordsys = 'mni';
cfg.inputcoord = 'mni';
% cfg.funparameter = 'brick0';
cfg.maskparameter = mas;
ft_sourceplot(cfg, grid);

%%
morph_map = mne_read_morph_map('j:/MNE_files/014_NJO-oct-6-src_old.fif')

bnd  = ft_read_headshape('j:/MNE_files/014_NJO-oct-6-src_old.fif', 'format', 'mne_source' );
sourcespace = ft_convert_units(bnd, 'cm');

bnd_norm  = ft_read_headshape('j:/MNE_files/014_NJO-oct-6-src.fif', 'format', 'mne_source' );
sourcespace_norm = ft_convert_units(bnd_norm, 'cm');

figure, ft_plot_mesh(sourcespace_norm, 'edgecolor', 'r'), hold on
ft_plot_vol(hdm_mri, 'facecolor', 'skin');alpha 0.5;
ft_plot_sens(grad)

figure, ft_plot_mesh(sourcespace, 'edgecolor', 'none'), hold on
ft_plot_vol(hdm_mri, 'facecolor', 'none');alpha 0.5;
ft_plot_sens(grad)

%% MNI warping
load('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\sourcemodel/standard_sourcemodel3d10mm.mat');
% 
% cfg = [];
% cfg.grid.warpmni   = 'yes';
% cfg.grid.template  = sourcemodel
% cfg.grid.nonlinear = 'yes'; % use non-linear normalization
% grid               = ft_prepare_sourcemodel(cfg, vol);

cfg = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.xgrid  = sourcemodel.xgrid;
cfg.grid.ygrid  = sourcemodel.ygrid;
cfg.grid.zgrid  = sourcemodel.zgrid;
cfg.grid.template = sourcemodel;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.vol        = vol;
grid  = ft_prepare_sourcemodel(cfg);

%% Atlas lookup

cfg=[];
cfg.atlas = 'C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip\template\atlas\afni\TTatlas+tlrc.HEAD'
atlas = ft_prepare_atlas(cfg)
atlas_cm = ft_convert_units(atlas, 'cm');

cfg = [];
cfg.atlas = atlas_cm;
cfg.roi='Postcentral Gyrus';
cfg.inputcoord = 'mni';
mas = ft_volumelookup(cfg,grid);


cfg                 = [];
% cfg.atlas           = atlas_cm;
% cfg.roi             ='Brodmann area 3';
% cfg.coordsys        = 'tal';
% cfg.inputcoord      = 'tal';
% cfg.funparameter    = 'brick1';
% cfg.maskparameter = mas;
ft_plot_mesh(grid.pos(grid.inside,:));
hold on,
ft_plot_vol(vol_cm, 'facealpha', 0.4, 'edgecolor', 'none', 'facecolor', 'skin')
hold on,
ft_plot_mesh(grid.pos(find(mas==1),:), 'vertexcolor', 'r');


