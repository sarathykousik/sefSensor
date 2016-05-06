%% 
% NOTE: the path to the template file is user-specific
template = ft_read_mri('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip-20150222\external\spm8\templates\T1.nii');
template.coordsys = 'spm'; % so that FieldTrip knows how to interpret the coordinate system
 
% segment the template brain and construct a volume conduction model (i.e. head model): this is needed
% for the inside/outside detection of voxels.
cfg          = [];
template_seg = ft_volumesegment(cfg, template);
 
cfg          = [];
cfg.method   = 'singleshell';
template_vol = ft_prepare_headmodel(cfg, template_seg);
template_vol = ft_convert_units(template_vol, 'cm'); % Convert the vol to cm, since the grid will also be expressed in cm
 
% construct the dipole grid in the template brain coordinates
% the source units are in cm
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.vol        = template_vol;
template_grid  = ft_prepare_sourcemodel(cfg);
 
% make a figure with the template head model and dipole grid
figure
hold on
ft_plot_vol(template_vol, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% Check Seg

cfg = [];
cfg.funparameter = 'csf';
ft_sourceplot(cfg,seg_cm); %only mri

%% MNI warping
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, seg);

cfg = [];
% cfg.grid.warpmni   = 'yes';
% cfg.grid.template  = template_grid;
% cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.vol            = vol;
grid               = ft_prepare_sourcemodel(cfg);
 
% make a figure of the single subject headmodel, and grid positions
figure;
ft_plot_vol(vol, 'edgecolor', 'none'); alpha 0.4;
ft_plot_mesh(grid.pos(grid.inside,:));

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

%%
mri = ft_read_mri('j:\MNE_files\subjects\014_NJO.nii');


cfg           = [];
cfg.output    = {'brain','skull','scalp'};
cfg.coordsys  = 'neuromag';
segmentedmri  = ft_volumesegment(cfg, mri);
seg_cm = ft_convert_units(segmentedmri, 'cm');


cfg=[];
cfg.tissue={'brain', 'skull','scalp'};
cfg.numvertices = [3000 2000 1000];
bnd=ft_prepare_mesh(cfg,seg_cm);

cfg        = [];
cfg.method ='bemcp';
vol        = ft_prepare_headmodel(cfg, bnd);

