%% No morph
bnd  = ft_read_headshape('j:/MNE_files/014_NJO-oct-6-src-nomorph.fif', 'format', 'mne_source' );
sourcespace = ft_convert_units(bnd, 'cm');

cfg                    = [];
cfg.method             = 'singleshell';
hdm                    = ft_prepare_headmodel(cfg, sourcespace);

cfg                    = [];
cfg.spheremesh         = length(hdm.bnd.pnt);
cfg.vol                = hdm;
grid_vol               = ft_prepare_sourcemodel(cfg);

figure, ft_plot_mesh(sourcespace, 'edgecolor', 'k'), hold on
ft_plot_vol(vol_cm, 'facecolor', 'skin');alpha 0.5;
ft_plot_mesh(sourcespace.pos(sourcespace.inside,:));
ft_plot_sens(grad)

%% Morphed
bnd_norm  = ft_read_headshape('j:/MNE_files/014_NJO-oct-6-src.fif', 'format', 'mne_source');
sourcespace_norm = ft_convert_units(bnd_norm, 'cm');

cfg                    = [];
cfg.method             = 'singleshell';
hdm_norm               = ft_prepare_headmodel(cfg, sourcespace_norm);

cfg                    = [];
cfg.spheremesh         = length(hdm_norm.bnd.pnt);
cfg.vol                = hdm_norm;
grid_vol_norm          = ft_prepare_sourcemodel(cfg);

figure, ft_plot_mesh(sourcespace_norm), hold on
ft_plot_vol(hdm_norm, 'facecolor', 'none');alpha 0.5;
ft_plot_mesh(grid_vol_norm.pos(grid_vol_norm.inside,:));
ft_plot_sens(grad)


ft_plot_mesh(grid_vol_norm.pos(grid_vol_norm.inside,:), 'vertexcolor', 'r');
hold on,ft_plot_mesh(grid_vol.pos(grid_vol.inside,:), 'vertexcolor', 'b');

%%

bnd_fsav       = ft_read_headshape('j:/MNE_files/fsaverage-ico-6-src.fif', 'format', 'mne_source');
sourcespace_fsav = ft_convert_units(bnd_norm, 'cm');

cfg                    = [];
cfg.method             = 'singleshell';
hdm_fsav              = ft_prepare_headmodel(cfg, sourcespace_fsav);
% 
% cfg                    = [];
% cfg.spheremesh         = length(hdm_fsav.bnd.pnt);
% cfg.vol                = hdm_fsav;
% grid_vol_fsav          = ft_prepare_sourcemodel(cfg);

figure, ft_plot_mesh(sourcespace_fsav), hold on
ft_plot_vol(hdm_fsav, 'facecolor', 'none');alpha 0.5;
ft_plot_mesh(grid_vol_fsav.pos(grid_vol_fsav.inside,:));
ft_plot_sens(grad)


ft_plot_mesh(grid_vol_fsav.pos(grid_vol_fsav.inside,:));
hold on,ft_plot_mesh(grid_vol.pos(grid_vol.inside,:), 'vertexcolor', 'r');


%% Try to warp

cfg = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = grid_vol_fsav;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.spheremesh = length(hdm.bnd.pnt);
cfg.vol = hdm;
grid_vol               = ft_prepare_sourcemodel(cfg);



%%
cfg   = [];
cfg.parameter = 'pnt';
cfg.interpmethod = 'smudge';
[interp] = ft_sourceinterpolate(cfg, sourcespace , sourcespace_fsav);















