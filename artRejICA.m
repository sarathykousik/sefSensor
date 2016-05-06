
%%
cfg                     = [];
cfg.parameter = 'avg';
cfg.colorbar = 'no';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.03];
% cfg.zlim                = [0 4.5e-12];
cfg.channel = {'all', '-MEG0242+0243', '-MEG1512+1513', '-MEG2612+2613'};
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.colorbar = 'no';

figure
for loop = 1:6
    subplot(2,3, loop)
    ft_topoplotER(cfg,tlckCmb{2,loop})
end

%%

 set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...      % figure handle
    ['Subj-2 N30'],... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', ...           % file format
    '-r250' );             % resolution in dpi
%%

data_IC = ft_appenddata([],proc.data_epoched, proc.data_epoched_ECG)

%%
cfg_sel = []; cfg_sel.channel = 'MEGGRAD'; 
cfg                  = [];
cfg.channel          = 'MEG';
cfg.method           = 'runica';
cfg.numcomponent     = 80;
comp                 = ft_componentanalysis(cfg,ft_selectdata(cfg_sel,data_IC));

%%
% plot the components for visual inspection
cfg                 = [];
cfg.component       = [1:40];
cfg.layout          = 'neuromag306planar.lay';
ft_topoplotIC(cfg, comp)

%%

cfg                 = [];
cfg.layout          = 'neuromag306mag.lay';
cfg.viewmode        = 'component'
ft_databrowser(cfg, comp)

%%

cfg = [];
cfg.component = [1]
data = ft_rejectcomponent(cfg, comp, dataArtRej)

%
cfg = [];
% cfg.channel = 'MEGGRAD';
tlck_data = ft_timelockanalysis(cfg, data);

%%
cfg                     = [];
cfg.parameter           = 'avg';
cfg.colorbar            = 'no';
cfg.layout              = 'neuromag306cmb.lay';
% cfg.channel             = {'all','-MEG0242+0243'};
cfg.xlim                = [0 0.1];
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
cfg.colorbar            = 'yes';
ft_multiplotER(cfg, ft_combineplanar([],tlck{2,6}))










