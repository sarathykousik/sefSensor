%% Beamforming pipeline
% Works subject-wise
% Repeat after anatomical pipeline for each subject/condition

%%      %%     Data processing pipeline - for beaming Gamma        %%    %%
% Load clean data, compute Gamma (Fourier coeff) for beaming
% Split data: baseline = -0.90 to -0.10 s; data = 0.20 to 0.60 s 
% Subject-wise TFR ROI? - else: t= 0.20-0.60 s; Fc=65; BW=25

%%  Loop over each condition  %%

cfg                     = [];
cfg.channel             = {'MEG*2', 'MEG*3'};
data                    = ft_selectdata(cfg, visClean); 

cfg                     = [];                                           
cfg.toilim              = [-0.09 -0.01];                       
data_bsl                = ft_redefinetrial(cfg, data);
     
cfg.toilim              = [0.020 0.060];                       
data_exp                = ft_redefinetrial(cfg, data);

cfg                     = [];
data_cmb                = ft_appenddata(cfg, data_bsl, data_exp);
data_cmb.trialinfo      = [zeros(length(data_bsl.trial), 1); ...
                                ones(length(data_exp.trial), 1)];

cfg                     = [];
cfg.paramter            = 'trial';
cfg.keeptrials          = 'yes';
cfg.taper               = 'dpss';
cfg.output              = 'fourier';
cfg.channel             = {'MEG*2', 'MEG*3'};
cfg.keeptrials          = 'yes';
cfg.method              = 'mtmfft';
cfg.pad                 = 1;
cfg.foi                 = 65;                          
cfg.tapsmofrq           = 25;
freq_cmb                = ft_freqanalysis(cfg, data_cmb);

cfg                     = [];
cfg.trials              = freq_cmb.trialinfo == 0;
freq_bsl                = ft_selectdata(cfg, freq_cmb);
freq_bsl.cumtapcnt      = freq_cmb.cumtapcnt(cfg.trials);
freq_bsl.cumsumcnt      = freq_cmb.cumsumcnt(cfg.trials);

cfg.trials              = freq_cmb.trialinfo == 1;
freq_exp                = ft_selectdata(cfg, freq_cmb);
freq_exp.cumtapcnt      = freq_cmb.cumtapcnt(cfg.trials);
freq_exp.cumsumcnt      = freq_cmb.cumsumcnt(cfg.trials);

%% Leadfield - calculated per subject
cfg                     = [];
cfg.grid.pos            = sourcespace.pnt;
cfg.grid.inside         = 1:length(sourcespace.pnt);  
cfg.vol                 = vol;
cfg.channel             = {'MEG*2', 'MEG*3'};
% cfg.grad                = grad;
% cfg.inwardshift         = -1;
cfg.reducerank          = 2;
sourcemodel_lf          = ft_prepare_leadfield(cfg,freq_cmb);

%%
figure;hold on;
% ft_plot_vol(vol, 'facecolor', 'none');alpha 0.5;
ft_plot_mesh(sourcespace, 'edgecolor', 'k','facecolor', 'skin' ), camlight('left')
% ft_plot_mesh(sourcemodel_lf.pos(sourcemodel_lf.inside,:), 'vertexcolor', 'r');  
% ft_plot_sens(freq_bsl.grad, 'style', '*b')

%%
cfg                     = [];
cfg.frequency           = 65;
cfg.method              = 'dics';
% cfg.keeptrials          = 'yes';
cfg.grid.leadfield      = sourcemodel_lf.leadfield;
cfg.grid.pos            = sourcemodel_lf.pos(sourcemodel_lf.inside,:);
cfg.grid.inside         = 1:length(sourcemodel_lf.inside);  
cfg.vol                 = vol;
cfg.dics.lambda         = '80%';
cfg.dics.keepfilter     = 'yes';
cfg.dics.fixedori       = 'yes';
cfg.dics.realfilter     = 'yes';
cfg.keepleadfield       = 'yes';
source_cmb              = ft_sourceanalysis(cfg,freq_cmb);

% Beam pre- and poststim by using the common filter
%  Loop over each conditions
cfg.grid.filter  = source_cmb.avg.filter;
source_bsl       = ft_sourceanalysis(cfg, freq_bsl);
source_exp       = ft_sourceanalysis(cfg, freq_exp);

source_diff = source_exp;
source_diff.avg.pow = (source_exp.avg.pow ./ source_bsl.avg.pow) ;

bnd.pnt = sourcespace.pnt;
bnd.tri = sourcespace.tri;
figure, ft_plot_mesh(bnd, 'vertexcolor', source_diff.avg.pow), title(num2str(cfg.dics.lambda));

%%
[maxval, maxGammaindx] = max(source_diff.avg.pow);

figure;hold on;
ft_plot_vol(vol, 'facecolor', 'none');alpha 0.5;
% ft_plot_mesh(sourcespace, 'edgecolor', 'k','facecolor', 'skin' ), camlight
ft_plot_mesh(sourcemodel_lf.pos(maxGammaindx,:), 'vertexcolor', 'r');  
ft_plot_sens(freq_cmb.grad, 'style', '*b')

%% Cov/Avg data for LCMV

cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.keeptrials        = 'yes';
tlock                 = ft_timelockanalysis(cfg, data);

%% LCMV on maxGamma location

cfg                 = [];
cfg.method          = 'lcmv';
cfg.lcmv.fixedori   = 'yes' 
cfg.lcmv.lambda     = '5%';
cfg.vol             = vol;
cfg.grid.pos        = sourcemodel_lf.pos([maxGammaindx], :);
cfg.grid.inside     = 1:size(cfg.grid.pos, 1);
cfg.grid.outside    = [];
cfg.keepfilter      = 'yes';
source_maxGamma     = ft_sourceanalysis(cfg, tlock);

maxGammaVirtSens    = [];
for i=1:length(data.trial)
    maxGammaVirtSens.trial{i} = source_maxGamma.avg.filter{1} * data.trial{i};
end
maxGammaVirtSens.label = {'1'};
maxGammaVirtSens.time  = data.time;
maxGammaVirtSens.fsample  = data.fsample;

cfg                   = [];
tlockVirtSens          = ft_timelockanalysis(cfg, maxGammaVirtSens);

plot(tlockVirtSens.time, tlockVirtSens.avg)

%% TFR of virtualchannel
cfg              = [];
cfg.paramter     = 'trial';
cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.pad          = 1;
cfg.taper        = 'hanning';
cfg.foi          = 1:2:90;                          
cfg.t_ftimwin    = 2./cfg.foi;   
cfg.toi          = -0.5:0.010:0.5;           
% cfg.tapsmofrq    = 25;
% cfg.keeptapers   = 'yes';
virtTFR          = ft_freqanalysis(cfg, maxGammaVirtSens);

cfg = [];
cfg.baseline                = [-0.090 -0.01];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';
virtTFR_bsl               = ft_freqbaseline(cfg, virtTFR);

cfg                         = [];
cfg.parameter               = 'powspctrm';
cfg.zlim                    = 'maxmin';
cfg.xlim                    = [-0.01 0.4]; 
% cfg.ylim                    = [60 100];    
cfg.maskstyle               = 'saturation';	
cfg.shading                 = 'interp';
cfg.masknans                = 'yes';
figure,ft_singleplotTFR(cfg, virtTFR_bsl);


