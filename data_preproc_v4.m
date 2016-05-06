clc;    close all;    clear all;
ft_defaults;

proc.data_folder                 = 'j:\MEG_Research\SEF\SEF-TSSS\014_NJO' ;
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              = 10;
proc.save_folder                 = 'j:\MEG_Research\SEF\results\014_NJO';

%% Declarations
cd(proc.data_folder);
filename                    = dir('*.fif');
loop                        = 2;
% Import data, Define trials, Cut trials
cfg                         = [];
cfg.dataset                 = filename(loop).name;
proc.data_import            = ft_preprocessing(cfg);

cfg = [];
cfg.trialdef.prestim        = 0.5;
cfg.trialdef.poststim       = 0.5;
cfg.trialdef.eventtype      = 'STI001'; 
cfg.trialdef.eventvalue     = 5;
cfg.dataset                 = filename(loop).name;
proc.data_trial             = ft_definetrial(cfg);

proc.blank_seq          = ones(size(proc.data_import.trial{1},2),1);
for ind = proc.pre_stim_time:proc.post_stim_time
    proc.blank_seq(proc.data_trial.trl(:,1)+(cfg.trialdef.prestim*1000)+ind) = 0;
end

proc.MEG_channel = find(strncmp(proc.data_import.label, 'MEG',3)==1);                    % Pick only MEG channels
proc.data_import.trial{1}(proc.MEG_channel,:) = ...
        bsxfun(@times, proc.data_import.trial{1}(proc.MEG_channel,:), proc.blank_seq');  % Multiply blanking sequence with MEG channels

% Filtering, Baseline correction and extraction of MEG channels only
cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfreq                  =  100;
cfg.hpfilter                = 'yes';
cfg.hpfreq                  =  2;
cfg.dftfilter               = 'yes';
cfg.channel                 =  {'MEG'}; 
cfg.detrend                 = 'yes';
proc.preproc_data_MEG       = ft_preprocessing(cfg, proc.data_import);

cfg.channel                 =  {'EOG'}; 
proc.preproc_data_EOG       = ft_preprocessing(cfg, proc.data_import);

proc.data_epoched           = ft_redefinetrial(proc.data_trial,proc.preproc_data_MEG);
proc.data_epoched_EOG       = ft_redefinetrial(proc.data_trial,proc.preproc_data_EOG);

%% Jump artefact
cfg                    	= [];
cfg.trl                 = proc.data_trial.trl;
cfg.continuous          = 'no';
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel    = 'MEG';
cfg.artfctdef.zvalue.cutoff     = 30;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
% algorithmic parameters
cfg.artfctdef.zvalue.cumulative    = 'yes';
cfg.artfctdef.zvalue.medianfilter  = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff       = 'yes';
% make the process interactive
temp.artifact_jump = ft_artifact_zvalue(cfg, proc.data_epoched);

% EOG
cfg            = [];
cfg.trl        = proc.data_trial.trl;
cfg.continuous = 'no'; 
cfg.trl                 = proc.data_trial.trl;
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = 'EOG*';
cfg.artfctdef.zvalue.cutoff      = 4;
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;
% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';
cfg.artfctdef.zvalue.interactive = 'yes';
[~, temp.artifact_EOG] = ft_artifact_zvalue(cfg, proc.data_epoched_EOG);

cfg = [];
cfg.artfctdef.jump.artifact = temp.artifact_jump.artfctdef.zvalue.artifact;
cfg.artfctdef.eog.artifact  = temp.artifact_EOG;
cfg.artfctdef.reject    = 'complete';
proc.dataArtRej = ft_rejectartifact(cfg,proc.data_epoched);

%%
% cfg = [];
% cfg.channel = {'MEG***2', 'MEG***3'};
% proc.dataArtRej = ft_preprocessing(cfg, proc.dataArtRej);
% proc.dataArtRej = save_data;
% proc.dataArtRej = ft_selectdata(proc.dataArtRej, 'channel', {'MEG***2', 'MEG***3'} );

%%
cfg                         = [];
cfg.output                  = 'pow';
cfg.paramter                = 'trial';
cfg.pad                     = 1;
cfg.method                  = 'mtmconvol';
cfg.taper                   = 'hanning';
cfg.foi                     = 4:4:100;                          
cfg.t_ftimwin               = 5./cfg.foi;  % 7 cycles per time window
cfg.toi                     = -0.5:0.005:0.5;                 
cfg.keeptrials              = 'yes';
TFR_hann                    = ft_freqanalysis(cfg,proc.dataArtRej);

cfg = [];
cfg.baseline                = [-0.090 -0.01];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';
TFR_hann_bsl = ft_freqbaseline(cfg, TFR_hann);

cfg                         = [];
TFR_hann_bsl_cmb                = ft_combineplanar(cfg, TFR_hann_bsl);

cfg                         = [];
cfg.zlim                    = 'maxmin';
cfg.xlim                    = [-0.01 0.4]; 
cfg.ylim                    = [55 100];    
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
figure,ft_multiplotTFR(cfg, TFR_hann_bsl_cmb);

%%  Leadfield

cfg             = [];
cfg.grad        = save_data.grad;                      % sensor positions
cfg.grid.pos    = sourcespace.pnt;              % source points
cfg.grid.inside = 1:size(sourcespace.pnt,1); % all source points are inside of the brain
% cfg.grid.inside = 1:size(grid.inside,2);
% cfg.grid.pos    = grid.pos(grid.inside,:);
cfg.vol         = vol_cm;                               % volume conduction model
leadfield       = ft_prepare_leadfield(cfg);

%%  DICS Beamforming 

cfg           = [];                                           
cfg.toilim    = [-0.09 0.05];           
cfg.minlength = 'maxperlen'; % this ensures all resulting trials are equal length
data          = ft_redefinetrial(cfg, saveData.dataArtRej);

cfg        = [];                                           
cfg.toilim = [-0.09 -0.01];                       
data_bsl   = ft_redefinetrial(cfg, data);
     
cfg.toilim = [0.01 0.050];                       
data_exp   = ft_redefinetrial(cfg, data);

cfg      = [];
data_cmb = ft_appenddata(cfg, data_bsl, data_exp);

% give a number to each trial: 0 = baseline, 1 = experimental condition
data_cmb.trialinfo = [zeros(length(data_bsl.trial), 1); ones(length(data_exp.trial), 1)];

% Gamma - 0.02-0.04 sec; 60-90 Hz; Center Freq = 75 Hz - tapsmofrq = 15
cfg            = [];
cfg.method     = 'mtmfft';
cfg.taper      = 'hanning';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 10;
cfg.foi        = 75;
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

%% Check forward model files
figure; hold on;
ft_plot_vol(vol_cm);alpha 0.5;
ft_plot_mesh(sourcespace, 'edgecolor', 'none'); camlight 
ft_plot_sens(data_cmb.grad, 'style', '*b');

%%
cfg             = [];
cfg.grid        = sourcespace;
cfg.vol         = vol_cm;
cfg.channel     = {'MEG'};
cfg.grad        = data_cmb.grad;
sourcemodel_lf  = ft_prepare_leadfield(cfg);

%%
cfg                   = [];
cfg.frequency         = freq_cmb.freq;
cfg.grad              = freq_cmb.grad;
cfg.method            = 'dics';
cfg.keeptrials        = 'yes';
cfg.grid              = sourcemodel_lf;
cfg.vol               = vol_cm;
cfg.keeptrials        = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
% cfg.dics.fixedori     = 'yes';
cfg.dics.realfilter   = 'yes';
source                = ft_sourceanalysis(cfg, freq_cmb);

cfg.grid.filter  = source.avg.filter;
source_bsl       = ft_sourceanalysis(cfg, freq_bsl);
source_exp       = ft_sourceanalysis(cfg, freq_exp);

source_diff = source_exp;
source_diff.avg.pow = (source_exp.avg.pow./ source_bsl.avg.pow);

bnd.pnt = sourcespace.pnt;
bnd.tri = sourcespace.tri;

% mean_source = nanmean(source_diff.avg.pow);
% sd_source   = nanstd(source_diff.avg.pow);
% up_lt = mean_source+2*sd_source;
% lo_lt = mean_source-2*sd_source;
% ind   = find(source_diff.avg.pow>up_lt);
% source_diff.avg.pow(ind) = 0;

mean_source = nanmean(source_diff.avg.pow);
sd_source   = nanstd(source_diff.avg.pow);
up_lt = mean_source+2*sd_source;
lo_lt = mean_source-2*sd_source;
ind   = find(source_diff.avg.pow>up_lt);

m=zeros(8196,1);
% m=source_diff.avg.pow;
m(ind)=source_diff.avg.pow(ind);
figure,ft_plot_mesh(bnd, 'vertexcolor', m, 'facecolor', 'cortex')

pl = zeros(8196,1);
pl(1860:2100) = m_11(1860:2100);
figure,ft_plot_mesh(bnd, 'vertexcolor', source_diff.avg.pow, 'facecolor', 'cortex')

grid_gamma_pow = find(pl~=0);
plot(grid_gamma_pow, m(grid_gamma_pow), '*')
hold on, plot(pl)

[maxval, maxcohindx] = max(source_exp.avg.pow(280:1068));

%%

cfg                         = [];
cfg.lpfilter                = 'yes';
cfg.lpfreq                  =  60;
cfg.hpfilter                = 'yes';
cfg.hpfreq                  =  85;
for_beaming                 = ft_preprocessing(cfg, saveData.dataArtRej);

cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.keeptrials        = 'yes';
tlock                 = ft_timelockanalysis(cfg, saveData.dataArtRej);

%%
cfg                 = [];
cfg.method          = 'lcmv';
cfg.lcmv.fixedori   = 'yes' 
cfg.lcmv.lambda     = '5%';
cfg.vol             = vol_cm;
cfg.grid.pos        = sourcemodel_lf.pos([grid_gamma_pow], :);
cfg.grid.inside     = 1:size(cfg.grid.pos, 1);
cfg.grid.outside    = [];
cfg.keepfilter      = 'yes';
source_idx          = ft_sourceanalysis(cfg, tlock);

chansel = ft_channelselection('MEG', data.label); % find MEG sensor names
chansel = match_str(data.label, chansel);         % find MEG sensor indices

gam_pow_data = saveData.dataArtRej;
gam_pow_data.label = {'1'};     %, '2', '3', '4','5','6','7','8','9','10','11','12','13'};
% gam_pow_data.time = data_cmb.time;
% for channel = 1:length(beamformer_gam_pow)
% beamformer_gam_pow = source_idx.avg.filter{10};

% end

% timeseries = cat(2, gam_pow_data.trial{:});
% [u, s, v] = svd(timeseries, 'econ');
% timeseriesmaxproj = u(:,1)' * timeseries;


for i=1:length(saveData.dataArtRej.trial)
    gam_pow_data.trial{i} = ...
      source_idx.avg.filter{1} * saveData.dataArtRej.trial{i};
end

cfg                         = [];
cfg.output                  = 'pow';
cfg.paramter                = 'trial';
cfg.pad                     = 1;
cfg.method                  = 'mtmconvol';
cfg.taper                   = 'hanning';
cfg.foi                     = 4:4:100;                          
cfg.t_ftimwin               = 4./cfg.foi;  % 7 cycles per time window
cfg.toi                     = -0.1:0.01:0.3;                 
cfg.keeptrials              = 'yes';
gam_pow_recon               = ft_freqanalysis(cfg,gam_pow_data);

cfg                         = [];
cfg.baseline                = [-0.090 -0.01];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';
gam_pow_bsl                 = ft_freqbaseline(cfg, gam_pow_recon);

cfg                         = [];
cfg.parameter               = 'powspctrm';
cfg.zlim                    = 'maxmin';
% cfg.xlim                    = [-0.01 0.4]; 
% cfg.ylim                    = [60 100];    
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
figure,ft_singleplotTFR(cfg, gam_pow_bsl);
title('3924')
maxval
%%
cfg                         = [];
cfg.removemean              = 'yes';
EP_recon                    = ft_timelockanalysis(cfg, gam_pow_data);

%%
figure; hold on     % plot all objects in one figure
ft_plot_vol(vol_cm, 'facecolor', 'cortex', 'edgecolor', 'none'), alpha 0.4;
ft_plot_sens(proc.dataArtRej.grad, 'style', '*b');

% ft_plot_mesh(sourcespace), hold on
ft_plot_mesh(sourcemodel_lf.pos([564],:), 'vertexcolor', 'r');

%%
source.pos = sourcemodel_lf.pos;
source.tri = sourcemodel_lf.tri;  

template_mri = ft_read_mri('J:\old MNE files\014_NJO\mri\014_NJO.nii');
template_mri = ft_convert_units(template_mri, 'cm');

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'avg.pow';
cfg.interpmethod = 'nearest';
cfg.coordsys     = 'mni';
source_coh_int   = ft_sourceinterpolate(cfg, source, template_mri);

source_cm = ft_convert_units(source_coh_int, 'cm');

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'avg.pow';
cfg.coordsys     = 'mni';
ft_sourceplot(cfg, source_diff);

%%

% cfg = []
% proc.tlck = ft_timelockanalysis(cfg, proc.dataArtRej);

%%
% % fit a dipole to the M50 and M100 components
% cfg = [];
% cfg.latency = [0.020 0.060];  % specify latency window around M50 peak
% cfg.numdipoles = 1;
% % cfg.hdmfile = 'bauer_m.hdm';
% cfg.feedback = 'textbar';
% cfg.grid.resolution = 2;
% cfg.grid.unit = 'cm';
% cfg.vol             = vol_cm;
% cfg.grid.pos        = sourcemodel_lf;
% % cfg.grid.inside     = 1:size(cfg.grid.pos, 1);
% % cfg.grid.outside    = [];
% cfg.keepfilter      = 'yes';
% dipM50 = ft_dipolefitting(cfg, proc.tlck);


%%
cfg = [];
% note that it is actually better (i.e. more stable) to copy the atlas to 
% your folder and use cfg.atlas = 'TTatlas+tlrc.HEAD';
cfg.atlas='C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip/template/atlas/afni/TTatlas+tlrc.HEAD';
cfg.roi={'Brodmann area 17','Brodmann area 18'};
cfg.inputcoord = 'mni';
grid = ft_volumelookup(cfg, grid);






