clc;    close all;    clear all;

ft_defaults;
% Get data blanked
proc.data_folder                 = ...
    'j:\MEG_Research\SEF\SEF-TSSS\014_NJO' ;
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              = 10;
proc.save_folder                 = ...
    'j:\MEG_Research\SEF\SEF-TSSS\014_NJO';
% blank_ft2fif(data_folder, stim_chan, pre_stim_time, post_stim_time, save_folder);

%% Declarations

proc.blanked_folder         = proc.save_folder;                                  %'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO\blanked';
% cd(proc.blanked_folder)
filename                    = []; 
filename                    = dir('*.fif');
loop                        = 2;

% Import data, Define trials, Cut trials
proc.data_import            = [];
cfg                         = [];
cfg.dataset                 = filename(loop).name;
proc.data_import            = ft_preprocessing(cfg);

proc.data_trial = [];
cfg = [];
cfg.trialdef.prestim        = 0.1;
cfg.trialdef.poststim       = 0.35;
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
proc.data_epoched_EOG           = ft_redefinetrial(proc.data_trial,proc.preproc_data_EOG);

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
% cfg.artfctdef.zvalue.interactive = 'yes';
[artifact_jump] = ft_artifact_zvalue(cfg, proc.data_epoched);


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
[artefact_EOG_cfg, artifact_EOG] = ft_artifact_zvalue(cfg, proc.preproc_data_EOG );


cfg = [];
cfg.artfctdef.jump.artifact = artifact_jump.artfctdef.zvalue.artifact;
cfg.artfctdef.eog.artifact  = artifact_EOG;
cfg.artfctdef.reject    = 'complete';
% cfg.artfctdef.crittoilim = [-0.090 0.25];
% cfg.artfctdef.minaccepttim  = '0.1';
dataArtRej = ft_rejectartifact(cfg,proc.data_epoched);

%%
cfg            = [];
cfg.method = 'runica';
comp           = ft_componentanalysis(cfg, proc.preproc_data_MEG);

cfg = [];
cfg.component = [1:20];       % specify the component(s) that should be plotted
cfg.layout    = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
% cfg.channel = 'MEG***1';
ft_topoplotIC(cfg, comp);

cfg = [];
% cfg.channel = [1 5 9 11 12]; % components to be plotted
cfg.viewmode = 'component';
cfg.layout    = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
figure, ft_databrowser(cfg, comp)

proc.data_epoched           = [];
proc.data_epoched           = ft_redefinetrial(proc.data_trial,proc.preproc_data_MEG);

%% Timelock-Averaging

cfg                         = [];
cfg.removemean              = 'yes';
% cfg.trials                  = 5:length(proc.data_epoched.trial)-5;   
cfg.keeptrials              = 'yes';
proc.timelocked_data        = ft_timelockanalysis(cfg, proc.data_epoched);

cfg                         = [];
cfg.parameter               = 'trial';
cfg.baseline                = [-0.090 -0.005];
proc.timelocked_data_bsl   = ft_timelockbaseline(cfg, proc.timelocked_data);

% 
% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.trials                  = 10:length(proc.data_epoched.trial)-10;   
% cfg.channel                 =  'MEG***1'; 
% proc.timelocked_data_1      = ft_timelockanalysis(cfg, proc.data_epoched);
% 
% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.trials                  = 10:length(proc.data_epoched.trial)-10;   
% cfg.channel                 =  'MEG***2'; 
% proc.timelocked_data_2      = ft_timelockanalysis(cfg, proc.data_epoched);
% 
% cfg                         = [];
% cfg.removemean              = 'yes';
% cfg.trials                  = 10:length(proc.data_epoched.trial)-10;   
% cfg.channel                 =  'MEG***3'; 
% proc.timelocked_data_3      = ft_timelockanalysis(cfg, proc.data_epoched);
% 

cfg                         = [];
cfg.parameter               = 'avg';    
cmb_tlkdata                 = ft_combineplanar(cfg, proc.timelocked_data);

cfg                         = [];
cfg.showlabels              = 'no';
% cfg.parameter               = 'avg';    
cfg.layout                  = 'neuromag306cmb.lay';
% cfg.ylim                    = 'maxmin';
ft_multiplotER(cfg, cmb_tlkdata ); 

%% Moving average

cfg                         = [];
% cfg.demean                  = 'yes';
% cfg.baseline                = [-0.05 -0.005];
cmb_tlkdata                 = ft_combineplanar(cfg, proc.timelocked_data);

cfg                         = [];
cfg.showlabels              = 'no';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.ylim                    = 'maxmin'; %[-30e-13 100e-13];
% cfg.xlim                    = [-0.005 0.06];
ft_multiplotER(cfg, cmb_tlkdata); 

%% Moving average

cfg                         = [];
cfg.removemean              = 'yes';


tlck = [];
trl_loop = 1;

for loop = 10:15:length(proc.data_epoched.trial)-35
    disp(['**********  ', num2str(loop)]);
    cfg.trials                 = loop:loop+30;
    tlck.mov_ave(trl_loop)     = ft_timelockanalysis(cfg, proc.data_epoched);
    trl_loop = trl_loop+1;
end

% cfg                         = [];
% cfg.showlabels              = 'no';
% cfg.layout                  = 'neuromag306all.lay';
% cfg.channel                 =  'MEG***1'; 
% 
% cfg.ylim                    = [-3e-13 3e-13];
% % cfg.ylim                    = [-80e-13 80e-13];
% % cfg.ylim                    = [-80e-13 80e-13];
% cfg.xlim                    = [-0.005 0.06];
% ft_multiplotER(cfg, tlck.mov_ave(1)); 

%% Extract channels of interest

MEG_grad = {'MEG0413', 'MEG0443', 'MEG1813'};
MEG_mag  = {'MEG0221', 'MEG0231', 'MEG0411', 'MEG0441'};

for loop = 1:length(tlck.mov_ave)

    maps.avg_mag(loop,:) = ...
        mean(tlck.mov_ave(1,loop).avg(match_str(tlck.mov_ave(1).label, MEG_mag),:),1);
    
    maps.avg_grad(loop,:) = ...
        mean(tlck.mov_ave(1,loop).avg(match_str(tlck.mov_ave(1).label, MEG_grad),:),1);
    
end

figure, imagesc(tlck.mov_ave(1).time, 1:length(tlck.mov_ave) ,maps.avg_grad)
xlabel('Time in seconds'), ylabel('Epoch windows')
axis xy, title('Average of grad data')

figure, imagesc(tlck.mov_ave(1).time, 1:length(tlck.mov_ave) ,maps.avg_mag)
xlabel('Time in seconds'), ylabel('Epoch windows')
axis xy, title('Average of mag data')

% Peak amplitude, latency  

N20_time = find(tlck.mov_ave(1).time>0.021 & tlck.mov_ave(1).time<0.030);
N30_time = find(tlck.mov_ave(1).time>0.028 & tlck.mov_ave(1).time<0.036);

for loop = 1:size(maps.avg_grad,1)

    [N20_amp(loop) N20_lat(loop)] = findpeaks(maps.avg_grad(loop,N20_time));
    [N30_amp(loop) N30_lat(loop)] = findpeaks(-maps.avg_grad(loop,N30_time));
    
end

N20_lat = N20_lat + min(N20_time);
N30_lat = N30_lat + min(N30_time);

figure,
subplot 221, plot(tlck.mov_ave(1).time(N20_lat), '-o'),axis([1 30 0.02 0.03])
title('N20 Latency'),  ylabel('Latency in seconds')
subplot 222, plot(tlck.mov_ave(1).time(N30_lat), '-o'),axis([1 30 0.028 0.036])
title('N30 Latency'), 
subplot 223, plot(N20_amp, '-o')
title('N20 Amplitude'),  xlabel('Epochs'), ylabel('Amplitude a.u.')
subplot 224, plot(N30_amp, '-o')
title('N30 Amplitude'), xlabel('Epochs')


%%
cfg              = [];
cfg.paramter     = 'trial';
cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.pad          = 2;
cfg.taper        = 'hanning';
cfg.foi          = 1:5:100;                          
cfg.t_ftimwin    = 2./cfg.foi;   
cfg.toi          = -0.5:0.01:0.5;                 
TFRhann          = ft_freqanalysis(cfg, saveData.dataArtRej);

cfg = [];
cfg.baseline     = [-0.090 -0.01];
cfg.baselinetype = 'absolute';
cfg.parameter    = 'powspctrm';
TFRhann_bsl = ft_freqbaseline(cfg,TFRhann);

cfg                         = [];
TFRhann_bsl                 = ft_combineplanar(cfg, TFRhann_bsl);

cfg                         = [];
cfg.zlim                    = 'maxmin';
% cfg.baseline                = [-0.090 -0.01];
% cfg.baselinetype            = 'absolute';
% cfg.parameter               = 'powspctrm';
cfg.xlim                    = [0.01 0.34]; 
cfg.ylim                    = [1 100];    
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
ft_multiplotTFR(cfg, TFRhann_bsl);

%%
cfg              = [];
cfg.output       = 'pow';
cfg.paramter     = 'trial';
cfg.pad          = 1;
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials   = 'yes';
cfg.foi          = 30:5:100;                          
cfg.t_ftimwin    = 2./cfg.foi;  % 7 cycles per time window
cfg.toi          = -0.1:0.015:0.35;                 
TFRhann7 = ft_freqanalysis(cfg,proc.data_epoched);

cfg.baseline     = [-0.090 -0.001];
cfg.baselinetype = 'absolute';
cfg.parameter    = 'powspctrm';
TFRhann_bsl_7 = ft_freqbaseline(cfg, TFRhann7);

cfg                         = [];
TFRhann_abs_7                 = ft_combineplanar(cfg, TFRhann_bsl_7);

cfg                         = [];
cfg.colorbar                = 'yes';
cfg.maskstyle               = 'saturation';	
cfg.layout                  = 'neuromag306all.lay';
cfg.zlim = [-4e-26 2e-26 ]
ft_multiplotTFR(cfg, TFRhann_abs_7);

%%
cfg                         = [];
plot_1                      = ft_combineplanar(cfg, proc.HGamma.tl_bsl);

cfg                         = [];
cfg.showlabels              = 'no';
cfg.parameter               = 'avg';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.ylim                    = 'maxmin';
cfg.xlim                    = [-0.09 0.4];
ft_multiplotER(cfg, plot_1); 

%%
freq_chg = freq_bsl;
freq_chg.powspctrm = (freq_exp.powspctrm-freq_bsl.powspctrm);

cfg = []
cfg.parameter = 'powspctrm';
cfg.layout                  = 'neuromag306cmb.lay';
figure,ft_topoplotTFR(cfg, ft_combineplanar([],freq_chg))


