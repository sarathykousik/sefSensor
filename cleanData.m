%% cleanData.m
%  Import, epoch, filter, artefact data
%  Save Data - patient-wise

cd(proc.data_folder);
filename                    = dir('*.fif');
loop                        = 2;

%% Import data, Define trials, Cut trials
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

%% Filtering, Baseline correction and extraction of MEG channels only
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
