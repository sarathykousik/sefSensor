function output = fif_file_preproc(file_name,subj, condN ,proc)
    
    [a b c]     =   fileparts(file_name);
    
    disp('#######################################')
    disp(['****    ', b ,'        ****'])
    disp('#######################################')
    
    % Import file
    cfg                         = [];
    cfg.dataset                 = file_name;
    proc.data_import            = ft_preprocessing(cfg);
    
    % Find epochs
    cfg = [];
    cfg.trialdef.prestim        = 0.5;
    cfg.trialdef.poststim       = 0.5;
    cfg.trialdef.eventtype      = proc.stim_chan ; 
    cfg.trialdef.eventvalue     = 5;
    cfg.dataset                 = file_name;
    proc.data_trial             = ft_definetrial(cfg);
    
    % Blanking
    proc.blank_seq              = ones(size(proc.data_import.trial{1},2),1);
    for ind = proc.pre_stim_time:proc.post_stim_time
        proc.blank_seq(proc.data_trial.trl(:,1)+(cfg.trialdef.prestim*1000)+ind) = 0;
    end

    proc.MEG_channel            = find(strncmp(proc.data_import.label, 'MEG',3)==1);                    % Pick only MEG channels
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
    cfg                                 = [];
    cfg.trl                             = proc.data_trial.trl;
    cfg.continuous                      = 'no';
    % cutoff and padding
    cfg.artfctdef.zvalue.channel        = 'MEG';
    cfg.artfctdef.zvalue.cutoff         = 30;
    cfg.artfctdef.zvalue.trlpadding     = 0;
    cfg.artfctdef.zvalue.artpadding     = 0;
    cfg.artfctdef.zvalue.fltpadding     = 0;
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative     = 'yes';
    cfg.artfctdef.zvalue.medianfilter   = 'yes';
    cfg.artfctdef.zvalue.medianfiltord  = 9;
    cfg.artfctdef.zvalue.absdiff        = 'yes';
    [~, artifact_jump]                  = ft_artifact_zvalue(cfg, proc.data_epoched);

    % EOG
    cfg                                 = [];
    cfg.trl                             = proc.data_trial.trl;
    cfg.continuous                      = 'no'; 
    cfg.trl                             = proc.data_trial.trl;
    % cutoff and padding
    cfg.artfctdef.zvalue.channel        = 'EOG*';
    cfg.artfctdef.zvalue.cutoff         = 100;
    cfg.artfctdef.zvalue.trlpadding     = 0;
    cfg.artfctdef.zvalue.artpadding     = 0.1;
    cfg.artfctdef.zvalue.fltpadding     = 0;
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter       = 'yes';
    cfg.artfctdef.zvalue.bpfilttype     = 'but';
    cfg.artfctdef.zvalue.bpfreq         = [1 15];
    cfg.artfctdef.zvalue.bpfiltord      = 6;
    cfg.artfctdef.zvalue.hilbert        = 'yes';
    [~, artifact_EOG]                   = ft_artifact_zvalue(cfg, proc.data_epoched_EOG);

    cfg                                 = [];
    cfg.artfctdef.jump.artifact         = artifact_jump;
    cfg.artfctdef.eog.artifact          = artifact_EOG;
    cfg.artfctdef.reject                = 'complete';
    cfg.artfctdef.crittoilim            = [-0.100 0.300];
    dataArtRej                          = ft_rejectartifact(cfg,proc.data_epoched);

    saveData.file                       = file_name;
    saveData.subject                    = subj; 
    saveData.condName                   = condN; 
    saveData.dataArtRej                 = dataArtRej;
    saveData.origTrialLen               = length(proc.data_epoched.trial);
    
    % Saving 
%     save ([proc.save_folder,'\', subj,'\',b, '-SaveData'], 'saveData', '-v7.3');
%     disp('*********************************')
%     disp(['Saved data to:   ', [proc.save_folder,'\', subj,'\',b, '-SaveData']]);
%     disp('*********************************')

    save ([proc.save_folder,'\',b, '-SaveData'], 'saveData', '-v7.3');
    disp('*********************************')
    disp(['Saved data to:   ', [proc.save_folder,'\', subj,'\',b, '-SaveData']]);
    disp('*********************************')
return