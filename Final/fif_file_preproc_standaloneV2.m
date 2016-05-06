% function output = fif_file_preproc(file_name,subj, condN ,proc)
clc;    clear all;     close all;
ft_defaults;
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')
proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\controlNew\xscan30';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -10;
proc.post_stim_time              =  8;
proc.save_folder                 = proc.data_folder;% 'J:\MEG_Research\ArtRejtSSS\results';

cd(proc.data_folder)
filenames      = dir('*.fif');
mkdir(proc.save_folder)

%%

for condLoop  = 2% 1:length(filenames)
 %%   
    [a b c]     =   fileparts(filenames(condLoop).name);

    disp('#######################################')
    disp(['****    ', b ,'        ****'])
    disp('#######################################')
    
    % Import file
    cfg                         = [];
    cfg.dataset                 = filenames(condLoop).name;%filenames(condLoop).name;
    proc.data_import            = ft_preprocessing(cfg);
    
    % Find epochs
    cfg = [];
    cfg.trialdef.prestim        = 0.5;
    cfg.trialdef.poststim       = 0.5;
    cfg.trialdef.eventtype      = proc.stim_chan ; 
    cfg.trialdef.eventvalue     = 5;
    cfg.dataset                 = filenames(condLoop).name;%filenames(condLoop).name;
    proc.data_trial             = ft_definetrial(cfg);

    % Epoch data
%     time_tot    = floor(proc.data_import.time{1}(end));
%     fs          = proc.data_import.fsample;
%     sample_tot  = time_tot*fs;
%     trl         = [];
%     trl         = [[0:fs:sample_tot-fs]' [fs:fs:sample_tot]'];
%     trl(:,3)    = 0;%[fs.*size(trl,1)]';
%     trl(:,1)    = trl(:,1)+1;


    % Blanking
%     proc.blank_seq              = ones(size(trl,2),1);
    proc.blank_seq              = ones(size(proc.data_import.trial{1},2),1);
    for ind = proc.pre_stim_time:proc.post_stim_time
        proc.blank_seq(proc.data_trial.trl(:,1)+(cfg.trialdef.prestim*1000)+ind) = 0;
%          proc.blank_seq(trl(:,1)+ind) = 0;
    end

    proc.MEG_channel            = find(strncmp(proc.data_import.label, 'MEG',3)==1);                    % Pick only MEG channels
    proc.data_import.trial{1}(proc.MEG_channel,:) = ...
            bsxfun(@times, proc.data_import.trial{1}(proc.MEG_channel,:), proc.blank_seq');  % Multiply blanking sequence with MEG channels

    % Filtering, Baseline correction and extraction of MEG channels only
    cfg                         = [];
    cfg.lpfilter                = 'yes';
    cfg.lpfreq                  =  100;
    cfg.hpfilter                = 'yes';
    cfg.hpfreq                  =  1;
    cfg.dftfilter               = 'yes';
    cfg.channel                 =  {'MEG'}; 
    cfg.detrend                 = 'yes';
    proc.preproc_data_MEG       = ft_preprocessing(cfg, proc.data_import);
    
    cfg = [];
    cfg.trl = proc.data_trial.trl;
    proc.data_epoched = ft_redefinetrial(cfg, proc.preproc_data_MEG);
    
%%
    cfg.channel                 =  {'EOG*'}; 
    proc.preproc_data_EOG       = ft_preprocessing(cfg, proc.data_import);
    cfg.channel                 =  {'ECG*'}; 
    proc.preproc_data_ECG       = ft_preprocessing(cfg, proc.data_import);

    proc.data_epoched           = ft_redefinetrial(proc.data_trial,proc.preproc_data_MEG);
    proc.data_epoched_EOG       = ft_redefinetrial(proc.data_trial,proc.preproc_data_EOG);
    proc.data_epoched_ECG       = ft_redefinetrial(proc.data_trial,proc.preproc_data_ECG);
    proc.data_epoched           = ft_appenddata([],...
            proc.data_epoched, proc.data_epoched_EOG,proc.data_epoched_ECG);
    %
    %Jump artefact
    cfg                                 = [];
    cfg.trl                             = proc.data_trial.trl;
    cfg.continuous                      = 'no';
    % cutoff and padding
    cfg.artfctdef.zvalue.channel        = 'MEG*';
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
    cfg.artfctdef.zvalue.cutoff         = 6;
    cfg.artfctdef.zvalue.trlpadding     = 0;
    cfg.artfctdef.zvalue.artpadding     = 0.1;
    cfg.artfctdef.zvalue.fltpadding     = 0;
    % algorithmic parameters
    cfg.artfctdef.zvalue.bpfilter       = 'yes';
    cfg.artfctdef.zvalue.bpfilttype     = 'but';
    cfg.artfctdef.zvalue.bpfreq         = [1 15];
    cfg.artfctdef.zvalue.bpfiltord      = 4;
    cfg.artfctdef.zvalue.hilbert        = 'yes';
    [~, artifact_EOG]                   = ft_artifact_zvalue(cfg, proc.data_epoched_EOG);
    disp(['##########################################'])
    cfg                                 = [];
    cfg.artfctdef.jump.artifact         = artifact_jump;
    cfg.artfctdef.eog.artifact          = artifact_EOG;
    cfg.artfctdef.reject                = 'complete';
    cfg.artfctdef.crittoilim            = [-0.100 0.300];
    dataArtRej                          = ft_rejectartifact(cfg,proc.data_epoched);
    disp(['##########################################'])
    tok_name = tokenize(b, '_');
%     dataArtRej=proc.data_epoched;
    
    saveData.file                       = filenames(condLoop).name;
%     saveData.subject                    = [tok_name{1},'_',tok_name{1}]; 
%     saveData.condName                   = condLoop; 
    saveData.dataArtRej                 = dataArtRej ;
    saveData.origTrialLen               = length(dataArtRej.trial);
    
    % Saving 
%     save ([proc.save_folder,'\', subj,'\',b, '-SaveData'], 'saveData', '-v7.3');
%     disp('*********************************')
%     disp(['Saved data to:   ', [proc.save_folder,'\', subj,'\',b, '-SaveData']]);
%     disp('*********************************')

%     save ([proc.save_folder,'\',b, '-SaveData'], 'saveData', '-v7.3');
    disp('*********************************')
    disp(['Saved data to:   ', [proc.save_folder,'\saveData\',b, '-SaveData']]);
    disp('*********************************')

    % Calc Ind gamma
    gammaFreq = calcGamma(saveData.dataArtRej); %calcGamma(visClean);%
    save ([proc.save_folder, '\',b, '-gammaFreq'], 'gammaFreq');

    proc.data_import = [];
    proc.data_trial= [];
    proc.preproc_data_MEG= [];
    proc.preproc_data_EOG= [];
    proc.data_epoched    = [];
    proc.data_epoched_EOG= [];
    artifact_jump= [];
    artifact_EOG= [];
    dataArtRej= [];
    saveData= [];
end

%%

plot(proc.data_import.time{1}, proc.blank_seq.*5e-13, 'k*')
hold on,
plot((proc.data_trial.trl(:,1)+(.5*1000))./1000, ...
    ones(1,length(proc.data_trial.trl(:,1))).*10e-13, 'bo')
plot(proc.preproc_data_MEG.time{1}, proc.preproc_data_MEG.trial{1}(15,:),'r')

