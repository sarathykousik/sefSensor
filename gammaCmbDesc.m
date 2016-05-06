clear all;   clc;    close all;
 
%
PD_folder = 'J:\MEG_Research\SEF\SEFGamma';
Control_folder = 'J:\MEG_Research\SEF\SEFGammaControl';

%% PD
cd(PD_folder)
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 7, []);
file_sub = file_sub';
[row col] = size(file_sub);

for subj = 1:row
    disp(' ')
    disp(' ')
    disp(['######   ', num2str(subj)])
    disp(' ')
    disp(' ')
    
    for cond = 1:col

        tok = tokenize(char(file_sub(subj, cond)),'-');
        load(char(file_sub(subj, cond)));
        cfg = [];
        cfg.baseline                = [-0.15 -0.1];
        cfg.baselinetype            = 'relchange';
        cfg.parameter               = 'powspctrm';
        gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);
       
        cfg                         = [];
        PDgammaCmb{subj,cond}       = ft_combineplanar(cfg, gammaFreq_bsl);
        PDgammaCmb{subj,cond}.name  = tok{1};
        
%         cfg                             = [];
%         cfg.jackknife                   = 'yes';
%         PDgammaDesc{subj,cond}          = ft_freqdescriptives(cfg, PDgammaCmb{subj,cond});
%         PDgammaDesc{subj,cond}.name     = tok{1};

    end
end

%% Control
cd(Control_folder)
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 7, []);
file_sub = file_sub';
[row col] = size(file_sub);

for subj = 1:row
    disp(' ')
    disp(' ')
    disp(['######   ', num2str(subj)])
    disp(' ')
    disp(' ')
    
    for cond = 1:col

        tok = tokenize(char(file_sub(subj, cond)),'-');
        load(char(file_sub(subj, cond)));
        cfg = [];
        cfg.baseline                = [-0.15 -0.1];
        cfg.baselinetype            = 'relchange';
        cfg.parameter               = 'powspctrm';
        gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);
       
        cfg                         = [];
        ControlgammaCmb{subj,cond}  = ft_combineplanar(cfg, gammaFreq_bsl);
        ControlgammaCmb{subj,cond}.name  = tok{1};
        
%         cfg                             = [];
%         cfg.jackknife                   = 'yes';
%         ControlgammaDesc{subj,cond}     = ft_freqdescriptives(cfg, ControlgammaCmb{subj,cond});
%         ControlgammaDesc{subj,cond}.name     = tok{1};
%         
    end
end

%%
cfg             = [];
gammaDesc       = ft_freqdescriptives(cfg, ControlgammaCmb{1,1});

%%
for subj = 1:10
    for cond = 1:7
        
        cfg_toi               = [];
        cfg_toi.latency       = [0.02];
        topoGammaControl{subj,cond}  = ft_selectdata(cfg_toi, ControlgammaCmb{subj,cond});
        topoGammaControl{subj,cond}.powspctrm = ...
            mean(max(ControlgammaCmb{subj,cond}.powspctrm(:,:,:,find(ControlgammaCmb{subj,cond}.time>=0.02 & ControlgammaCmb{subj,cond}.time<=0.09)),[],4))';       

        topoGammaPD{subj,cond}  = ft_selectdata(cfg_toi, PDgammaCmb{subj,cond});
        topoGammaPD{subj,cond}.powspctrm = ...
            mean(max(PDgammaCmb{subj,cond}.powspctrm(:,:,:,find(PDgammaCmb{subj,cond}.time>=0.02 & PDgammaCmb{subj,cond}.time<=0.09)),[],4))';       

        
    end
end

%%

for subj = 1:10
    for cond = 1:7
   
        topoGammaPD{subj,cond}.dimord = 'chan_freq';
%         topoGammaPD{subj,cond}.powspctrm = topoGammaPD{subj,cond}.powspctrm';

        
        topoGammaControl{subj,cond}.dimord = 'chan_freq';
%         topoGammaC{subj,cond}.powspctrm = topoGammaC{subj,cond}.powspctrm';
        
    end
end

%%
cfg = [];
cfg.parameter = 'powspctrm';
for cond = 1:7
    
       GAgammaPD{cond}      = ft_freqgrandaverage(cfg, topoGammaPD{:,cond});
       GAgammaControl{cond} = ft_freqgrandaverage(cfg, topoGammaControl{:,cond});
        
end



