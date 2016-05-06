
clc;    clear all;      close all;

%%
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 6, []);
file_sub = file_sub';
[row col] = size(file_sub);

for subj = 1:row
    disp(' ')
    disp(' ')
    disp(['######   ', num2str(subj)])
    disp(' ')
    disp(' ')
    
    for cond = 1:col
        data = load(char(file_sub(subj, cond)));
        
        
        cfg = [];
        cfg.baseline                = [-0.1 -0.04];
        cfg.baselinetype            = 'absolute';
        cfg.parameter               = 'powspctrm';
        ITPC_bsl                    = ft_freqbaseline(cfg,data.TFRfourier);
     
        cfg                                   = [];
        cfg.method                            = 'itc' ;
        Control_ITPC_mtm{subj, cond}          = ft_combineplanar_itc(cfg, ITPC_bsl);
        Control_ITPC_mtm{subj, cond}.filename = char(file_sub(subj, cond));
        
    end
end

%%
cfg_toi               = [];
cfg_toi.latency       = [0.02 0.06];
 
% cfg_sing               = [];
% cfg_sing.latency       = [0.02];

for subj = 1:10
    for cond = 1:6
        
        topoITPC{subj,cond}         = ft_selectdata(cfg_toi, Control_ITPC_mtm{subj,cond});
%         ITPC{subj,cond}             = ft_selectdata(cfg_sing, Control_ITPC_mtm{subj,cond});
        maxITPCTopo(subj,cond,:)    = mean(topoITPC{subj,cond}.powspctrm,3);
%         ITPC{subj,cond}.powspctrm   = maxITPCTopo;

    end
end

% load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20.mat')
load('J:\MEG_Research\SEF\SEFepresults\control\maxPosN20Control.mat')
load('J:\MEG_Research\SEF\neighbours.mat')

for subj = 1:10
    for cond = 1:6
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(topoITPC{subj,cond}.label, chan_sel);
        ControlITPC_mtm(subj, cond) = mean(maxITPCTopo(subj,cond,chanPos),3);
    end
end

%%
contS = [1:8 9 10];
subjP = [1  3:10];

[p, table] = anova_rm({ControlITPC_mtm(contS,:), PDITPC_mtm(subjP,:)});


  
