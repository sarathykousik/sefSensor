clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFcleanData';
proc.save_folder                 = 'J:\MEG_Research\SEF\hilbERD-v2';
mkdir(proc.save_folder)

%%
cd(proc.data_folder)
filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)

    cfg = []; 
    cfg.hilbert = 'abs';
    cfg.bpfilter  = 'yes';
    cfg.bpfreq   = [40 90];    
    cfg.bpfiltdir = 'twopass';
    cfg.bpfilttype = 'but';
    cfg.bpfiltord     = 5;
    % cfg.trials = 1:100;
    hilbData = ft_preprocessing(cfg, saveData.dataArtRej);

    %
    cfg                     = [];
    hilbComb                = ft_combineplanar(cfg, hilbData);

    cfg = [];
    tlckHilb = ft_timelockanalysis(cfg, hilbComb);

    cfg = [];
    cfg.baseline = [-0.09 -0.01];
    cfg.baselinetype = 'absolute';
    % cfg.parameter = 'trial';
    hilbBase = ft_timelockbaseline(cfg, tlckHilb);

    save ([proc.save_folder,'\',b, '-hilbBase'], 'hilbBase', '-v7.3');
end

%% Result matrix - loading and saving
cd(proc.save_folder)
filenames_ERD      = dir('*.mat');

for loop = 1: length(filenames_ERD)
    
    file_sub(loop) = {filenames_ERD(loop).name};

end

file_sub = reshape(file_sub, 6, [])';
[row col] = size(file_sub);

hilbERD=[];
for subjLoop = 1:row

    disp(['######   ', num2str(subjLoop)])
    for condLoop = 1:col
        disp(['******** ', char(file_sub(subjLoop, condLoop))])
        data = load(char(file_sub(subjLoop, condLoop)));
        hilbERD{subjLoop, condLoop} =  data.hilbBase;
        hilbERD{subjLoop, condLoop}.name = char(file_sub(subjLoop, condLoop));
    end
end

hilbERD = hilbERD';         
%%

toi = find(hilbERD{1}.hilbBase.time>0.02 & hilbERD{1}.hilbBase.time<0.061);

% for subjLoop = 1:row
% 
%     disp(['######   ', num2str(subjLoop)])
%     for condLoop = 1:col
%         hilbERD{subjLoop, condLoop}.hilbBase. 
%         
%         
%     end
% end

%%
% for subj = 1:6
subj = 4;    
cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
%     cfg.channel  = {'MEG0432+0433', 'MEG0442+0443'};
cfg.xlim = [-0.1 0.1];
cfg.ylim = 'maxmin';
h1 = figure, 
ft_multiplotER(cfg, hilbERD{1,subj},hilbERD{2,subj},hilbERD{3,subj},...
                hilbERD{5,subj},hilbERD{subj,6})
title(num2str(subj));
%     waitfor(h1)
% end

%% topoplots
for subj= 1:6

    for cond = 1:6
        cfg = [];
        cfg.parameter = 'avg';
        cfg.layout    = 'neuromag306cmb.lay';
    %     cfg.channel  = {'MEG0432+0433', 'MEG0442+0443'};
        cfg.xlim = [0.02 0.06];
        cfg.zlim = [-4.43e-13 14.7e-13];
        cfg.colorbar = 'yes';
        figure
        ft_topoplotER(cfg, hilbERD{subj,cond})
        title(['Subj: ' num2str(subj), '  Cond: ' num2str(cond)']);
        saveas(gcf,['J:\MEG_Research\SEF\hilbTopo','\',num2str(subj),'-',num2str(cond),'-gammaHilb'],'jpg')
    end
end

%% Statistics


cfg             = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.template    = 'neuromag306cmb_neighb.mat';               % specify type of template
cfg.layout      = 'neuromag306cmb.lay';
cfg.feedback    = 'yes';                             % show a neighbour plot 
neighbours      = ft_prepare_neighbours(cfg, hilbERD{1,1}); % define neighbouring channels

%% get act and bsl

cfg_bsl = [];
cfg_bsl.latency = [-0.09 -0.05];

cfg_act = [];
cfg_act.latency = [0.02 0.06];

for condLoop = 1:6
    
hilbERD_act{condLoop} = ft_selectdata(cfg_act,hilbERD{condLoop,6});
hilbERD_bsl{condLoop} = ft_selectdata(cfg_bsl,hilbERD{condLoop,6});

end


%
cfg = [];
cfg.channel     = 'MEG';
% cfg.latency     = [0.02 0.06];
cfg.avgovertime = 'yes';

cfg.method = 'montecarlo';
cfg.statistic = 'actvsblT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.01;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.tail = 1;
cfg.clustertail = 1;
cfg.alpha = 0.01;
cfg.numrandomization = 500;

Nsub = 6;
Ncond = 2;
cfg.design(1,1:Ncond*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];% 3*ones(1,Nsub) ...
%                         4*ones(1,Nsub) 5*ones(1,Nsub) 6*ones(1,Nsub)];
cfg.design(2,1:Ncond*Nsub)  = [1:Nsub 1:Nsub ];%1:Nsub 1:Nsub 1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,  hilbERD_act{:},hilbERD_bsl{:});

% stat_pval = [stat.posclusters(:).prob];
% plot(stat.mask);figure(gcf);


cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.highlight = 'labels';
cfg.highlightchannel = stat.label(stat.mask==1);
ft_topoplotER(cfg, hilbERD_bsl{1,1})

%% Extract only 20-60 ms - mean of all trials and max - topos
% t-test on topos
% close all
subj = 6;
for cond = 1:6

cfg               = [];
cfg.latency       = [0.02 0.06];
topoGamma         = ft_selectdata(cfg, hilbERD{cond,subj});

maxGamma_ch = mean(topoGamma.avg,2);
% 
% [h p ci stats] = ttest(maxG_Z,mean(maxG_Z),0.05, 'right')


cfg               = [];
cfg.latency       = [0.02];
topoGamma_thresh  = ft_selectdata(cfg, hilbERD{cond,subj});
topoGamma_thresh.avg = maxGamma_ch.* (maxGamma_ch>(mean(maxGamma_ch)+2*std(maxGamma_ch))); 
% figure,plot(topoGamma_thresh.avg)

cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.colorbar = 'yes';
figure(cond),
subplot 121, ft_topoplotER(cfg, topoGamma_thresh), title(['Thresh', num2str(cond)])
subplot 122, ft_topoplotER(cfg, topoGamma), title('True')
suptitle(num2str(cond))
end

