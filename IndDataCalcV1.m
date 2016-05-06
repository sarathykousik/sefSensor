clc;        clear all;      close all;

% load('j:\MEG_Research\SEF\SEFepresults\updated\tlck.mat');

%%
% cd('j:\MEG_Research\SEF\SEFVisClean')
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 6, []);
file_sub = file_sub';
[row col] = size(file_sub);

%% Induced data

% for subj = 1:row   %subj
subj=11
    disp(['######  subj ', num2str(subj)])

for cond = 1:col   %col
        [a b c]     =   fileparts(char(file_sub(subj, cond)));
        disp(['   **** visClean   ',b,'     ****    '])
        disp(['   #### visClean   ',tlck{subj,cond}.name,'     ####'])
        load(char(file_sub(subj, cond)));
        visCleanGrad = [];  ind_dat = [];

        cfg= []; cfg.channel = 'MEGGRAD';
        visCleanGrad = ft_selectdata(cfg,visClean);

        for trialLoop = 1: length(visCleanGrad.trial)
            ind_dat{1,trialLoop} = ...
                visCleanGrad.trial{1,trialLoop} - tlck{subj,cond}.avg;
        end

        visCleanInd = visCleanGrad;
        visCleanInd.trial = ind_dat;
        visCleanInd.filename =char(file_sub(subj, cond));
        
       save (['J:\MEG_Research\SEF\SEFCleanInd\',b, '-visCleanInd'], 'visCleanInd', '-v7.3');

end
% end

%% AVG of neigh per subj/cond; TFR plot; 

cfg_plot                     = [];
cfg_plot.parameter           = 'powspctrm';
cfg_plot.colorbar            = 'yes';
cfg_plot.layout              = 'neuromag306cmb.lay';
cfg_plot.xlim                = [0 0.2];
cfg_plot.shading             = 'interp';
cfg_plot.maskstyle           = 'saturation';	

cfg_sel                      = [];
cfg_sel.parameter            = 'powspctrm';
% cfg_sel.avgoverchan          = 'yes';

for subj = 1:row
    for cond = 1:col
        
        % load data
        [a b c]     =   fileparts(char(file_sub(subj, cond)));
        disp(['   ****    ',b,'     ****    '])
        load(char(file_sub(subj, cond)));  
        
        cfg_desc = [];
        cfg_desc.jackknife   = 'yes';
        TFRwave_desc = ft_freqdescriptives(cfg_desc, TFRwave_bsl);
   
        
        % Select neigh data
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                    neighbours(1,maxPosN20(subj,cond)).label};
        cfg_sel.channel =  chan_sel;
        SEF_IND_sel{subj,cond} = ft_selectdata(cfg_sel, ft_combineplanar([],TFRwave_desc));
        SEF_IND_sel{subj,cond}.label = {'meanNeigh'};
        SEF_IND_sel{subj,cond}.filename   = b;
        
        
        % TFR plot - mean of neighbours
        hFig = figure(subj);
        set(hFig, 'Position', [20 60 1824 1050]);
        subplot(2,3,cond),
        ft_singleplotTFR(cfg_plot, SEF_IND_sel{subj,cond})
        title([SEF_IND_sel{subj,cond}.filename(1:3), '-', num2str(cond)])
        
        
    end
    
    set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...      % figure handle
    ['single gradITC stat Cluster'],... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', ...           % file format
    '-r250' );             % resolution in dpi

end

%% Grand average TFR

cfg = [];
cfg.parameter = 'powspctrm';

cfg_ind = [];
cfg_ind.parameter = 'powspctrm';
cfg_ind.keepindividual = 'yes';

for cond = 1:6
   
    SEF_TFR_cond_avg{cond} = ft_freqgrandaverage(cfg, SEF_TFR_sel{:,cond});
    SEF_TFR_cond_ind{cond} = ft_freqgrandaverage(cfg_ind, SEF_TFR_sel{:,cond});
    
end


%% TFR condition-wise

cfg_plot                     = [];
cfg_plot.parameter           = 'powspctrm';
cfg_plot.colorbar            = 'yes';
cfg_plot.layout              = 'neuromag306cmb.lay';
cfg_plot.xlim                = [0.01 0.2];
cfg_plot.ylim                = [10 100];
cfg_plot.channel             = {'all', '-MEG1722+1723', '-MEG1812+1813'};
% cfg_plot.zlim                = [0 12];
cfg_plot.shading             = 'interp';
cfg_plot.maskstyle           = 'saturation';

% for cond = 1:col
%     hFig = figure(cond);
%     set(hFig, 'Position', [20 60 1824 1050]);
    ft_multiplotTFR(cfg_plot, TFRInd_cmb{1,2})
%     title(num2str(cond))
%     
%       set(gcf, 'Color', 'white'); % white bckgr
%     export_fig( gcf, ...      % figure handle
%     ['UPDRS vsGamma-Med ON-StimON'],... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', ...           % file format
%     '-r250' );             % resolution in dpi
% end        

%%
for cond = 1:6
    cfg                 = [];
    cfg.parameter       = 'powspctrm';
    cfg.operation       = 'log10';
    TFR_logpow{cond}   = ft_math(cfg, SEF_ITC_cond{1,cond});
end

%%
chan_left = {'MEG0112+0113', 'MEG0122+0123', 'MEG0132+0133',...
    'MEG0142+0143', 'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233',...
    'MEG0242+0243', 'MEG0312+0313', 'MEG0322+0323', 'MEG0332+0333', ...
    'MEG0342+0343', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433',...
    'MEG0442+0443', 'MEG0512+0513', 'MEG0522+0523', 'MEG0532+0533', ...
    'MEG0542+0543', 'MEG0612+0613', 'MEG0622+0623', 'MEG0632+0633', ...
    'MEG0642+0643', 'MEG0712+0713', 'MEG0742+0743', 'MEG0812+0813', ...
    'MEG0822+0823', 'MEG1012+1013', 'MEG1512+1513', 'MEG1522+1523',...
    'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623',...
    'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723',...
    'MEG1732+1733', 'MEG1742+1743', 'MEG1812+1813', 'MEG1822+1823', ...
    'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923',...
    'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2042+2043',...
    'MEG2112+2113', 'MEG2122+2123', 'MEG2142+2143'};

%%
load('J:\MEG_Research\SEF\neighbours.mat')

%%

ntrials = 120;
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];



%% Statistics over grand average
% combos = nchoosek([1:6],2);

% for loop = 1:length(combos)
%     loop
    cfg                     = [];
    cfg.latency             = [0 0.1];
    cfg.frequency           = [10 100];
    cfg.parameter           = 'powspctrm';
    % cfg.channel             = {'all', '-MEG1612+1613', '-MEG2642+2643'};
    % cfg.channel = {'MEG0232+0233', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433',...
    %     'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
    %     'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
    cfg.statistic           = 'ft_statfun_depsamplesT';
%     cfg.statistic           = 'ft_statfun_depsamplesFunivariate';
    % cfg.statistic           = 'ft_statfun_depsamplesregrT';
    % cfg.method              = 'analytic';
    % cfg.correctm            = 'no';
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    cfg.clusterstatistic    = 'maxsum';
    cfg.minnbchan           = 2;
    cfg.tail                = 0;
    cfg.clustertail         = 0;
    cfg.alpha               = 0.05;
    cfg.clusteralpha        = 0.05;
    cfg.numrandomization    = 500;
    % cfg_neighb.method       = 'distance';
    % cfg.neighbours          = ft_prepare_neighbours(cfg_neighb, SEFITC_cmb{1,1}.grad);
    cfg.neighbours = neighbours;

    design                  = [];
    subj                    = [1 3:10];
    nsubj                   = length(subj);
    design(1,:)             = [ones(1,nsubj), 2*ones(1,nsubj)];% 3*ones(1,nsubj)...
%                                 4*ones(1,nsubj) 5*ones(1,nsubj) 6*ones(1,nsubj)];
    design(2,:)             = [1:nsubj 1:nsubj];% 1:nsubj 1:nsubj 1:nsubj 1:nsubj ];
    cfg.design              = design;
    cfg.ivar                = 1;
    cfg.uvar                = 2;
    stat_log         = ft_freqstatistics(cfg,ITC_bsl{subj,1},ITC_bsl{subj,6});
%                                 SEFITC_cmb{subj,combos(loop,1)}, ...
%                                     SEFITC_cmb{subj,combos(loop,2)});
%     stat_log{loop}.cond          = [combos(loop,:)];

%%

% temp_stat = stat_log;
% temp_stat.stat  = stat_log.stat.*stat_log.stat>18;

cfg                     = [];
cfg.xlim                = [0.02 0.03];
cfg.ylim                = [60 100];
cfg.parameter           = 'stat';
cfg.maskparameter       = 'mask';
cfg.maskstyle           = 'saturation';
cfg.shading             = 'interp';
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
figure
ft_topoplotTFR(cfg, stat_log);


%%


cfgact                     = [];
cfgact.parameter           = 'powspctrm';
% cfgact.avgoverchan         = 'yes';
cfgact.foilim              = [60 100];
% cfgact.latency             = [0.02 0.08];

cfgbsl                     = [];
cfgbsl.parameter           = 'powspctrm';
% cfgbsl.avgoverchan         = 'yes';
cfgbsl.foilim              = [60 100];
cfgbsl.latency             = [-0.08 -0.02];

% for subj = 1:120
%     disp(['###### ',num2str(subj)])
%     for cond = 1:6
%         disp(['###### ',num2str(cond)])
% 
%         SEF_ITC_act{subj} = ft_selectdata(cfgact, SEF{subj});
%         SEF_ITC_bsl{subj} = ft_selectdata(cfgbsl, SEF{subj});

SEF_ITC_act = ft_selectdata(cfgact, ft_combineplanar([],TFRwave));
SEF_ITC_bsl = ft_selectdata(cfgbsl, ft_combineplanar([],TFRwave));
SEF_ITC_bsl.time = SEF_ITC_act.time;
        
%     end
% end

cfg                     = [];
cfg.parameter           = 'powspctrm';
cfg.statistic           = 'ft_statfun_actvsblT';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 2;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.numrandomization    = 500;
cfg_neighb.method       = 'distance';
%     cfg.neighbours          = ft_prepare_neighbours(cfg_neighb, SEFITC_cmb{1,1}.grad);
cfg.neighbours          = neighbours;

design                  = [];
%     subj                    = [1:120];
nsubj                   = length(TFRwave.trialinfo);
design(1,:)             = [ones(1,nsubj), 2*ones(1,nsubj)];% 3*ones(1,nsubj)...
%                                 4*ones(1,nsubj) 5*ones(1,nsubj) 6*ones(1,nsubj)];
design(2,:)             = [1:nsubj 1:nsubj];% 1:nsubj 1:nsubj 1:nsubj 1:nsubj ];
cfg.design              = design;
cfg.ivar                = 1;
cfg.uvar                = 2;
stat_log                = ft_freqstatistics(cfg,SEF_ITC_act, SEF_ITC_bsl);
    
%%

% Then take the difference of the averages using ft_math
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
raweffectFICvsFC = ft_math(cfg,ft_freqgrandaverage([],SEFITC_cmb{:,1}),...
                        ft_freqgrandaverage([],SEFITC_cmb{:,1}));

%%
% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals = [stat.posclusters(:).prob];
% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep = 0.01;		% timestep between time windows for each subplot (in seconds)
sampling_rate = 1000;	% Data has a temporal resolution of 300 Hz
sample_count = length(stat.time);
					% number of temporal samples in the statistics object
j = [0:timestep:0.1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
for k = 1:10
     subplot(3,2,k);
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
     cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel reaches this significance, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).
   
   % Next, check which channels are significant over the
   % entire time interval of interest.
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     neg_int = all(neg(:, m(k):m(k+1)), 2);

     cfg.highlight = 'on';
   % Get the index of each significant channel
     cfg.highlightchannel = find(pos_int | neg_int);
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     cfg.layout              = 'neuromag306cmb.lay';
     ft_topoplotTFR(cfg, raweffectFICvsFC);   
end

%%
hFig = figure(1)
set(hFig, 'Position', [30 -100 1800 1200]);

%%
cfg = [];
cfg.frequency = [70 100];
cfg.avgoverfreq = 'yes';
freq_plot = ft_selectdata(cfg,stat_log)

cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.channel             = {'all','-MEG1822+1823','-MEG0732+0733'};
% cfg.xlim = [0.02 0.2];
cfg.maskparameter       = 'mask';
cfg.shading             = 'interp';
cfg.maskstyle           = 'box';
% cfg.frequency   = [40 90];
cfg.layout              = 'neuromag306cmb.lay';
% ft_multiplotER(cfg, freq_plot);

freq_plot.posclusters  = stat_log.posclusters;

cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
% cfg.toilim = [0.02 0.2];
cfg.zlim   = [-14 14];
cfg.layout              = 'neuromag306cmb.lay';
ft_clusterplot(cfg, freq_plot);

%%
loop=1;
for subj = [1 3:10]
    for cond = 1:6
   
        tlckBsl_10{loop,cond} = tlckBsl{subj, cond};
        
        
    end
    loop = loop+1;
end

figure,anova_rm(updrs)
