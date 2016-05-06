clear all;   clc;    close all;

%%
% cd('J:\MEG_Research\SEF\SEFGamma3tap')
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 6, []);
file_sub = file_sub';
[row col] = size(file_sub);

%

cfg_toi               = [];
cfg_toi.latency       = [0.02 0.08];


cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 11
    disp(['######   ', num2str(subj)])

    for cond = 1:col
        data = load(char(file_sub(subj, cond)));
        
        cfg = [];
        cfg.baseline                = [-0.10 -0.06];
        cfg.baselinetype            = 'relchange';
        cfg.parameter               = 'powspctrm';
        gammaFreq_bsl               = ft_freqbaseline(cfg,data.gammaFreq);
     
        cfg             = [];
        cfg.jackknife   = 'yes';
        gammaDesc       = ft_freqdescriptives(cfg, gammaFreq_bsl);
        
        cfg                         = [];
%         cfg.method = 'svd' ;
        gammaCmb               = ft_combineplanar(cfg, gammaDesc);
        
     
  
        topoGamma         = ft_selectdata(cfg_toi, gammaCmb);
        maxGamma(subj,cond,:)     = max(topoGamma.powspctrm,[],3);

        maxGammaStruct{subj,cond} = ...
                       ft_selectdata(cfg_sing, gammaCmb);

        maxGammaStruct{subj,cond}.powspctrm = max(topoGamma.powspctrm,[],3);
        maxGammaStruct{subj,cond}.name= char(file_sub(subj, cond));
                   
%         gammaValCh(subj,cond,:) = max(squeeze(data.gammaFreq_bsl.powspctrm(:,:,toi)),[],2);
        
        allGamma{subj, cond} = gammaCmb;
        
        
%         cfg                     = [];
%         cfg.parameter = 'powspctrm';
%         cfg.colorbar = 'yes';
%         cfg.layout              = 'neuromag306cmb.lay';
%         cfg.xlim                = [-0.02 0.16];
%         cfg.shading                 = 'interp';
%         cfg.maskstyle               = 'saturation';	
% %               figure(subj),
%         ft_multiplotER(cfg, allGamma{7,3})
        title(num2str(cond));

    end
end

blankForPlot = maxGammaStruct{1,1};
blankForPlot.powspctrm = zeros(204,1);


% chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
%     'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633',...
%     'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623', 'MEG1812+1813',...
%     'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
    'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG1612+1613',...
    'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823'};
chan_CMC_paper = {'MEG0222+0223', 'MEG0232+0233', 'MEG0412+0413', ...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
    'MEG0712+0713', 'MEG0742+0743', 'MEG1622+1623', 'MEG1632+1633', ...
    'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

for loop = 1:length(chan_CMC_paper) 
%     chanPos(loop) =  find(strcmp(blankForPlot.label,chan_CMC_paper(loop)));
    chanPos(loop) ={neighbours(1,60).neighblabel{:}, neighbours(1,60).label};
end

cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';

for subj = 11
    [maxVal maxPoschk] = max(maxGamma(subj, :,chanPos), [], 3);
    maxPos(subj,:) = chanPos(maxPoschk);
%     maxVal = max(maxGamma(subj, :,chanPos), [], 3);
    for cond = 1:col
        cfg.highlightchannel =  {maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}};
        cfg.highlight        =  'numbers';
        hFig = figure(subj);
        set(hFig, 'Position', [10 80 1824 968]);
        subplot(2,3,cond),ft_topoplotER(cfg, maxGammaStruct{subj,cond})
        title({maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}})
    end
   suptitle(num2str(subj))
end
 %%
for subj = 1:row
    for cond = 1:col
        maxGammaChExtract(subj, cond) = max(maxGamma(subj,cond,maxPos(subj,cond)));
    end
end
figure,boxplot(maxGammaChExtract([1 3:9],:))

%%

[p, table, stats] = anova_rm(maxGammaChExtract([1 2:9],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.02,'ctype','dunn-sidak');



%% 
cfg = [];
cfg.baseline                = [-0.09 -0.05];
cfg.baselinetype            = 'relative';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);


% chan_left = {'MEG0112+0113', 'MEG0122+0123', 'MEG0132+0133',...
%     'MEG0142+0143', 'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233',...
%     'MEG0242+0243', 'MEG0312+0313', 'MEG0322+0323', 'MEG0332+0333', ...
%     'MEG0342+0343', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433',...
%     'MEG0442+0443', 'MEG0512+0513', 'MEG0522+0523', 'MEG0532+0533', ...
%     'MEG0542+0543', 'MEG0612+0613', 'MEG0622+0623', 'MEG0632+0633', ...
%     'MEG0642+0643', 'MEG0712+0713', 'MEG0742+0743', 'MEG0812+0813', ...
%     'MEG0822+0823', 'MEG1012+1013', 'MEG1512+1513', 'MEG1522+1523',...
%     'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623',...
%     'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723',...
%     'MEG1732+1733', 'MEG1742+1743', 'MEG1812+1813', 'MEG1822+1823', ...
%     'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923',...
%     'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2042+2043',...
%     'MEG2112+2113', 'MEG2122+2123', 'MEG2142+2143'};
chan_left = {'MEG0112+0113', 'MEG0122+0123', 'MEG0132+0133', 'MEG0142+0143'};%,...
%     'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', 'MEG0312+0313',...
%     'MEG0322+0323', 'MEG0332+0333', 'MEG0342+0343', 'MEG0412+0413', 'MEG0422+0423',...
%     'MEG0432+0433', 'MEG0442+0443', 'MEG0512+0513', 'MEG0522+0523', 'MEG0532+0533',...
%     'MEG0542+0543', 'MEG0612+0613', 'MEG0622+0623', 'MEG0632+0633', 'MEG0642+0643',...
%     'MEG0712+0713', 'MEG0742+0743', 'MEG0812+0813', 'MEG0822+0823', 'MEG1012+1013',...
%     'MEG1512+1513', 'MEG1522+1523', 'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613',...
%     'MEG1622+1623', 'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723',...
%     'MEG1732+1733', 'MEG1742+1743', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833',...
%     'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923', 'MEG1932+1933', 'MEG1942+1943',...
%     'MEG2012+2013', 'MEG2042+2043', 'MEG2112+2113', 'MEG2122+2123', 'MEG2142+2143'};

cfg             = [];
cfg.jackknife   = 'yes';
%  cfg.variance      = 'yes';
gammaDesc       = ft_freqdescriptives(cfg, gammaFreq_bsl);

cfg                         = [];
gammaFreq_bsl               = ft_combineplanar(cfg, gammaFreq_bsl);

%%
[h p ci stats]= ttest(log(gammaFreq_bsl.powspctrm),0.01);
% imagesc(gammaFreq_bsl.time(26:40),1:102,squeeze(p(:,:,:,26:40)))
imagesc(gammaFreq_bsl.time(51:61),1:102,squeeze(stats.tstat(:,:,:,51:61))>38), colorbar

%%
tval = stats.tstat.*(stats.tstat(:,:,:,:)>36);
gammaFreq_bsl.tval = tval;

%%
% cfg = [];
% cfg.baseline                = [-0.09 -0.05];
% cfg.baselinetype            = 'relchange';
% cfg.parameter               = 'powspctrm';
% gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);

cfg                     = [];
% cfg.parameter           = 'tval';
% cfg.channel             = {'all', '-MEG0222+0223'}
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
% cfg.zlim                = [1 3.2];
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
% cfg.interplimits = 'electrodes';
figure,ft_topoplotER(cfg, gammaDesc)

%%
cfg           = [];
% cfg.channel   = {'MEG*2','MEG*3'};
cfg.method    = 'distance';
% cfg.template = 'neuromag306planar_neighb.mat';
cfg.neighbourdist = 6;
cfg.grad      = gammaCmb.grad;
cfg.feedback  = 'yes';
neighbours    = ft_prepare_neighbours(cfg);

cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.channel = {'all', '-MEG1512+1513'};
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';
cfg.highlightchannel =  {neighbours(1,60).neighblabel{:}, neighbours(1,60).label};
cfg.highlight        =  'numbers';
figure,ft_topoplotER(cfg, gammaFreq_bsl)

%%
cfg=[];
cfg.latency=[0 0.15];
cfg.method='stats';
% cfg.method    = 'analytic';
cfg.correctm  = 'bonferoni';
cfg.statistic='ttest';
cfg.alpha            = 0.01./(102*101);
cfg.design = [ones(1,size(gammaFreq_bsl.powspctrm,1))];
[stat] = ft_freqstatistics(cfg, ft_combineplanar([],gammaFreq_bsl));

imagesc(stat.time, 1:102, squeeze(stat.stat))

%%
% newGamma= [maxGammaChExtract(:,1)...
%     mean(maxGammaChExtract(:,[2,3]),2) ...
%     maxGammaChExtract(:,6)];

subj = [1 3:6 7];
baseline = repmat(updrs(subj,1),1,3);
updrsNorm = (updrs(subj,[1 2 7])-baseline);%./baseline;
updrs_res = reshape(updrsNorm', 1,[]);

baseline = repmat(maxGammaChExtract(subj,1),1,3);
gammaNorm = (maxGammaChExtract(subj,[1 3 6])-baseline);%./baseline;
gamma_res = reshape(gammaNorm', 1,[]);

baseline = repmat(currents(subj,1),1,3);
currNorm = (currents(subj,[1 2 7]));%-baseline)./baseline;
curr_res = reshape(currNorm', 1,[]);

% baseline = repmat(max77115LatExtract(subj,1),1,3);
% max77115N = (max77115LatExtract(subj,[1 2 6])-baseline);%./baseline;
% max77115N_res = reshape(max3040LatExtract(subj, [1 2 6]), 1,[]);
[R,p] = corrcoef(updrs_res,gamma_res,'rows','pairwise')
% R2 = R(1,2).^2
pvalue = str2num(sprintf('%1.3f',p(1,2)))

figure(1),scatter(updrs_res, gamma_res, 'jitter','on', 'jitterAmount',0.5), lsline;
legend('Data',['Rho: ', num2str(R(1,2)), char(10),' p=', num2str(pvalue)])
title('Gamma power vs UPDRS III score')
xlabel('UPDRS III - abs normalized'), ylabel('Gamma power - abs normalized') 
% 
figure(2), scatter(curr_res, gamma_res), lsline
legend('Data',['Rho: ', num2str(corr2(curr_res, gamma_res))])
title('Gamma power vs Sitmulation current')
xlabel('Current mA'), ylabel('Gamma power - abs normalized') 

% figure(3), scatter(updrs_res, max77115N_res), lsline
% legend('Data',['Rho: ', num2str(corr2(updrs_res, max77115N_res))])
% title('77-115 ms latency vs Sitmulation current')
% xlabel('UPDRS'), ylabel('77-115 Peak Latency ms') 

%%
figure,plot(maxGammaChExtract', 'o')
% for loop=  1:7
%     [R,p] = corrcoef([1:6],maxGammaChExtract(loop,:),'rows','pairwise');
% % R2 = R(1,2).^2
%     pvalue(loop) = p(1,2);
%     Rval(loop) = R(1,2);
% end
lsline
axis([0.5 6.5 0 3])

%%
% dataGamma(6,6).pow = [];

for subjLoop = 1:row
    disp(['######   ', num2str(subjLoop)])

    for condLoop = 1:col
        data_gamma_trials = [];
        disp(['******** ', char(file_sub(subjLoop, condLoop))])
        load(char(file_sub(subjLoop, condLoop)))
        trialsGamma     = [];
        trialsGamma     =    ...
            max(max(max(data_gamma_trials.powspctrm, [], 4), [], 3), [], 2);
        dataGamma(subjLoop, condLoop).pow = trialsGamma;
        maxGammaVal(subjLoop, condLoop) = ...
            max(max(squeeze(mean((data_gamma_trials.powspctrm),1))));
    end
end

%%

for subjLoop = 1:row
    for condLoop = 1:col
        
        maxGammaVal(subjLoop, condLoop) = ...
            max(dataGamma(subjLoop, condLoop).pow);
        
    end
end

%%
baseline = repmat(maxGammaChExtract(:,1),1,6);
maxGammaValNorm = (maxGammaChExtract -baseline);


%%
figure,
subj = [1 3:6];
errorbar(1:6,mean(maxGammaChExtract(subj,:),1), std(maxGammaChExtract(subj,:),0,1),...
    '-ks', 'LineWidth',1.5,'MarkerSize',10);
hold on
plot(maxGammaChExtract(subj,:)', 'bx', 'MarkerSize',6, 'LineWidth',1.5)

% boxplot(maxGammaVal)
set(gca,'XTickLabel',...
    {'Stim ON', 'Stim OFF (0min)', 'Stim OFF (60min)', 'Stim OFF (90min)',...
    'Stim OFF (120min)','Med ON'},'fontsize', 9, 'fontname', 'Georgia') ;

axis([0.5 6.5 0.5 3.5]);
    
hlegend = legend({'Mean $\pm$ SD'},'Interpreter','latex');
hlegend = legend('Mean');
set(hlegend,'Fontsize',14);
set(hlegend,'Fontangle','italic');
set(hlegend,'Fontname','Georgia');
legend boxoff

htitle = title('Longitudinal evolution of Gamma power');
set(htitle,'Fontsize',18);
% set(htitle,'Fontangle','italic')
set(htitle,'Fontname','Georgia');

hxlabel = xlabel('Timeline');
set(hxlabel,'Fontsize',16);
set(hxlabel,'Fontangle','italic');
set(hxlabel,'Fontname','Georgia');

hylabel = ylabel('Baseline-corrected Gamma power');
set(hylabel,'Fontsize',16);
set(hylabel,'Fontangle','italic');
set(hylabel,'Fontname','Georgia');

% groups={[1,2],[1,6]};
% H = sigstar(groups,[0.025,0.025],0);

% htext = text(3.5,3.5,' *   p<0.025',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(3.5,3.5,{'$\pm$'},'Interpreter','latex',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(3.5,3.5,'SD',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% box off
% 


%%
cfg           = [];
cfg.channel   = {'MEG*2','MEG*3'};
cfg.method    = 'triangulation';
cfg.grad      = gammaFreq.grad;
cfg.feedback  = 'yes';
neighbours    = ft_prepare_neighbours(cfg);

%%
filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
%     cfg = [];
%     cfg.baseline                = [-0.15 -0.09];
%     cfg.baselinetype            = 'relchange';
%     cfg.parameter               = 'powspctrm';
%     gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);
% 
%     cfg                         = [];
%     gammaFreq_bsl               = ft_combineplanar(cfg, gammaFreq);
    
%
    cfg = [];
cfg.baseline                = [-0.09 -0.05];
cfg.baselinetype            = 'relative';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);

    cfg = [];
    cfg.latency = [-0.10 -0.04];
    gamma_baseline = ft_selectdata(cfg, ft_combineplanar([],gammaFreq));
%     gamma_baseline.powspctrm = log10(gamma_baseline.powspctrm);

    cfg = [];
    cfg.latency = [0.02 0.08];
    gamma_activation = ft_selectdata(cfg, ft_combineplanar([],gammaFreq));
    gamma_baseline.time = gamma_activation.time;
%     gamma_activation.powspctrm = log10(gamma_activation.powspctrm);
    %
    cfg = [];
    cfg.latency          = [0.02 0.08];
    % cfg.parameter = 'trial';
    cfg.method           = 'montecarlo';
%     cfg.frequency        = [65];
    cfg.statistic        = 'ft_statfun_actvsblT';
%     cfg.avgovertime = 'yes';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 1;
    cfg.tail             = 1;
    cfg.clustertail      = 1;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 500;
    % prepare_neighbours determines what sensors may form clusters
    cfg.neighbours       = neighbours;

    ntrials = size(gamma_activation.powspctrm,1);
    design  = zeros(2,2*ntrials);
    design(1,1:ntrials) = 1;
    design(1,ntrials+1:2*ntrials) = 2;
    design(2,1:ntrials) = [1:ntrials];
    design(2,ntrials+1:2*ntrials) = [1:ntrials];

    cfg.design   = design;
    cfg.ivar     = 1;
    cfg.uvar     = 2;
    [stat] = ft_freqstatistics(cfg, gamma_activation, gamma_baseline);
    
    figure,imagesc(squeeze(stat.prob)), colorbar

    cfg = [];
    cfg.alpha  = 0.05;
    cfg.parameter = 'stat';
    cfg.zlim   = [-20 20];
    cfg.layout              = 'neuromag306cmb.lay';
    ft_clusterplot(cfg, stat);
%      suptitle(b)
    
%
end

%%

cfg                 =   [];
cfg.method          =   'stats';
cfg.tail            =    0;
cfg.alpha           =   0.05;
cfg.parameter       =   'powspctrm';
cfg.feedback        =   'no'; 
nsubj               =    size(gammaFreq_bsl.powspctrm,1);
cfg.design(1,:)     =   [ones(1,nsubj)];
cfg.constantvalue   =   0;
cfg.statistic       =   'ttest_samples_vs_const'; % compares the mean to zero
frstat              =   ft_freqstatistics(cfg,gammaFreq_bsl);

imagesc(squeeze(stats.tstat)>45)

%%

h= []; p =[]; ci =[]; stats = [];
for chan = 1:204,
    [h(:,:,:,:), p, ci, stats] = ...
        ttest(gammaFreq_bsl.powspctrm(:,chan, :,:), 0, 0.05, 0);
    mask(chan,:,:,:)   = h;
    tval(chan,:,:,:)   = stats.tstat;
end

%%
% now plot 1-probability (1 = sig, less than 0.95 not sig)
cfg=[];
cfg.layout              = 'neuromag306cmb.lay';
frstat.powspctrm=1-frstat.prob;
cfg.zlim=[0.999 1]
cfg.interactive='yes';
% fig3=figure;
% set(fig3,'Position',[0,0,800,800]);
ft_multiplotTFR(cfg, frstat);
