clear all;   clc;    close all;
% 
%%
% cd('J:\MEG_Research\SEF\SEFGammaControl\freq')
% load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20.mat')
load('J:\MEG_Research\SEF\SEFepresults\control\maxPosN20Control.mat')

load('J:\MEG_Research\SEF\neighbours.mat')
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 6, []);
file_sub = file_sub';
[row col] = size(file_sub);

cfg_toi               = [];
cfg_toi.latency       = [0.02 0.08];


cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 1:row
    disp(' ')
    disp(' ')
    disp(['######   ', num2str(subj)])
    disp(' ')
    disp(' ')
    for cond = 1:col
        data = load(char(file_sub(subj, cond)));
    %    
        cfg = [];
        cfg.baseline                = [-0.15 -0.1];
        cfg.baselinetype            = 'relchange';
        cfg.parameter               = 'powspctrm';
        gammaFreq_bsl               = ft_freqbaseline(cfg,data.gammaFreq);
    
        
        cfg                         = [];
        gammaCmb                    = ft_combineplanar(cfg, gammaFreq_bsl);
    
         
        cfg             = [];
        cfg.jackknife   = 'yes';
        gammaDesc       = ft_freqdescriptives(cfg, gammaCmb);
        
     %
  
        topoGamma                 = ft_selectdata(cfg_toi, gammaCmb);
        maxGamma(subj,cond,:)     = max(topoGamma.powspctrm,[],3);

        maxGammaStruct{subj,cond} = ...
                       ft_selectdata(cfg_sing, gammaCmb);

        maxGammaStruct{subj,cond}.powspctrm = max(topoGamma.powspctrm,[],3);
        maxGammaStruct{subj,cond}.name= char(file_sub(subj, cond));
                   
%         gammaValCh(subj,cond,:) = max(squeeze(data.gammaFreq_bsl.powspctrm(:,:,toi)),[],2);
        
        allGamma{subj, cond} = gammaCmb;
        tok = tokenize(char(file_sub(subj, cond)),'-');
        allGamma{subj,cond}  = tok{1};
        
        


    end
end

blankForPlot = maxGammaStruct{1,1};
blankForPlot.powspctrm = zeros(204,1);
%

% chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
%     'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633',...
%     'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623', 'MEG1812+1813',...
%     'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

% chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
%     'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG1612+1613',...
%     'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823'};
% chan_CMC_paper = {'MEG0222+0223', 'MEG0232+0233', 'MEG0412+0413', ...
%     'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
%     'MEG0712+0713', 'MEG0742+0743', 'MEG1622+1623', 'MEG1632+1633', ...
%     'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

% for loop = 1:length(chan_CMC_paper) 
% %     chanPos(loop) =  find(strcmp(blankForPlot.label,chan_CMC_paper(loop)));
%     chanPos(loop) ={neighbours(1,maxPosN20).neighblabel{:}, neighbours(1,60).label};
% end

%
gammaCmb = maxGammaStruct{1,1};
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
% cfg.highlightsize    = 8;
% cfg.highlightfontsize    = 8;
% cfg.highlightsymbol  = 'o';

for subj = 1:row
    for cond = 1:col
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(gammaCmb.label, chan_sel);
        [maxVal maxPoschk] = max(maxGamma(subj, cond,chanPos), [], 3);
%         [maxVal maxPoschk] = max(maxGammaStruct{subj,cond}.powspctrm(chanPos));
%         [maxVal] = mean(maxGamma(subj, cond,chanPos), 3);
        maxPos(subj,cond) = chanPos(maxPoschk);
%         cfg.highlightchannel = {maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}};
%         cfg.highlight        =  'numbers';
%         cfg.channel = {'all', '-MEG1612+1613', '-MEG2342+2343'};
%         cfg.zlim                = [0 1.5];
        hFig = figure(subj);
        set(hFig, 'Position', [10 80 1824 968]);
        subplot(2,3,cond),ft_topoplotER(cfg, maxGammaStruct{subj,cond})
        title({maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}})
    end
%    suptitle(maxGammaStruct{subj,1}.name(1:3))
end
 %
for subj = 1:12
    for cond = 1:6
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(maxGammaStruct{subj,cond}.label, chan_sel);
%        maxGammaChExtract(subj, cond) =  max(maxGammaStruct{subj,cond}.powspctrm(chanPos))
%        maxGammaChExtract(subj, cond) = max(maxGamma(subj,cond,maxPosN20(subj,cond)),[],3);
        maxGammaChExtract(subj, cond) = mean(maxGamma(subj,cond,chanPos),3);
    end
end
gammaPDmean_12 = maxGammaChExtract;
figure,boxplot(gammaPDmean_12([1 3:5 7:10 12],:))

%%
contS = [1:8 9 10]
subjP = [1  3:5 7:10 12]
[p, table] = anova_rm({gammaPDmean_12(subjP,:), gammaControlmean(contS,:)});

%% [p, table, stats] = anova_rm(gammaPD([1 3:10],:));
figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.05,'ctype','hsd');
% title('log')

%% Plot gamma
cfg                     = [];
cfg.parameter           = 'powspctrm';
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
cfg.xlim = [0.02 0.08]
cfg.graphcolor = 'rbgkcm';
cfg.zlim=[0.6 1.5];

% for subj = [1 2 3:row]
% subj = 2
% cond = 1
    hFig = figure(1);
    set(hFig, 'Position', [10 80 1824 968]);
%         
for cond = 1:col
   subplot(2,3,cond),

        ft_topoplotER(cfg, GA_gamma_Control{1,cond})
        title(['ControlTopo--', num2str(cond)])
%     end
%         suptitle(['Control--', num2str(subj)])

     set(gcf, 'Color', 'w'); % white bckgr
        export_fig( gcf, ...      % figure handle
            ['ControlTopoGA'],... % name of output file without extension
            '-painters', ...      % renderer
            '-png', ...           % file format
            '-r250' );             % resolution in dpi
end

%%

plot(mean(ControlSEF_ITC_summary_1([1:7],:)), '-ob'), hold on,
plot(mean(SEF_ITC_summary_2), '-or'), hold on,



%%

for  i = 1:6
   
    [h p ci stats] = ttest(gammaPD(:,i), ControlGamma(:,i));
    ci
end

%%

figure(20),
subplot 221, qqplot(1./(reshape(PDgammaFromTFR,[],1))),title('log')
subplot 222, qqplot((reshape(PDgammaFromTFR,[],1))),title('raw')
subplot 223, histfit(1./(reshape(PDgammaFromTFR,[],1)), 20,'normal'),title('log')
subplot 224, histfit(reshape(PDgammaFromTFR,[],1) ,20, 'normal'),title('raw')
suptitle('PD')

figure(2),
subplot 221, qqplot(1./(reshape(ControlGammaTFR,[],1))),title('log')
subplot 222, qqplot((reshape(ControlGammaTFR,[],1))),title('raw')
subplot 223, histfit(1./(reshape(ControlGammaTFR,[],1)), 10,'normal'),title('log')
subplot 224, histfit(reshape(ControlGammaTFR,[],1) ,10, 'normal'),title('raw')
suptitle('Control')

%%
figure
plot(SEF_ITC_summary_2(subj,:)', '-o')
axis([0.5 6.5 0 0.6])

figure, plot(gammaPD(subj, [1,6])','o')
lsline
axis([0.5 2.5 0.4 1.4])
title('DBS ON vs MED ON')
%%
figure, 
plot(updrs([3 5 7:10],[1 2 3 7])', '-o')
axis([0.5 4.5 0 45])

%%
figure,
scatter([gammaPD(subj,1)-gammaPD(subj,6)],[SEF_ITC_summary_2(subj,1)-SEF_ITC_summary_2(subj,6)] )
lsline

%%
hFig = figure(1);
set(hFig, 'Position', [10 80 1200 800]);
boxplot(maxGammaPD,'Labels',{'DBS ON','DBS OFF(0 min)','DBS OFF(60 min)',...
                     'DBS OFF(90min)', 'DBS OFF(120min)', 'MED ON'})
title('Early gamma augmentation - PD')
xlabel('Timeline'), ylabel('Relative increase in Gamma power a.u.')
axis([0.5 6.5 0.4 1.5])

hFig = figure(2),
set(hFig, 'Position', [10 80 1200 800]);
boxplot(maxGammaControl,'Labels',{'#1','#2','#3','#4','#5','#6',})
title('Early gamma augmentation - Control')
xlabel('Timeline'), ylabel('Relative increase in Gamma power a.u.')
axis([0.5 6.5 0.4 1.5])
%%

distributionPlot(maxGammaChExtract)
distributionPlot(maxGammaChExtract,'colormap',copper,'showMM',6,'histOpt',2); 

%%
plot(ControlGamma(:,:)','o',...
                'LineWidth',2,...
                'MarkerSize',10)
axis([0.5 6.5 0 2])
set(gca,'XTickLabel',{'DBS ON', 'DBS OFF (0)','DBS OFF (60)', 'DBS OFF(120)', 'MED ON' })
set(gca,'XTick',[1:5]) 
lsline

%%
for subj = 1:9
    [R,p] = corrcoef(1:2,maxGammaChExtract([subj],[1 2]),'rows','pairwise');
    % R2 = R(1,2).^2
    pval(subj) = str2num(sprintf('%1.3f',p(1,2)));
    Rval(subj) = R(1,2);
    
end
Rval
pval
%% save figures

% for subj=1:9
%     figure(subj)
 set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...      % figure handle
    ['clustT-cond2vs4'],... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', ...           % file format
    '-r250' );             % resolution in dpi
% end

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
[h p ci stats]= ttest(gamma_baseline.powspctrm,gamma_activation.powspctrm,0.01/(4*1e5));
% [h p ci stats]= ttest(log(gammaFreq_bsl.powspctrm),0.01);
% imagesc(gammaFreq_bsl.time(26:40),1:102,squeeze(p(:,:,:,26:40)))
imagesc(squeeze(stats.tstat)>8), colorbar

%%
tval = stats.tstat.*(stats.tstat(:,:,:,:)>36);
gammaFreq_bsl.tval = tval;

%%
cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relative';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);
%%
cfg                     = [];
% cfg.parameter           = 'tval';
cfg.channel             = {'all', '-MEG0632+0633'}
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.0 0.18];
cfg.zlim                = [0 1.5];
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
figure,ft_multiplotER(cfg, gammaCmb)
% title('ER')

%%

cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [-0.02 0.16];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
%               figure(subj),
ft_multiplotER(cfg, allGamma{1,:})

%%
cfg           = [];
cfg.channel   = {'MEG*2','MEG*3'};
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

subj = [1:10];
baseline = repmat(updrs(subj,1),1,2);
updrsNorm = 100*(updrs(subj,[1 2])-baseline)./baseline;
updrs_res = reshape(updrsNorm(:,[2]),1,[]);

% updrsM = median(updrsNorm,2);

baseline = repmat(gammaPDmean(subj,1),1,2);
gammaNorm = 100*(gammaPD(subj,[1 2])-baseline)./baseline;
gamma_res = reshape(gammaNorm(:,[2])', 1,[]);

% baseline = repmat(currents(subj,1),1,3);
% currNorm = (currents(subj,[1 2 7]));%-baseline)./baseline;
% curr_res = reshape(currNorm', 1,[]);

% baseline = repmat(max77115LatExtract(subj,1),1,3);
% max77115N = (max77115LatExtract(subj,[1 2 6])-baseline);%./baseline;
% max77115N_res = reshape(max3040LatExtract(subj, [1 2 6]), 1,[]);
disp('#######')
[R,p] = corr(updrs_res', gamma_res','type', 'pearson', 'tail', 'both',...
    'rows', 'pairwise')
% R2 = R(1,2).^2
% pvalue = str2num(sprintf('%1.3f',p))

figure,scatter(updrs_res, gamma_res, ...
    'jitter','on', 'jitterAmount',0.2,  'LineWidth',2, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 .1 .63]),
lsline;
legend('Data',['Rho: ', num2str(R), char(10),' p=', num2str(p)])
title('Gamma power vs UPDRS III score - DBS ON, DBS OFF(0), MED ON')
xlabel('Change(%) UPDRS III'), ylabel('Change (%) in Gamma power') 
% 
% figure(2), scatter(curr_res, gamma_res), lsline
% legend('Data',['Rho: ', num2str(corr2(curr_res, gamma_res))])
% title('Gamma power vs Sitmulation current')
% xlabel('Current mA'), ylabel('Gamma power - abs normalized') 

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

fig = figure(1);
time_vector = [1 2 4 5 6 8];

set(fig, 'Position', [10 80 1024 800]);
clf(fig)
hold on;
[AX,H1,H2] = plotyy([time_vector; time_vector+0.1]' , ...
    [mean(maxITC60_100([1:10],:),1); mean(conITC60_100([1:10],:),1)]',...
    [1 2 4 8], mean(updrs([3:5 7:10],[1 2 3 7]),1));

set(H1(1), 'linestyle', '-', 'Marker', 'x', 'color', 'r', 'linewidth', 2);
set(H1(2), 'linestyle', '-', 'Marker', 'o', 'color', 'b', 'linewidth', 2);
set(H2, 'linestyle', '--','Marker', 'v',  'color', 'k', 'linewidth', 2);

errorbar(time_vector, mean(PDITC_meanNeigh([1:10],:),1),  sem(PDITC_meanNeigh([1:10],:)), 'r.');
errorbar(time_vector+0.1, ...
    mean(conITC_meanNeigh([1:10],:),1),  sem(conITC_meanNeigh([1:10],:)), 'b.');

set(fig, 'CurrentAxes', AX(2));
hold on
errorbar([1 2 4 8], mean(updrs([3:5 7:10],[1 2 3 7]),1), ...
    sem(updrs([3:5 7:10],[1 2 3 7])), '.k');

set(AX(1),'YLim',[0.4 1.5]);
set(AX(1),'YTick',[0.4:0.1:1.1], 'fontsize', 14, 'fontname', 'Georgia');
set(AX(2),'YLim',[0 200]);
set(AX(2),'YTick',[0:5:40],'fontsize', 14, 'fontname', 'Georgia');
set(AX(2), 'YColor', [0 0 0]);

set(gca,'XLim',[0.5 8.5]);
set(gca,'XTick',[1 2 4 5 6 8]);
set(gca,'YLim',[0.9 2]);
set(gca,'YTick',[0.9:0.25:1.95]);

set(AX(2),'XLim',[0.5 8.5]);
set(AX(2),'XTick',[1:1:8]);

xTick = 1:8;
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, '',{'DBS OFF','(60min)'},...
    {'DBS OFF','(90min)'},{'DBS OFF','(120min)'}, '','Med ON'};
for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.15*(yTick(end)-yTick(1)),xTickLabel{k},...
        'HorizontalAlignment','center', 'fontsize', 14, 'fontname', 'Georgia');
end


% 
% set(AX(1),'XTickLabel',...
%     {{'DBS'; 'ON'}, 'DBS OFF (0min)', '','DBS OFF (60min)', 'DBS OFF (0min)',...
%     'DBS OFF (120min)', '','Med ON'},'fontsize', 14, 'fontname', 'Georgia') ;
set(AX(1),'XTickLabel',{''}) ;
set(AX(2),'XTickLabel',{''}) ;

hxlabel = xlabel('Conditions');
set(hxlabel,'Fontsize',20);
set(hxlabel,'Fontangle','italic');
set(hxlabel,'Fontname','Georgia');

set(fig, 'CurrentAxes', AX(1));
hylabel = ylabel('Mean of Gamma 20-80 ms');
set(hylabel,'Fontsize',20);
set(hylabel,'Fontangle','italic');
set(hylabel,'Fontname','Georgia');

set(fig, 'CurrentAxes', AX(2));
hylabel = ylabel('UPDRS III');
set(hylabel,'Fontsize',20);
set(hylabel,'color', 'k');
set(hylabel,'Fontangle','italic');
set(hylabel,'Fontname','Georgia');

htitle = title('Early somatosensory cortical processing');
set(htitle,'Fontsize',20);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Georgia');

hlegend = legend([H1(1) H1(2) H2], {'PD','Control','UPDRS III'},...
    'Interpreter','latex','FontSize',14,'FontAngle','italic',...
    'FontName','Georgia','FontWeight','bold');
% set(hlegend,'Fontsize',14);
% set(hlegend,'Fontangle','italic');
% set(hlegend,'FontWeight','bold');
% set(hlegend,'Fontname','Georgia');
legend boxoff

% htext = text(0.5,0.5,' *   p<0.05',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
htext = text(0.5,0.5,{'$\pm$'},'Interpreter','latex',...
   	'HorizontalAlignment','Center',...
   	'BackGroundColor','none','Fontsize',14,...
    'Fontangle','italic', 'Fontname','Georgia');
htext = text(0.5,0.5,'Mean',...
   	'HorizontalAlignment','Center',...
   	'BackGroundColor','none','Fontsize',14,...
    'Fontangle','italic', 'Fontname','Georgia');
htext = text(0.5,0.5,'SEM',...
   	'HorizontalAlignment','Center',...
   	'BackGroundColor','none','Fontsize',14,...
    'Fontangle','italic', 'Fontname','Georgia');
box off;

%%
figure
contS = [1:8 9 10];
subjP = [1  3:5 7:10 12];
time_vector = [1 2 4 5 6 8];

errorbar(time_vector, mean(conITC_meanNeigh(contS,:),1),  std(conITC_meanNeigh(contS,:)), '-r.'); hold on
errorbar(time_vector+0.2, mean(PDITC_meanNeigh(subjP,:),1),  std(PDITC_meanNeigh(subjP,:)), '-b.');
% set(fig, 'CurrentAxes', AX(2));
% errorbar([1 2 6 8], mean(updrs([3:5 7:10],[1 2 3 7]),1), ...
%     sem(updrs([3:5 7:10],[1 2 3 7])), '.k');



%%
[bX,b1,b2] = plotyy(time_vector+0.1, mean(ControlGamma([1:10],:),1),...
    [1 2 6 8], mean(updrs([1 3:5 7:10],[1 2 3 7]),1));

set(fig, 'CurrentAxes', AX(1));
hold on;
errorbar(time_vector, mean(gammaPD([1:10],:),1),  sem(gammaPD([1:10],:)), 'r.');

set(fig, 'CurrentAxes', AX(2));
hold on;
errorbar(time_vector+0.1, mean(gammaPD([1:10],:),1), sem(updrs([1 3:5 7:10],[1 2 3 7])), 'r.');

set(fig, 'CurrentAxes', AX(2));
xlim(H1,[0.4 1])

%%
x=[0:1:10];
[AX,H1,H2] = plotyy(x,x.^2,x,x+1)
set(AX(1),'YLim',[0 360])
set(AX(1),'YTick',[0:20:360])
set(AX(2),'YLim',[0 12])
set(AX(2),'YTick',[0:1:12])
ylabel('wind speed')
set(get(AX(2),'Ylabel'),'string','direction')

%%
figure
clf
conS = [1:8 9 10]
subjP = [1  3:10]
% time_vector = [1 2 4  8];
time_vector = [1 2 4 5 6 8];


% h = errorbar(time_vector+0.1,mean(gammaPD([1:10],:),1), sem(gammaPD([1:10],:)),...
%     '--rs', 'LineWidth',1.5,'MarkerSize',10);
% hold on
% h = errorbar(time_vector,mean(ControlGamma([1:10],:),1), sem(ControlGamma([1:10],:)),...
%     '--bs', 'LineWidth',1.5,'MarkerSize',10);

% plot([1 2 6 8],mean(updrs([1 3:5 7:10],[1 2 3 7]),1),...
%     '--ro', 'LineWidth',1.5,'MarkerSize',10)
% h = errorbar([1 2 6 8],mean(updrs([1 3:5 7:10],[1 2 3 7]),1), sem(updrs([1 3:5 7:10],[1 2 3 7])),...
%     '--ko', 'LineWidth',1);
dat_1 = conITC_MaxNeigh(conS,:);%  N20mAmpControl(conS,:)./1e-12.*10;
dat_2 = PDITC_maxNeigh(subjP,:);%N20mAmpPD(subjP,:)./1e-12.*10;
hold on
% ylim([5 40])
% plot(1:1:4,mean(updrs([1 3:5 7:10],[1 2 3 7]),1),...
%     'ro', 'LineWidth',1.5,'MarkerSize',10)
h = errorbar(time_vector+0.05,mean(dat_1,1), sem(dat_1,1),...
    '-bs', 'LineWidth',2.5,'MarkerSize',10);
hold on
h = errorbar(time_vector-0.05,mean(dat_2,1), sem(dat_2,1),...
    '-ro', 'LineWidth',2.5,'MarkerSize',10);
% 
% h = errorbar(time_vector+0.1,mean(ControlSEF_ITC_summary_1([1:9],:),1), sem(ControlSEF_ITC_summary_1([1:10],:),1),...
%     '--bs', 'LineWidth',2.5,'MarkerSize',10);
% hold on
% h = errorbar(time_vector+0.15,mean(ControlSEF_ITC_summary_2([1:10],:),1), sem(ControlSEF_ITC_summary_2([1:6],:),1),...
%     '-bs', 'LineWidth',2.5,'MarkerSize',10);

% hold on
% plot(maxGammaChExtract(subj,:)', 'rx', 'MarkerSize',10, 'LineWidth',1.5)

%% boxplot(maxGammaVal)
set(gca,'XTick',[1 2 4 5 6 8])
set(gca,'XTick',[])
set(gca,'XTickLabel',[])

set(gca,'XTickLabel',...
    {'Stim ON', 'Stim OFF (0min)', '', 'Stim OFF (60min)', 'Stim OFF (90min)'...
    'Stim OFF (120min)','Med ON'},'fontsize', 14, 'fontname', 'Georgia') ;
xTick = 1:8;
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, '',{'DBS OFF','(60min)'},...
    {'DBS OFF','(90min)'},{'DBS OFF','(120min)'}, '','Med ON'};
for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.15*(yTick(end)-yTick(1)),xTickLabel{k},...
        'HorizontalAlignment','center', 'fontsize', 24, 'fontname', 'Calibri',...
        'FontWeight','bold');
end


hxlabel = xlabel('Conditions')
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Calibri');

hylabel = ylabel('ITPC')
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Calibri');

ylim([0.15 0.55])
set(gca,'YTick',[0:0.1:1])
set(gca,'Fontsize', 24)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Calibri');

htitle = title('Inter-trial phase coherence');
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Calibri');


% axis([0.5 8.5 0.5 1]);
%     %
hlegend = legend({'Control (n=10)','PD (n=9)'  });
% hlegend = legend({'UPDRS III'},'Interpreter','latex');
set(hlegend,'Fontsize',18);
% set(hlegend,'Fontangle','italic');
set(hlegend,'Fontname','Calibri');
legend boxoff
% 
% htitle = title('Early cortical sensory processing - Gamma power');
% set(htitle,'Fontsize',20);
% set(htitle,'FontWeight','bold');
% set(htitle,'Fontname','Georgia');
% 
% hxlabel = xlabel('Conditions');
% set(hxlabel,'Fontsize',20);
% set(hxlabel,'Fontangle','italic');
% set(hxlabel,'Fontname','Georgia');

% hylabel = ylabel('Mean of ITC (ROI-PD F-clustStat)');
% hylabel = ylabel('Mean of Gamma - 20-80 ms');
% set(hylabel,'Fontsize',20);
% set(hylabel,'Fontangle','italic');
% set(hylabel,'Fontname','Georgia');

%
% groups={[1,2],[1,6], [5,6]};
% H = sigstar(groups,[0.05,0.05, 0.05],0);
% set(H,'color','r')

% htext = text(0.5,0.5,' *   p<0.05',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(0.5,0.5,{'$\pm$'},'Interpreter','latex',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(0.5,0.5,'Mean',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(0.5,0.5,'SEM',...
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
    
%%
    cfg = [];
    cfg.baseline                = [-0.09 -0.05];
    cfg.baselinetype            = 'relative';
    cfg.parameter               = 'powspctrm';
    gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);

    cfg = [];
    cfg.latency = [-0.10 -0.04];
    gamma_baseline = ft_selectdata(cfg, ft_combineplanar([],gammaFreq));
    gamma_baseline.powspctrm = log10(gamma_baseline.powspctrm);

    cfg = [];
    cfg.latency = [0.02 0.08];
    gamma_activation = ft_selectdata(cfg, ft_combineplanar([],gammaFreq));
    gamma_baseline.time = gamma_activation.time;
    gamma_activation.powspctrm = log10(gamma_activation.powspctrm);
    %%
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
    
%%
end

%%

cfg               = [];
cfg.method        = 'stats';
cfg.constantvalue = 0;

%     ntrials = size(gamma_activation.powspctrm,1);
%     design  = zeros(1,2*ntrials);
%     design(1,1:ntrials) = 1;
%     design(1,ntrials+1:2*ntrials) = 2;
%     design(2,1:ntrials) = [1:ntrials];
%     design(2,ntrials+1:2*ntrials) = [1:ntrials];
cfg.design = design(2,:);    
% cfg.latency       = [0.02 0.2];
cfg.statistic     = 'ttest_2samples_by_timepoint'; % compares the mean to zero
cfg.feedback      = 'no';
frstat            = ft_freqstatistics(cfg,TFRcmb_ev{[1 3:10],1}, TFRcmb_ev{[1 3:10],2});

imagesc(squeeze(frstat.prob)), colorbar

%%
hFig = figure(1);
set(hFig, 'Position', [10 80 1824 968]);
%         
for loop=1:6
    subplot(2,3,loop)
        cfg                     = [];
        cfg.parameter = 'powspctrm';
        cfg.colorbar = 'yes';
        cfg.layout              = 'neuromag306cmb.lay';
            cfg.graphcolor              = 'brgkycmbrgkycm';

        cfg.xlim                = [0 .1];
        cfg.shading                 = 'interp';
        cfg.maskstyle               = 'saturation';	
        ft_multiplotER(cfg, allGamma{1,:})
        title(num2str(loop))
        set(gcf, 'Color', 'w'); % white bckgr
  
end

        set(gcf, 'Color', 'w'); % white bckgr
        export_fig( gcf, ...      % figure handle
        ['updatedSEFGamma'],... % name of output file without extension
        '-painters', ...      % renderer
        '-png', ...           % file format
        '-r250' );             % resolution in dpi
        
 %%
 
figure, plot(mean(gammaPDmean(subjP,:),1))
plot(mean(gammaPDmean(contS,:),1), '-bo')
hold on,plot(mean(gammaControlmean,1), '-rs')
axis([0.5 6.5 1 1.65])
title('With 11 PD as on June 30')
legend('PD(n=9)', 'Control (n=10)')
        
        
        
        