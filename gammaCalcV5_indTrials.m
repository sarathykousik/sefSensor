clear all;   clc;    close all;
 
%
filenames = dir('*.mat');
% load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20PD.mat')
load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20PD-before newtSSS.mat')

maxPosN20 = maxPosN20PD;
save_folder = 'J:\MEG_Research\SEF\SEFGammaStatistics\';

% load('J:\MEG_Research\SEF\SEFepresults\control\maxPosN20Control.mat')
% maxPosN20 = maxPosN20Control;

load('J:\MEG_Research\SEF\neighbours.mat')

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
    %    
        cfg = [];
        cfg.baseline                = [-0.15 -0.1];
        cfg.baselinetype            = 'relchange';
        cfg.parameter               = 'powspctrm';
        gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);
       
        cfg                         = [];
        gammaCmb                    = ft_combineplanar(cfg, gammaFreq_bsl);
        
        cfg             = [];
        cfg.jackknife   = 'yes';
        gammaDesc       = ft_freqdescriptives(cfg, gammaCmb);
        
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(gammaCmb.label, chan_sel);
        gammaC_trials{subj, cond}.powspctrm = ...
            mean(...
            max(gammaCmb.powspctrm(:,chanPos,:,find(gammaCmb.time>=0.02 & gammaCmb.time<=0.09)),[],4),2);
        gammaC_trials{subj,cond}.name  = tok{1};
        

        gammaC(subj, cond) = ...
            mean(...
            max(gammaDesc.powspctrm(chanPos,:,find(gammaCmb.time>=0.02 & gammaCmb.time<=0.09)),[],3),1);
        

        % Plotting
        
        cfg_toi               = [];
        cfg_toi.latency       = [0.02];
        topoGammaControl_a{subj,cond}  = ft_selectdata(cfg_toi, gammaDesc);
%         topoGammaPD{subj,cond}.powspctrm   = max(gammaDesc.powspctrm(:,:,find(gammaCmb.time>=0.02 & gammaCmb.time<=0.09)),[],3);
        topoGammaControl_a{subj,cond}.powspctrm = mean(max(gammaCmb.powspctrm(:,:,:,find(gammaCmb.time>=0.02 & gammaCmb.time<=0.09)),[],4));       
%         cfg                     = [];
%         cfg.parameter = 'powspctrm';
%         cfg.colorbar = 'yes';
%         cfg.layout              = 'neuromag306cmb.lay';
% %         cfg.xlim                = [0.02 0.08];
%         cfg.shading                 = 'interp';
%         cfg.maskstyle               = 'saturation';	
%         
%         hFig = figure(subj);
%         set(hFig, 'Position', [10 80 1824 968]);
%         subplot(3,3,cond),ft_topoplotER(cfg, topoGammaPD{subj,cond})
%         title(num2str(cond))
%                 
    end
%     suptitle(tok{1}(1:3))
%     export_fig( gcf,[save_folder,'PDTopo', num2str(tok{1}(1:3))],'-painters','-png','-r250' );  % resolution in dpi
end

%%
load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20PD-before newtSSS.mat')
% load('J:\MEG_Research\SEF\SEFepresults\control\maxPosN20Control.mat')
maxPosN20 = maxPosN20PD;
load('J:\MEG_Research\SEF\neighbours.mat')

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
        
        topoGamma                 = ft_selectdata(cfg_toi, allGamma{subj,cond});
        maxGamma                  = max(topoGamma.powspctrm,[],3);
        
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(topoGamma.label, chan_sel);
        maxGammaChExtract(subj, cond) = mean(maxGamma(chanPos));

    end
end
gammaPD = maxGammaChExtract

%%
subj=[1:10];
diff_56 = maxGammaChExtract(subj,5)-maxGammaChExtract(subj,6);
diff_16 = maxGammaChExtract(subj,1)-maxGammaChExtract(subj,6);
diff_12 = maxGammaChExtract(subj,1)-maxGammaChExtract(subj,2);

Cdiff_56 = gammaControl_6(:,[5])-gammaControl_6(:,[6]);
Cdiff_16 = gammaControl_6(:,[1])-gammaControl_6(:,[6]);
Cdiff_12 = gammaControl_6(:,[1])-gammaControl_6(:,[2]);


figure,
meas = diff_56;
plot(maxGammaChExtract(subjP,[1 2])', '-ro'), hold on
bar(1,mean(meas), 'r'),
errorbar(1, mean(meas), sem(meas))
xlim([0.5 2.5]), ylim([-0.5 3])

% figure
meas = Cdiff_56;
plot(gammaControl_6(subjP,[1 2])', '-bx'), hold on
bar(2,mean(meas),'b'),
errorbar(2, mean(meas), sem(meas))
xlim([0.5 2.5]), ylim([-0.5 3])

%%
blankForPlot = maxGammaStruct{1,1};
blankForPlot.powspctrm = zeros(204,1);
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
for subj = 1:row
    for cond = 1:col
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(maxGammaStruct{subj,cond}.label, chan_sel);
%        maxGammaChExtract(subj, cond) =  max(maxGammaStruct{subj,cond}.powspctrm(chanPos))
%        maxGammaChExtract(subj, cond) = max(maxGamma(subj,cond,maxPosN20(subj,cond)),[],3);
        maxGammaChExtract(subj, cond) = mean(maxGamma(subj,cond,chanPos),3);
    end
end
gammaPDmean = maxGammaChExtract;
figure,boxplot(gammaPDmean)

%%
% contS = [1:10];
% subjP = [1:10];
[p, table] = anova_rm({gammaPDmean(subjP,:), gammaControlmean(contS,:)});


%% Plot gamma
cfg                     = [];
cfg.parameter           = 'powspctrm';
% cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
cfg.xlim = [0 0.2]
% cfg.graphcolor = 'rbgkcm';
% cfg.zlim=[0.9 1.5];

% for subj = [1 2 3:row]
% subj = 2
% cond = 1
%     hFig = figure(1);
% for subj = 1:10         
%     hFig = figure(2);
%     set(hFig, 'Position', [10 80 1824 968]);

%     for cond = 1:7
%         subplot(3,3,cond),

        ft_topoplotER(cfg, gammaCmb_bsl)
%         title(cond)
%     end
%     suptitle('C')
%     suptitle(['PD--', num2str(subj)])
% end

%%
     set(gcf, 'Color', 'w'); % white bckgr
        export_fig( gcf, ...      % figure handle
            ['ControlTopoGA'],... % name of output file without extension
            '-painters', ...      % renderer
            '-png', ...           % file format
            '-r250' );             % resolution in dpi



%% save figures

% for subj=1:9
%     figure(subj)
 set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...      % figure handle
    ['gamma-STD'],... % name of output file without extension
    '-painters', ...      % renderer
    '-jpg', ...           % file format
    '-r250' );             % resolution in dpi
% end


%%
cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relative';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);

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

%% multiplotER Gamma

cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [-0.02 0.16];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
%               figure(subj),
ft_multiplotER(cfg, allGamma{1,:})


%% Plot for paper - plotyy

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

%% Plot for paper - Mean=/- SEM
figure
contS = [1:10];
subjP = [1:10];
time_vector = [1 2 3 4 5 6 8];

errorbar(time_vector+0.2, mean(gammaMeanExtract(subjP,:),1),  sem(gammaMeanExtract(subjP,:)), '-r.');
hold on
errorbar(time_vector, mean(gammaCExtract(contS,:),1),  sem(gammaCExtract(contS,:)), '-b.')
ylim([2 2.8])

set(gca,'XTick',[1 2 4 5 6 8])
set(gca,'XTick',[])
set(gca,'XTickLabel',[])

% set(gca,'XTickLabel',...
%     {'Stim ON', 'StimOFF(0min)', 'StimOFF(30min)', 'StimOFF(60min)', 'StimOFF(90min)'...
%     'StimOFF(120min)', '','MedON'},'fontsize', 10, 'fontname', 'Georgia') ;
xTick = [1:6 8];
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, {'DBS OFF';'(30min)'},{'DBS OFF','(60min)'},...
    {'DBS OFF','(90min)'},{'DBS OFF','(120min)'}, 'Med ON'};
for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.15*(yTick(end)-yTick(1)-0.7),xTickLabel{k},...
        'HorizontalAlignment','center', 'Fontsize', 14, 'fontname', 'Calibri',...
        'FontWeight','bold');
end


hxlabel = xlabel('Conditions')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Calibri');

hylabel = ylabel('Gamma power % increase   (Mean \pm SEM)')
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Calibri');

set(gca,'YTick',[2:0.2:3])
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Calibri');

htitle = title('Early Somatosensory processing');
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Calibri');


% axis([0.5 8.5 0.5 1]);
hlegend = legend({'PD (n=10)','Control (n=10)'});
% hlegend = legend({'UPDRS III'},'Interpreter','latex');
set(hlegend,'Fontsize',18);
set(hlegend,'Fontangle','italic');
set(hlegend,'Fontname','Calibri');
set(hlegend,'Location','Northeast');
legend boxoff
% 
 set(gcf, 'Color', 'white'); % white bckgr
%     export_fig( gcf, ...      % figure handle
%     ['DBSOFF120vsMEDON'],... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', ...           % file format
%     '-r250' );             % resolution in dpi

% % PDF print

export_fig -painters -r600 -q101 Gamma-SEM.pdf


%% cluster stat
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


        
        
        