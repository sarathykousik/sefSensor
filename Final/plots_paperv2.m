

%% N20m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% N20m GA

cfg = [];
cfg.parameter = 'avg';

for condLoop = 1:7
    GAPD{condLoop} = ft_timelockgrandaverage(cfg, tlckCmb{:,condLoop});
    GAPD{condLoop}.name = tlckCmb{1,condLoop}.name;
end

%% Plots

% Topoplots PD
cfg                     = [];
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.03];
cfg.zlim                = [0 150]; 
cfg.shading             = 'interp';
% cfg.maskstyle           = 'saturation';	
% cfg.contournum          = 12;
% cfg.fontsize            = 20;

for loop = 1:7
    dat = GAPD{loop};
    dat.avg = GAPD{loop}.avg./1e-12.*10;
    
    hFig = figure(loop)
    set(hFig, 'Position', [10 80 968 968]);
    ft_topoplotER(cfg, dat)
    title(['PD',num2str(loop)])
    
    hcb = colorbar;
    set(hcb,'YTick',[0:50:100])
    set(hcb,'fontsize',15);
    
    set(gcf, 'Color', 'w'); % white bckgr
    export_fig( gcf, ['PDTopoGA', num2str(loop)],...
        '-painters','-png', '-r250' ); 
end

% Topoplots Control
cfg                     = [];
cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.03];
cfg.zlim                = [0 150]; 
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
% cfg.contournum          = 12;
% cfg.fontsize            = 20;

for loop = 1:7
    dat = GAControl{loop};
    dat.avg = GAControl{loop}.avg./1e-12.*10;
    
    hFig = figure(10)
    set(hFig, 'Position', [10 80 968 968]);
    ft_topoplotER(cfg, dat)
    title(['C',num2str(loop)])
    
    hcb = colorbar;
    set(hcb,'YTick',[0:50:100])
    set(hcb,'fontsize',15);
    
    set(gcf, 'Color', 'w'); % white bckgr
    export_fig( gcf, ['ControlTopoGA', num2str(loop)],...
        '-painters','-png', '-r250' ); 
end

%%  ERF single plot + topo of N20m

dat = tlckCmb{3,1};
dat.avg = tlckCmb{3,1}.avg./1e-12.*10;
figure
cfg                     = [];
cfg.channel             = {'MEG0442+0443'};
cfg.colorbar            = 'yes';
cfg.linewidth           = 4;
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0 0.1];
cfg.ylim                = [0 75];
cfg.graphcolor          = clrOp1;
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
ft_singleplotER(cfg, dat)

set(gca,'XTick',[0:0.02:0.1])
set(gca,'XTickLabel',[])
% set(gca,'Fontsize', 18, 'fontname', 'Helvetica','FontWeight','bold');


% hylabel = ylabel(['Amplitude fT-cm ^{-1}'])
% set(hylabel,'Fontsize',24);
% set(hylabel,'FontWeight','bold');
% set(hylabel,'Fontname','Helvetica');
% 
set(gca,'YTick',[0:25:100])
set(gca,'YTickLabel',[])
% set(gca,'Fontsize', 18)
% set(gca,'Fontweight','bold');
% set(gca,'Fontname','Helvetica');
% 
% hxlabel = xlabel('Time s')
% set(hxlabel,'Fontsize',24);
% set(hxlabel,'FontWeight','bold');
% set(hxlabel,'Fontname','Helvetica');
% 
% 
% htitle = title('');
% set(htitle,'Fontsize',24);`
% set(htitle,'FontWeight','bold');
% set(htitle,'Fontname','Helvetica');
% 
% set(gcf, 'Color', 'white'); % white bckgr
% export_fig -painters -r600 -q101 singleplotER.pdf
set(gcf, 'Color', 'None'); % white bckgr
set(gca, 'Color', 'None'); % white bckgr
export_fig( gcf, 'ERF-singleplotER','-transparent', ...
        '-painters','-pdf', '-r250' ); 

%%
% dat = tlckCmb{3,1};
% dat.avg = tlckCmb{3,1}.avg./1e-12.*10;
figure
cfg                     = [];
cfg.colorbar = 'yes';
cfg.linewidth = 2;
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.03];
% cfg.zlim = [0 75];
% cfg.graphcolor = clrOp1;
cfg.shading                 = 'flat';
cfg.maskstyle               = 'saturation';	
% cfg.gridscale  = 200;
% ft_topoplotER(cfg, dat)
ft_topoplotER(cfg, tlckCmb{1:10})

cfg.colormap = 'jet';


%%
hcb = colorbar;
set(hcb,'YTick',[0:25:75])
set(hcb,'fontsize',15);


% set(gcf, 'Color', 'None'); % white bckgr
% export_fig -painters -r600 -q101 topoN20.pdf
export_fig( gcf, 'ERF-topoplotER','-transparent', ...
        '-painters','-pdf', '-r250' ); 

%% ERF amplitude- Mean=/- SEM
% Load: maxPosN20Amp-PD, ControlN20mAmp

clrOp1 = rgb('dodgerBlue');%[0.7 0.2 0.2];
clrOp2 = rgb('DimGray');%[0.4 0.4 0.4];

load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20Amp.mat')
load('J:\MEG_Research\SEF\SEFepresults\Control\ControlN20mAmp.mat')
figure
contS = [1:10];
subjP = [1:10];
time_vector = [1 2 3 4 5 6 8];
PDERF       = maxPosN20Amp./1e-12.*10;
ControlERF  = ControlN20mAmp./1e-12.*10;

e2 = errorbar(time_vector, mean(ControlERF(contS,:),1),  sem(ControlERF(contS,:)), '-s')
set(e2,'LineWidth', 4.5 )
set(e2, 'MarkerSize', 18,'Color', clrOp2 , 'MarkerFaceColor', clrOp2 , ...
    'MarkerEdgeColor', clrOp2);
hold on
e1 = errorbar(time_vector, mean(PDERF(subjP,:),1),  sem(PDERF(subjP,:)), '-o');
set(e1,'LineWidth', 4.5 )
set(e1, 'MarkerSize', 18,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , ...
    'MarkerEdgeColor', clrOp1);


ylim([0 150])
xlim([0.5 8.5])

set(gca,'XTick',[1 2 3 4 5 6 8])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
% 
set(gca,'XTickLabel',...
    {'Stim ON', 'StimOFF(0min)', 'StimOFF(30min)', 'StimOFF(60min)', 'StimOFF(90min)'...
    'StimOFF(120min)', '','MedON'},'fontsize', 10, 'fontname', 'Georgia') ;
xTick = [1:6 8];
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, {'DBS OFF';'(30min)'},{'DBS OFF','(60min)'},...
    {'DBS OFF','(90min)'},{'DBS OFF','(120min)'}, 'Med ON'};
for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.04*(yTick(end)-yTick(1)-0.5),xTickLabel{k},...
        'HorizontalAlignment','center', 'Fontsize', 14, 'fontname', 'Helvetica',...
        'FontWeight','bold');
end
% 
% 
hxlabel = xlabel('Conditions')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Helvetica');
% 
hylabel = ylabel(['Amplitude fT-cm ^{-1}'])
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Helvetica');
% 
set(gca,'YTick',[0:50:150])
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Helvetica');
% 
htitle = title(['Amplitude of ',title_str]);
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');


% axis([0.5 8.5 0.5 1]);
hlegend = legend({'Control','PD'});
% hlegend = legend({'UPDRS III'},'Interpreter','latex');
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'Location','Northeast');
legend boxoff
box off

% set(0,'defaultAxesFontName', 'Helvetica')
set(gcf, 'Color', 'None'); % white bckgr
set(gca, 'Color', 'None'); % white bckgr

% saveas(gcf,'ERFAmp_SEMsavas.pdf')
% export_fig -painters -r600 -q101 ERFAmp-SEM.pdf
export_fig( gcf, 'ERF-trans','-transparent', ...
        '-painters','-pdf', '-r250' ); 

%% UPDRS- Mean=/- SEM
% Load: gammaPD, gammaControlmean
cd('J:\MEG_Research\SEF\sefgammaStatistics\')
load updrs.mat; 
figure
contS = [1:10];
subjP = [1:10];
time_vector = [0.7 2.3 4.3 7.7];
clrOpUpdrs = rgb('dodgerBlue');
e1=errorbar(time_vector, mean(updrs(subjP,[1 2 3 7]),1), ...
    sem(updrs(subjP,[1 2 3 7])), '-o','Color', clrOpUpdrs );
set(e1, 'MarkerSize', 18,'Color', clrOpUpdrs , 'MarkerFaceColor', clrOpUpdrs , ...
    'MarkerEdgeColor', clrOpUpdrs);
set(e1,'LineWidth', 4.5 )
% hold on
% e2=errorbar(time_vector, mean(gammaControlmean(contS,:),1),  sem(gammaControlmean(contS,:)), '-b.')
% set(e2,'LineWidth', 2.5 ) 
xlim([0.5 8.5])
ylim([0 40])

set(gca,'XTick',[1 2 3 4 5 6 8])
set(gca,'XTickLabel',[])
set(gca,'YTick',[0:10:40])
set(gca,'YTickLabel',[])
% set(gca,'XTickLabel',...
%     {'Stim ON', 'StimOFF(0min)', 'StimOFF(30min)', 'StimOFF(60min)', 'StimOFF(90min)'...
%     'StimOFF(120min)', '','MedON'},'fontsize', 10, 'fontname', 'Georgia') ;
% xTick = time_vector;
% set(gca,'xtick',xTick);
% yTick = get(gca,'ytick');
% set(gca,'xticklabel',[]);
% xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, {'DBS OFF';'(60min)'},...
%     'Med ON'};
% for k = 1:length(xTick)
%     text(xTick(k),yTick(1)-0.04*(yTick(end)-yTick(1)+0.5),xTickLabel{k},...
%         'HorizontalAlignment','center', 'Fontsize', 14, 'fontname', 'Helvetica',...
%         'FontWeight','bold');
% end


% hxlabel = xlabel('Conditions')
% set(hxlabel,'Fontsize',24);
% set(hxlabel,'FontWeight','bold');
% set(hxlabel,'Fontname','Helvetica');
% 
% hylabel = ylabel('UPDRS III  (Mean \pm SEM)')
% set(hylabel,'Fontsize',24);
% set(hylabel,'FontWeight','bold');
% set(hylabel,'Fontname','Helvetica');
% 
% set(gca,'YTick',[0:10:40])
% set(gca,'Fontsize', 18)
% set(gca,'Fontweight','bold');
% set(gca,'Fontname','Helvetica');
% 
% htitle = title('Clinical motor state');
% set(htitle,'Fontsize',24);
% set(htitle,'FontWeight','bold');
% set(htitle,'Fontname','Helvetica');


% axis([0.5 8.5 0.5 1]);
% hlegend = legend('PD (n=10)');
% % hlegend = legend({'UPDRS III'},'Interpreter','latex');
% set(hlegend,'Fontsize',18);
% set(hlegend,'Fontname','Helvetica');
% set(hlegend,'Location','Northeast');
% legend boxoff

set(0,'defaultAxesFontName', 'Calibri')
set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
% export_fig -painters -r600 -q101 UPDRS-SEM.pdf
box off
export_fig( gcf, 'UPDRS_new','-transparent', ...
        '-painters','-pdf', '-r250' ); 
% saveas(gcf, 'UPDRS-trans.pdf', 'pdf') 

%% Calc gamma from gammaControl_trials

for subj = 1:10
    for cond = 1:7
   
        gammaPD(subj,cond) = mean(gammaPD_trials{subj,cond}.powspctrm);
        gammaControl(subj,cond) = mean(gammaC_trials{subj,cond}.powspctrm);
    
    end
end
% a
%%



%% Calc GA gamma
% Control
for subj = 1:10
    for cond = 1:7
   
%         topoGammaPD{subj,cond}.dimord = 'chan_freq';
%         topoGammaPD{subj,cond}.powspctrm = topoGammaPD{subj,cond}.powspctrm';
% 
%         
        topoGammaPD{subj,cond}.dimord = 'chan_freq';
        topoGammaPD {subj,cond}.powspctrm = topoGammaPD{subj,cond}.powspctrm';
        
    end
end
%%
cfg = [];
cfg.parameter = 'powspctrm';
for cond = 1:7
    cond
%       GAgammaPD{cond}      = ft_freqgrandaverage(cfg, dat{:,cond});
%        GAgammaPD{cond}      = ft_freqgrandaverage(cfg, topoGammaPD{:,cond});
       GAgammaPD{cond} = ft_freqgrandaverage(cfg, dat{[1:7 9 10],cond});
        
end

%% Plots

% Topoplots PD
cfg                     = [];
cfg.colorbar            = 'no';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.03];
cfg.zlim                = [100 300]; 
cfg.shading             = 'interp';
% cfg.colormap = colormap('gray');
% cfg.maskstyle           = 'saturation';	
cfg.contournum          = 4;
% cfg.fontsize            = 20;

for loop = 1:1
    dat = GAgammaControl{loop};
    dat.powspctrm = GAgammaControl{loop}.powspctrm.*100;
    
    hFig = figure(1)
    subplot(3,3,loop)
    set(hFig, 'Position', [10 80 968 968]);
    ft_topoplotER(cfg, dat)
    title(['PD',num2str(loop)])
    
%     hcb = colorbar;
%     set(hcb,'YTick',[0:50:300])
%     set(hcb,'fontsize',15);
    
%     set(gcf, 'Color', 'w'); % white bckgr
%     export_fig( gcf, ['PDTopoGammaGA', num2str(loop)],...
%         '-painters','-png', '-r250' ); 
 
   
end
 export_fig( gcf, ['PDTopoGammaGA ', num2str(cfg.zlim)],...
        '-painters','-png', '-r250' ); 
%
% Topoplots Control
% cfg                     = [];
% cfg.colorbar            = 'no';
% cfg.layout              = 'neuromag306cmb.lay';
% cfg.xlim                = [0.02 0.03];
% cfg.zlim                = [0 150]; 
% cfg.shading             = 'interp';
% cfg.maskstyle           = 'saturation';	
% cfg.contournum          = 12;
% cfg.fontsize            = 20;

for loop = 1:7
    dat = GAgammaControl{loop};
    dat.powspctrm = GAgammaControl{loop}.powspctrm.*100;
    
    hFig = figure(8)
    subplot(3,3,loop)
    set(hFig, 'Position', [900 80 968 968]);
    ft_topoplotER(cfg, dat)
    title(['C',num2str(loop)])
    
%     hcb = colorbar;
%     set(hcb,'YTick',[0:50:300])
%     set(hcb,'fontsize',15);
    
    set(gcf, 'Color', 'w'); % white bckgr
%     export_fig( gcf, ['ControlTopoGammaGA', num2str(loop)],...
%         '-painters','-png', '-r250' ); 

end
    export_fig( gcf, ['ControlTopoGammaGA ', num2str(cfg.zlim)],...
        '-painters','-png', '-r250' ); 
    
%% Plot for paper Gamma- Mean=/- SEM
% Load: gammaPD, gammaControlmean
% cd('J:\MEG_Research\SEF\sefgammaStatistics\')
% load gammaControl.mat; load gammaPD.mat
figure
contS = [1:10];
subjP = [1:10];
time_vector = [1 2 3 4 5 6 8];
PDnorm       = gammaPD.*100;
Controlnorm  = gammaC.*100;
clrOp1 = rgb('dodgerBlue');%[0.7 0.2 0.2];
clrOp2 = rgb('DimGray');%[0.4 0.4 0.4];


% c = get(e2,'Children')
% yd = get(c(2),'Ydata');
% yd([7:9:end 8:9:end])=nan;
% yd(2:9:end)=mean(Controlnorm(contS,:),1);
% set(c(2),'YData',yd)
hold on

e1=errorbar(time_vector, mean(PDnorm(subjP,:),1), sem(PDnorm(subjP,:)), '-o');
set(e1, 'MarkerSize', 18,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , ...
    'MarkerEdgeColor', clrOp1);
set(e1,'LineWidth', 4.5 )

e2=errorbar(time_vector, mean(Controlnorm(contS,:),1),sem(Controlnorm(contS,:)),   '-s')
set(e2, 'MarkerSize', 18,'Color', clrOp2 , 'MarkerFaceColor', clrOp2 , ...
    'MarkerEdgeColor', clrOp2);
set(e2,'LineWidth', 4.5 )

xlim([0.5 8.5])
%ylim([190 300])

set(gca,'XTick',[1 2 4 5 6 8])
set(gca,'YTick',[175:25:400])
set(gca,'YTickLabel',[175:25:400])
% set(gca,'XTickLabel',...
%     {'Stim ON', 'StimOFF(0min)', 'StimOFF(30min)', 'StimOFF(60min)', 'StimOFF(90min)'...
%     'StimOFF(120min)', '','MedON'},'fontsize', 10, 'fontname', 'Georgia') ;
% xTick = [1:6 8];
% set(gca,'xtick',xTick);
% yTick = get(gca,'ytick');
% set(gca,'xticklabel',[]);
% xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, {'DBS OFF';'(30min)'},{'DBS OFF','(60min)'},...
%     {'DBS OFF','(90min)'},{'DBS OFF','(120min)'}, 'Med ON'};
% for k = 1:length(xTick)
%     text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)-0.05),xTickLabel{k},...
%         'HorizontalAlignment','center', 'Fontsize', 18, 'fontname', 'Helvetica',...
%         'FontWeight','bold');
% end


% hxlabel = xlabel('Conditions')
% set(hxlabel,'Fontsize',24);
% set(hxlabel,'FontWeight','bold');
% set(hxlabel,'Fontname','Helvetica');
% 
% hylabel = ylabel('Gamma power % increase   (Mean \pm SEM)')
% set(hylabel,'Fontsize',24);
% set(hylabel,'FontWeight','bold');
% set(hylabel,'Fontname','Helvetica');
% 
% set(gca,'YTick',[0:10:300])
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Helvetica');
% 
% htitle = title('Early somatosensory processing');
% set(htitle,'Fontsize',24);
% set(htitle,'FontWeight','bold');
% set(htitle,'Fontname','Helvetica');


% axis([0.5 8.5 0.5 1]);
hlegend = legend({'PD','Control'});
% hlegend = legend({'UPDRS III'},'Interpreter','latex');
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'Location','Northeast');
legend boxoff

% set(0,'defaultAxesFontName', 'Calibri')
% set(gcf, 'Color', 'None'); % white bckgr
% export_fig -painters -r600 -q101 Gamma-SEM_new.pdf
% export_fig( gcf, 'Gamma-SEM_new.pdf','-transparent', ...
%         '-painters','-pdf', '-r250' ); 
set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'Gamma-magbad-vs-xscan',...
        '-painters','-pdf', '-r250' ); 
% saveas(gcf, 'gamma_try.emf')    

%% Gamma - multiplotER, singleplotER, topoplot

cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relchange';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);

cfg                         = [];
cfg.jackknife               = 'yes';
gammaDesc                   = ft_freqdescriptives(cfg, gammaFreq_bsl);

cfg                         = [];
gammaCmb                    = ft_combineplanar(cfg, gammaDesc);

%%
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.linewidth = 2.2;
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0 0.2];
cfg.ylim = [0 2.5];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.graphcolor = clrOp1;

ft_multiplotER(cfg, gammaCmb)

set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'gammaMultiPlot.pdf','-transparent', ...
        '-painters','-pdf', '-r250' ); % export_fig( gcf, 'multiplotGamma',...
%         '-painters','-png', '-r250' ); 

%%    
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.channel = {'MEG1812+1813'};
cfg.colorbar = 'yes';
cfg.linewidth = 4;
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0 0.2];
cfg.ylim = [0 2];
cfg.graphcolor = clrOp1;
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
ft_singleplotER(cfg, gammaCmb)

%%
figure
% dat      = gammaDesc.powspctrm([133, 134],:,:,:);
% dat_sem  = gammaDesc.powspctrmsem([133, 134],:,:,:);

chan = match_str(gammaCmb.label,'MEG0442+0443'); 
toi = find(gammaCmb.time>=-0.01 & gammaCmb.time<=0.2);
% time_vector = 
e2 = boundedline(gammaCmb.time(toi), squeeze(gammaCmb.powspctrm(chan,:,toi)).*100,...
    squeeze(gammaCmb.powspctrmsem(chan,:,toi).*100), 'cmap',  clrOp1)
% set(e2, 'MarkerSize', 10,'Color', clrOp2 , 'MarkerFaceColor', clrOp1 , ...
%     'MarkerEdgeColor', clrOp2);
set(e2,'LineWidth', 2.5)
ylim([0 210])

set(gca,'XTick',[0:0.05:0.2])
set(gca,'XTickLabel',[0:0.05:0.2])
set(gca,'Fontsize', 18, 'fontname', 'Helvetica','FontWeight','bold');

set(gca,'YTick',[0:50:200])
set(gca,'YTickLabel',[0:50:200])
set(gca,'Fontsize', 18, 'fontname', 'Helvetica','FontWeight','bold');

hylabel = ylabel('Gamma power % increase')
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Helvetica');

hxlabel = xlabel('Time s')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Helvetica');

htitle = title('Mean of selected channels');
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');

set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'gamma_single_new.pdf','-transparent', ...
        '-painters','-pdf', '-r250' ); 
% export_fig -painters -r600 -q101 singleplotGamma.pdf
% export_fig( gcf, 'singleGamma',...
%         '-painters','-png', '-r250' ); 

%%
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.linewidth = 2;
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.zlim = [0 2];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
ft_topoplotER(cfg, gammaCmb)
hcb = colorbar;
set(hcb,'YTick',[0:0.5:2])
set(hcb,'fontsize',15);


set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'gamma_topo_new.pdf','-transparent', ...
        '-painters','-pdf'); 
% export_fig -painters -r600 -q101 topoplotGamma.pdf
% export_fig( gcf, 'topoGamma',...
%         '-painters','-png', '-r250' ); 

%% Trend plots - Gamma
% 6-7
figure
PDnorm       = gammaPD.*100;
Controlnorm  = gammaControl.*100;
conditions_str = {'DBS ON','DBS OFF-0min','DBS OFF-30min)','DBS OFF(60min)',...
    'DBS OFF(90min)','DBS OFF(120min)','Med ON'};
cond1= 1; cond2 = 7;
p1 = yyplot([1,2],PDnorm(:,[cond1 cond2])')
set(p1,'MarkerSize', 10,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , 'LineWidth', 2,...
    'MarkerEdgeColor', clrOp1);
% alpha(p1,.5)

bar([1.5],mean([PDnorm(:,cond2)-PDnorm(:,cond1)]))
errorbar(1.5, mean([PDnorm(:,cond2)-PDnorm(:,cond1)])+sem([PDnorm(:,cond2)-PDnorm(:,cond1)]), ...
    sem([PDnorm(:,cond2)-PDnorm(:,cond1)]))
xlim([0.5 2.5]); ylim([0 300])

set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'DBS ON', 'DBS OFF(0min)'})
set(gca,'Fontsize',18); set(gca,'FontWeight','bold');

set(gca,'YTick',[0:50:300])
set(gca,'Fontsize',18); set(gca,'FontWeight','bold');
set(gca,'Fontname','Calibri');

hylabel = ylabel('Gamma power % increase')
set(hylabel,'Fontsize',20);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Calibri');

htitle = title('Trend');
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Calibri');

%%

figure
[ax,b,p] = plotyy(1.5,mean([PDnorm(:,cond2)-PDnorm(:,cond1)]),[1 2],PDnorm(:,[cond1 cond2])','bar','plot');
hold on
errorbar(ax(1),1.5, mean([PDnorm(:,cond2)-PDnorm(:,cond1)])+sem([PDnorm(:,cond2)-PDnorm(:,cond1)]), ...
    sem([PDnorm(:,cond2)-PDnorm(:,cond1)]))
xlim(ax(2),[0.5 2.5])
ylim(ax(1), [0 60])


%%
diff12=[gammaPD(:,1)-gammaPD(:,2)];
diff17=[gammaPD(:,1)-gammaPD(:,7)];
diff67=[gammaPD(:,6)-gammaPD(:,7)];

%%
scatter(dbsMonths, diff12)
ylim([-1 0.4]); xlim([8 20])
lsline
[rho p]=corrcoef(dbsMonths, diff12);
disp(['p: ',num2str(p(1,2))])
disp(['rho: ',num2str(rho(1,2))])

% PDYrs-diff67: rho: -0.37 p: 0.29
% PDYrs-diff17: rho: -0.32 p: 0.37
% PDYrs-diff12: rho: -0.06 p: 0.87

% dbsMonths-diff67: rho:0.36 p:0.30
% dbsMonths-diff17: rho:0.24 p:0.51
% dbsMonths-diff12: rho:0.22 p:0.55

%% Between -  condition wise
p=[]; h=[]; 
clear stat;
% [p, table] = anova_rm(gammaPD);

for loop= 1:7

    dat = (gammaPD(:,loop) - gammaControlmean(:,loop));
    dispSummary(dat)
    [h(loop) p(loop) ci(loop,:) stat(loop)] = ttest(dat);
   
end

% p
% h
% abs([stat(:).tstat])
result_PDvsC= [h' p'];

logPD = log(gammaPD);
logC = log(gammaControlmean);
[p, table] = anova_rm({logPD,logC});

%% Does the degree of rebound predict the washout? - Yes
dat_1 = [gammaPD(:,2)-gammaPD(:,1)];
dat_2 = [gammaPD(:,2)-gammaPD(:,6)];%diseaseYrs;
scatter(dat_1, dat_2), lsline
axis([-0.4 1.2 -0.4 1.2])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['p: ',num2str(p)])
disp(['rho: ',num2str(rho)])
% p: 0.002
% rho: 0.838

%% 1. Does the degree of rebound correlate with worsening motor symptoms? - No
dat_1 = [gammaPD(:,2)-gammaPD(:,1)];
dat_2 = [updrs(:,2)-updrs(:,1)];%diseaseYrs;
% scatter(dat_1, dat_2), lsline
% axis([0 1 0 25])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['Sp-p: ',num2str(p)])
disp(['Sp-rho: ',num2str(rho)])

% Sp-p: 0.62287
% Sp-rho: 0.17793
% p: 0.36998
% rho: 0.31836

%% 2. DBS ON to untreated? UPDRS  vs gamma - NO
dat_1 = [gammaPD(:,6)-gammaPD(:,1)];
dat_2 = [updrs(:,6)-updrs(:,1)];
% scatter(dat_1, dat_2), lsline
% axis([0 1 0 25])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['Sp-p: ',num2str(p)])
disp(['Sp-rho: ',num2str(rho)])

% p: 0.11276
% rho: 0.53286
% Sp-p: 0.1854
% Sp-rho: 0.45593

%% 3. DBS ON to MED ON? UPDRS  vs gamma - No
dat_1 = [gammaPD(:,1)-gammaPD(:,7)];
dat_2 = [updrs(:,1)-updrs(:,7)];%diseaseYrs;
scatter(dat_1, dat_2), lsline
% axis([0 1 0 25])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['Sp-p: ',num2str(p)])
disp(['Sp-rho: ',num2str(rho)])

% p: 0.16816
% rho: 0.47224
% Sp-p: 0.41333
% Sp-rho: 0.29179

%% 4. DBS OFF(0min) to DBS OFF(120min)? UPDRS  vs gamma - No
dat_1 = [gammaPD(:,2)-gammaPD(:,7)];
dat_2 = [updrs(:,2)-updrs(:,7)];%diseaseYrs;
scatter(dat_1, dat_2), lsline
% axis([0 1 0 25])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['Sp-p: ',num2str(p)])
disp(['Sp-rho: ',num2str(rho)])

% p: 0.22946
% rho: 0.4179
% Sp-p: 0.39015
% Sp-rho: 0.30582

%% 5. DBS OFF-0 to MED ON? UPDRS  vs gamma - No
dat_1 = [gammaPD(:,2)-gammaPD(:,7)];
dat_2 = [updrs(:,2)-updrs(:,7)];%diseaseYrs;
scatter(dat_1, dat_2), lsline
% axis([0 1 0 25])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['Sp-p: ',num2str(p)])
disp(['Sp-rho: ',num2str(rho)])

% p: 0.22946
% rho: 0.4179
% Sp-p: 0.39015
% Sp-rho: 0.30582

%% 6. DBS OFF 120min to MED ON? UPDRS  vs gamma - No
dat_1 = [gammaPD(:,6)-gammaPD(:,7)];
dat_2 = [updrs(:,6)-updrs(:,7)];%diseaseYrs;
scatter(dat_1, dat_2), lsline
% axis([0 1 0 25])

[rho p]=corr(dat_1, dat_2,'type','Spearman');
disp(['Sp-p: ',num2str(p)])
disp(['Sp-rho: ',num2str(rho)])

% p: 0.81021
% rho: -0.087425
% Sp-p: 0.6096
% Sp-rho: -0.18464


%%
% Does the degree of rebound predict the washout? - Yes (p:0.002; rho:0.838)
% DBS ON to untreated? UPDRS  vs gamma - YES  (p:0.04, rho: 0.64)

% Does the degree of rebound correlate with worsening motor symptoms? - No
%                                                                (p:0.37 ; rho:0.32 )
% DBS OFF 120min to MED ON? UPDRS  vs gamma - No (p:0.42; rho: -0.28)




