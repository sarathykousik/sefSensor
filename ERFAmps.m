

chan_sel = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

%%
toi_1=[0.02 0.03];
toi_2=[0.032 0.042];
toi_3=[0.043 0.05];
toi_4=[0.059 0.068];
toi_5=[0.07 0.08];

for subj = 1:20
    disp(['######   ', num2str(subj)])

    for cond = 1:7
        chanPos = match_str(allTlck{subj,cond}.label, chan_sel);
        [maxERFChExtract_1(subj, cond) maxERFPos_1(subj, cond)]= ...
            max(max(allTlck{subj,cond}.avg(chanPos,find(allTlck{subj,cond}.time>=toi_1(1) & allTlck{subj,cond}.time<=toi_1(2)))));
        
        [maxERFChExtract_2(subj, cond) maxERFPos_2(subj, cond)]= ...
            max(max(allTlck{subj,cond}.avg(chanPos,find(allTlck{subj,cond}.time>=toi_2(1) & allTlck{subj,cond}.time<=toi_2(2)))));
        
        [maxERFChExtract_3(subj, cond) maxERFPos_3(subj, cond)]= ...
            max(max(allTlck{subj,cond}.avg(chanPos,find(allTlck{subj,cond}.time>=toi_3(1) & allTlck{subj,cond}.time<=toi_3(2)))));
        
        [maxERFChExtract_4(subj, cond) maxERFPos_4(subj, cond)]= ...
            max(max(allTlck{subj,cond}.avg(chanPos,find(allTlck{subj,cond}.time>=toi_4(1) & allTlck{subj,cond}.time<=toi_4(2)))));
        
        [maxERFChExtract_5(subj, cond) maxERFPos_5(subj, cond)]= ...
            max(max(allTlck{subj,cond}.avg(chanPos,find(allTlck{subj,cond}.time>=toi_5(1) & allTlck{subj,cond}.time<=toi_5(2)))));
        
    end
end
%%
for subj = 1:10
    disp(['######   ', num2str(subj)])

    for cond = 1:7

    evokedgammaC(subj,cond)= evokedgammaC_trials{subj,cond}.powspctrm;
    evokedgammaPD(subj,cond)= evokedgammaPD_trials{subj,cond}.powspctrm;

    end
end

%%

% subjs = 1:10;
clf(gcf)
clrOp1 = rgb('dodgerBlue');%[0.7 0.2 0.2];
clrOp2 = rgb('DimGray');%[0.4 0.4 0.4];
% plot_param = amp20./1e-12.*10;
plot_param = [evokedgammaC; evokedgammaPD]; %amp20./1e-12.*10; %
% plot_param_2 = [indgammaC; indgammaPD].*100;%[ratioC;ratioPD].*100;% %amp40./1e-12.*10;
% plot_param_3 = amp70./1e-12.*10;

title_str='Evoked gamma'
% load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20Amp.mat')
% load('J:\MEG_Research\SEF\SEFepresults\Control\ControlN20mAmp.mat')
% figure
contS = [1:10];
subjP = [11:20];
time_vector = [1 2 3 4 5 6 8];
% PDERF       = maxPosN20Amp./1e-12.*10;
% ControlERF  = ControlN20mAmp./1e-12.*10;

e2 = errorbar(time_vector, mean(plot_param(contS,:),1),  sem(plot_param(contS,:)), '-s')
set(e2,'LineWidth', 4.5 )
set(e2, 'MarkerSize', 18,'Color', clrOp2 , 'MarkerFaceColor', clrOp2 , ...
    'MarkerEdgeColor', clrOp2);
hold on
e1 = errorbar(time_vector, mean(plot_param(subjP,:),1),  sem(plot_param(subjP,:)), '-o');
set(e1,'LineWidth', 4.5 )
set(e1, 'MarkerSize', 18,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , ...
    'MarkerEdgeColor', clrOp1);

% e4 = errorbar(time_vector, mean(plot_param_2(contS,:),1),  sem(plot_param_2(contS,:)), '-.s')
% set(e4,'LineWidth', 4.5 )
% set(e4, 'MarkerSize', 18,'Color', clrOp2 , 'MarkerFaceColor', clrOp2 , ...
%     'MarkerEdgeColor', clrOp2);
% hold on
% e3 = errorbar(time_vector, mean(plot_param_2(subjP,:),1),  sem(plot_param_2(subjP,:)), '-.o');
% set(e3,'LineWidth', 4.5 )
% set(e3, 'MarkerSize', 18,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , ...
%     'MarkerEdgeColor', clrOp1);

% e5 = errorbar(time_vector, mean(plot_param_3(contS,:),1),  sem(plot_param_3(contS,:)), '-.s')
% set(e5,'LineWidth', 4.5 )
% set(e5, 'MarkerSize', 18,'Color', clrOp2 , 'MarkerFaceColor', clrOp2 , ...
%     'MarkerEdgeColor', clrOp2);
% hold on
% e6 = errorbar(time_vector, mean(plot_param_3(subjP,:),1),  sem(plot_param_3(subjP,:)), '-.o');
% set(e6,'LineWidth', 4.5 )
% set(e6, 'MarkerSize', 18,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , ...
%     'MarkerEdgeColor', clrOp1);

% ylim([0.9 1.3])
% ylim([89 101])

xlim([0.5 8.5])

set(gca,'XTick',[1 2 3 4 5 6 8])
set(gca,'XTickLabel',[])

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
    text(xTick(k),yTick(1)-.1*(yTick(end)-yTick(1)-0.2),xTickLabel{k},...
        'HorizontalAlignment','center', 'Fontsize', 14, 'fontname', 'Helvetica',...
        'FontWeight','bold');
end


hxlabel = xlabel('Conditions')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Helvetica');

% hylabel = ylabel(['Relative gamma power % increase ','(Mean ',177, 'SEM)'])
% hylabel = ylabel(['Amplitude fT-cm ^{-1}']);
hylabel = ylabel(['Ratio of induced to total Gamma power (%)']);
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Helvetica');

set(gca,'YTick',[0:50:400])
% set(gca,'YTick',unique(sort([[1.2],[0:50:150]])))
% set(gca,'YTick',[0:0.1:2])
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Helvetica');
% 
htitle = title(title_str);
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');


% axis([0.5 8.5 0.5 1]);
hlegend = legend({'Control','PD'})
% hlegend = legend({'Control-total','PD-total','Control-induced','PD-induced'});
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
export_fig( gcf, title_str,'-transparent', ...
        '-painters','-pdf', '-r250' ); 

%% global mean field power

% tlckAllDat = [tlck_PD; tlck_control];

cfg = [];
cfg.method ='amplitude';

for subjLoop=1:20
%     for condLoop=
    gmf{subjLoop} = ft_globalmeanfield(cfg, grandAvgSubj{subjLoop});
    
end

%%

% toi_n20 = [20,30];
% toi_n40 = [40,60];
% Alltlck - 1:10 - PD;  11-20-Ctrl
for subj=1:20
    for cond=1:7

        toiWin_20  = [toi_PD_C(subj,1)+499:toi_PD_C(subj,1)+503];
        toiWin_40  = [toi_PD_C(subj,2)+499:toi_PD_C(subj,2)+503];
        toiWin_70  = [toi_PD_C(subj,3)+499:toi_PD_C(subj,3)+503];

        avgERF     = mean(allTlck{subj,cond}.avg,1);
        
        amp20(subj,cond) = mean(avgERF(toiWin_20));
        amp40(subj,cond) = mean(avgERF(toiWin_40));
        amp70(subj,cond) = mean(avgERF(toiWin_70));

   
   
    end
end











