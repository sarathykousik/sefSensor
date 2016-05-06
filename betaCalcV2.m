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

%%

cfg               = [];
cfg.latency       = [0.08 0.14];


cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 1:6
    disp(['######   ', num2str(subj)])

    for cond = 1:6
        load(char(file_sub(subj, cond)));
        
        cfg = [];
        cfg.baseline                = [-0.100 -0.040];
        cfg.baselinetype            = 'relchange';
        cfg.parameter               = 'powspctrm';
        gammaFreq_bsl               = ft_freqbaseline(cfg,gammaFreq);
        
        cfg                         = [];
        gammaFreq_bsl               = ft_combineplanar(cfg, gammaFreq_bsl);
        
        cfg             = [];
        cfg.jackknife   = 'yes';
        gammaDesc       = ft_freqdescriptives(cfg, gammaFreq_bsl);
  
        topoGamma         = ft_selectdata(cfg, gammaDesc);
        maxGamma(subj,cond,:)     = max(topoGamma.powspctrm,[],3);

        maxGammaStruct{subj,cond} = ...
                       ft_selectdata(cfg_sing, gammaDesc);

        maxGammaStruct{subj,cond}.powspctrm = max(topoGamma.powspctrm,[],3);
        maxGammaStruct{subj,cond}.name= char(file_sub(subj, cond));
                   
%         gammaValCh(subj,cond,:) = max(squeeze(data.gammaFreq_bsl.powspctrm(:,:,toi)),[],2);
        
        allGamma{subj, cond} = gammaDesc;
        
        
%         cfg                     = [];
%         cfg.parameter = 'powspctrm';
%         cfg.colorbar = 'yes';
%         cfg.layout              = 'neuromag306cmb.lay';
%         cfg.xlim                = [0.02 0.06];
%         cfg.shading                 = 'interp';
%         cfg.maskstyle               = 'saturation';	
%               figure(subj),
%         subplot(2,3,cond),ft_topoplotER(cfg, data.gammaFreq_bsl)
%         title(num2str(cond));

    end
end

blankForPlot = maxGammaStruct{1,1};
blankForPlot.powspctrm = zeros(204,1);

chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

chanSel_2                 = {'MEG0432+0433', 'MEG0442+0443'};

for loop = 1:length(chanSel_1) 
    chanPos(loop) =  find(strcmp(blankForPlot.label,chanSel_1(loop)));
end

cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.06];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';

for subj = 1:6
     maxVal = max(maxGamma(subj, :,chanPos), [], 3);
    for cond = 1:6
        maxPos(subj,cond) = find(maxGamma(subj, cond,:)==maxVal(cond));
%     end

    cfg.highlightchannel =  {maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}};
    cfg.highlight        =  'numbers';
    figure(subj),
    subplot(2,3,cond),ft_topoplotER(cfg, maxGammaStruct{subj,cond})
    title({maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}})
   end
end
 
for subj = 1:6
    for cond = 1:6
        maxGammaChExtract(subj, cond) = max(maxGamma(subj,cond,maxPos(subj,cond)));
    end
end
figure,boxplot(maxGammaChExtract([1,4:6],:))

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
maxGammaValNorm = (maxGammaChExtract-baseline);

%% Plots

figure,
subplot 221
plot(maxGammaValNorm([1:4,6],:)', '-o'), legend('009','010', '011', '012', '014')
set(gca,'XTickLabel',...
    {'Stim ON(30)', 'Stim OFF 0(60)', 'Stim OFF 60(120)', 'Stim OFF 90(150)','Stim OFF 120(180)','Med ON(210)'}, ...
        'XTick', 1:6)
title('Norm')
xlabel('Timeline'), ylabel('Normalized Gamma power')

subplot 222
plot(maxGammaVal([1:4,6],:)', '-o'), legend('009','010', '011', '012',  '014')
set(gca,'XTickLabel',...
    {'Stim ON(30)', 'Stim OFF 0(60)', 'Stim OFF 60(120)', 'Stim OFF 90(150)','Stim OFF 120(180)','Med ON(210)'}, ...
        'XTick', 1:6)
title('No Norm')
xlabel('Timeline'), ylabel('Baseline-corrected Gamma power')

subplot 223
boxplot(maxGammaValNorm([1:4,6],:)), title('Norm')
xlabel('Timeline'), ylabel('Baseline-corrected Gamma power')

subplot 224
boxplot(maxGammaVal([1:4,6],:)), title('No Norm')
xlabel('Timeline'), ylabel('Baseline-corrected Gamma power')

suptitle('Change in Gamma power')
figure
plot(maxGammaVal(:,2:6)', 'o'), lsline
legend('009','010', '011', '012', '013', '014')
set(gca,'XTickLabel',...
    { 'Stim OFF 0(60)', 'Stim OFF 60(120)', 'Stim OFF 90(150)','Stim OFF 120(180)','Med ON(210)'}, ...
        'XTick', 2:6)
title('No Norm')
xlabel('Timeline'), ylabel('Baseline-corrected Gamma power')

%% t-test between groups

for loop = 2:6
    [h(loop-1), p(loop-1), ci, stats(loop-1)] = ...
        vartest2(maxGammaVal([1:2,3],1), maxGammaVal([1:2,3],loop),0.05);
end

disp(h)
disp(p)

%%

[p, table, stats] = anova1(maxGammaValNorm([1 2 3 4 5 6],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.05,'ctype','lsd');

%%
figure,
subj = [1:6];
errorbar(1:6,mean(maxGammaValNorm(subj,:),1), min(maxGammaValNorm(subj,:)),...
    max(maxGammaValNorm(subj,:)),'-ks', 'LineWidth',1.5,'MarkerSize',10);
hold on
plot(maxGammaValNorm(subj,:)', 'bx', 'MarkerSize',6, 'LineWidth',1.5)

% boxplot(maxGammaVal)
set(gca,'XTickLabel',...
    {'Stim ON(30)', 'Stim OFF 0(60)', 'Stim OFF 60(120)', 'Stim OFF 90(150)','Stim OFF 120(180)','Med ON(210)'},...
    'fontsize', 9, 'fontname', 'Georgia') ;

% axis([0.5 6.5 -1e-12 2.5e-12]);
    
% hlegend = legend({'Mean $\pm$ SD'},'Interpreter','latex');
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

% groups={[1,2],[1,3],[1,4], [1,5], [1,6]};
% H = sigstar(groups,[0.01,0.01,0.01,0.01,0.01]);
% 
% htext = text(6.08,1.459e-22,' **   p<0.01',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(6.08,1.459e-22,{'$\pm$'},'Interpreter','latex',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% htext = text(6.08,1.459e-22,'SD',...
%    	'HorizontalAlignment','Center',...
%    	'BackGroundColor','none','Fontsize',14,...
%     'Fontangle','italic', 'Fontname','Georgia');
% box off

