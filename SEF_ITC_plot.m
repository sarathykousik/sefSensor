clear all;   clc;    close all;

% cd('J:\MEG_Research\SEF\SEFGamma3tap')
filenames = dir('*.mat');

% for loop = 1: length(filenames)
%     
%     file_sub(loop) = {filenames(loop).name};
% 
% end
% 
% file_sub = reshape(file_sub, 6, []);
% file_sub = file_sub';
% [row col] = size(file_sub);

filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
    cfg              = [];
    cfg.keeptrials   = 'yes';
    cfg.output       = 'fourier';
    cfg.channel      ='MEGGRAD';
    % cfg.trials       = 1:400;
    cfg.method       = 'wavelet';
    cfg.pad          = 2;
    cfg.foi          = 6:2:100;                          
    cfg.width        = 5; 
    cfg.toi          = -0.5:0.010:0.5;        
    TFRfourier    = ft_freqanalysis(cfg, visClean);

    % cfg             = [];
    % cfg.variance   = 'yes';
    % TFRDesc       = ft_freqdescriptives(cfg, TFRfourier);

    tmpdat = TFRfourier.fourierspctrm;
    tmpdat = tmpdat./(abs(tmpdat)); % this will normalise each trial for its
    ITC = abs(mean((tmpdat))); % this will give the itc

    % TFRfourier = rmfield(TFRfourier, 'fourierspctrm');
    TFRfourier = rmfield(TFRfourier, {'fourierspctrm','trialinfo', 'cumtapcnt'});
    TFRfourier.dimord = 'chan_freq_time';
    TFRfourier.powspctrm = squeeze(ITC);
    save (['J:\MEG_Research\SEF\SEF_ITC\',b, '-ITC'], 'TFRfourier', '-v7.3');

end

%%

filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 6, []);
file_sub = file_sub';
[row col] = size(file_sub);


for subj = 1:row
    disp(['######   ', num2str(subj)])

    for cond = 1:col
        
      load(char(file_sub(subj, cond)));        
      SEFITC_PD{subj, cond} = TFRfourier;
      SEFITC_PD{subj, cond}.filename = char(file_sub(subj, cond));
      
    end
end

%% control


% roi_1 = {'MEG0222+0223', 'MEG0412+0413', 'MEG0422+0423', 'MEG0442+0443', ...
%             'MEG0632+0633', 'MEG0642+0643', 'MEG0712+0713', 'MEG0222+0223'};
% roi_2 = {'MEG0432+0433', 'MEG0742+0743', 'MEG1812+1813', 'MEG1822+1823', ...
%             'MEG1832+1833', 'MEG1842+1843'};

conITC_60_100  =[];
PDITC_60_100  =[]

cfg_toi               = [];
cfg_toi.latency       = [0.02 0.08];
cfg_toi.foilim        = [55 95];


cfg_sing               = [];
cfg_sing.latency       = [0.02];

cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'absolute';
cfg.parameter               = 'powspctrm';


cfg_cmb = [];
% cfg_cmb.combinemethod = 'itc';
% ITC_bsl                 = ft_freqbaseline(cfg,...
%          ft_combineplanar_itc(cfg_cmb,  TFRfourier));


for subj = 1:10
    disp(['######   ', num2str(subj)])

    for cond = 1:6
        
         ITC_bsl                 = ft_freqbaseline(cfg,...
             ft_combineplanar_itc(cfg_cmb,  SEFITC_control{subj, cond}));

        topoGamma60_100         = ft_selectdata(cfg_toi, ITC_bsl);
       
        conITC_60_100(subj,cond,:)     = max(max(topoGamma60_100.powspctrm,[],2),[],3);

    end
end

% Now select the ROI and mean ITC in ROI
load('J:\MEG_Research\SEF\neighbours.mat')
% load('J:\MEG_Research\SEF\maxPosN20.mat')
load('J:\MEG_Research\SEF\SEFepresults\control\maxPosN20Control.mat')

for subj = 1:10
    for cond = 1:6

        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};

        conITC_N20Pos(subj,cond) = conITC_60_100(subj, cond, maxPosN20(subj,cond));
        
        chanPos = match_str(topoGamma60_100.label, chan_sel);        
        conITC_meanNeigh(subj,cond) = mean(conITC_60_100(subj, cond,chanPos),3);

    end
end

%%%%%%%%%%%%%%%%%%%%%% PD


for subj = 1:10
    disp(['######   ', num2str(subj)])

    for cond = 1:6
        
         ITC_bsl                 = ft_freqbaseline(cfg, ...
                        ft_combineplanar_itc(cfg_cmb, SEFITC_PD{subj, cond}));

        topoGamma60_100         = ft_selectdata(cfg_toi, ITC_bsl);
        
        
        PDITC_60_100(subj,cond,:)     = max(max(topoGamma60_100.powspctrm,[],2),[],3);

    end
end

% Now select the ROI and mean ITC in ROI
load('J:\MEG_Research\SEF\neighbours.mat')
load('J:\MEG_Research\SEF\SEFepresults\updated\maxPosN20.mat')

for subj = 1:10
    for cond = 1:6

        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};

        PDITC_N20Pos(subj,cond) = PDITC_60_100(subj, cond, maxPosN20(subj,cond));

        chanPos = match_str(topoGamma60_100.label, chan_sel);
        PDITC_meanNeigh(subj,cond) = mean(PDITC_60_100(subj, cond,chanPos),3);
        
    end
end

contS = [1:10];
subjP = [1:10];
[p, table] = anova_rm({gammaPD(subjP,:), gammaControlmean(contS,:)});
% [p, table, stats] = anova_rm({PDITC_N20Pos(subj,:),  conITC_N20Pos(contS,:)});

% h1= figure,
% clf(h1)
% contS = [1:8 9 10]
% subj = [1  3:10]
% errorbar(mean(PDITC_roi1(subj,:)), sem(PDITC_roi1(subj,:)), '-ro'),hold on
% errorbar(mean(conITC_roi1(contS,:)), sem(PDITC_roi1(subj,:)), '-bo')
% % errorbar(mean(PDITC_roi2(subj,:)), sem(PDITC_roi2(subj,:)), '--ro'),hold on
% % errorbar(mean(conITC_roi2(subj,:)),sem(conITC_roi2(subj,:)), '--bo')
% legend('PD ROI #1', 'Control ROI #1');%, 'PD ROI #2', 'Control ROI #2')
% conITC = conITC_N20Pos(contS, :)';
% PDITC  = PDITC_N20Pos(subj, :)';
% 
% [p, table, stats] = anova_rm({N20mAmpPD(subj,:), N20mAmpControl(contS,:)});

%
% contS = [1:8 9 10];
% subjP = [1  3:5 7:10 12];

%%
cond = [6 7];
[p, table] = anova_rm({gammaMeanExtract(:,cond), gammaCExtract(:,cond)});


%%
% clc
gammaPD_log = (gammaMeanExtract);
gammaControl_log = (gammaCExtract);
contS = [1:10];
subjP = [1:10];
p=[]; h=[]; clear stat;
combos = nchoosek([1:7],2);
% combos 
% 1     2
% 1     7
% 6     7
% 2     6
% 3     6
somevar = 1;
for loop= [1 6 21] 
%     disp(combos(loop,:))
    [h(somevar) p(somevar) ci stat(somevar)] = ttest2(...
        (gammaPD_log(subjP,combos(loop,1))- gammaPD_log(subjP,combos(loop,2))),...
        (gammaControl_log(contS,combos(loop,1))- gammaControl_log(contS,combos(loop,2))),0.05);
        somevar = somevar+1;
    
end
disp(p)
disp(h)
% abs([stat(:).tstat])
% [p, table] = anova_rm({ControlITPC_mtm(contS,:), PDITPC_mtm(subjP,:)});

% Final: 

% ==================================================================================
% Gamma-mean around N20m:  0.0084    0.0851    0.0116    0.1012    0.6642    0.0241
%                          2.9594    1.8221    2.8081    1.7277    0.4414    2.4628    
% Gamma - at N20m Pos:      0.0044    0.0967    0.0012    0.0676    0.6959    0.1314
% Gamma - max of N20 neigh: 0.0104    0.0987    0.0064    0.0599    0.9137    0.0872
% ITC max of N20 neigh:     0.0032    0.3944    0.0220    0.3164    0.5632    0.1319
% ITC at N20m Pos:          0.0037    0.4934    0.0088    0.1118    0.4428    0.0475

% Gamma - mean around N20:  0.0049    0.1326    0.0035    0.0420    0.9354    0.0463
%                           0.0032    0.1076    0.0027    0.0409    0.9480    0.0476
%  tstat:                   
% % PD:1,3-10;Cont:1-8,10:    *0.0021    0.1011    *0.0028    *0.0335    0.9780    *0.0543
% tstat:                       3.4271    1.6986     3.5049    2.2120     0.0663   2.1353
%  PD:1,3-5,7-10,12;contall    0.0061    0.1732    0.0103    0.0499    0.6731    0.0083
%  tstat:                      3.1314    1.4217    2.8840    2.1107    0.4293    2.9885 

% ITC :                     0.0099    0.3383    0.0053    0.1778    0.3832    0.0444
% tstat:                    2.9280    0.9871    3.2272  1.4095      0.8967   2.1814

% 54-94 Hz
% p:                        0.0347    0.4120    0.0095    0.2215    0.3913    0.0492
% tstat:                    2.2951    0.8410    2.9204    1.2692    0.8797    2.1184

%%
% contS = [1:10];
subjP = [1:10];
p=[]; h=[]; 
clear stat;
for loop= 1:21
% for loop= [1:15] 

    [h(loop) p(loop) ci stat(loop)] = ttest(...
            (gammaMeanExtract(subjP,combos(loop,1))- gammaMeanExtract(subjP,combos(loop,2))));

%         maxGammaChExtract(subjP,loop),maxGammaChExtract(subjP,loop+1),0.05);
        
    
end

p
h
abs([stat(:).tstat])
stat_ans =[combos(:,1)';combos(:,2)'; h; p; abs([stat(:).tstat])];

%%
% contS = [1:10];
subjP = [1:10];
p=[]; h=[]; 
clear stat;
for loop= 1:7
% for loop= [1:15] 

    [h(loop) p(loop) ci stat(loop)] = ttest(...
            [gammaMeanExtract(:,loop)- gammaCExtract(:,loop)]);

%         maxGammaChExtract(subjP,loop),maxGammaChExtract(subjP,loop+1),0.05);
        
    
end

p
h
abs([stat(:).tstat])
% stat_ans =[combos(:,1)';combos(:,2)'; h; p; abs([stat(:).tstat])];
%% UPDRS stat
sel_cols = [1 2 3 7];
[p, table, stats] = anova_rm(updrs(:,sel_cols))
figure
[c, m, h, nms] = multcompare(stats,'alpha',0.05,'ctype','hsd');

%
p=[]; h=[]; clear stat;
combos = nchoosek(sel_cols,2);
% combos 
% 1     2
% 1     7
% 6     7
% 2     6
% 3     6
somevar = 1;
for loop= [1:6] 
%     disp(combos(loop,:))
    [h(somevar) p(somevar) ci stat(somevar)] = ttest2(...
        updrs(:,combos(loop,1)), updrs(:,combos(loop,2)),0.05);
        somevar = somevar+1;
    
end
disp(p)
disp(h)
stat_ans =[combos(:,1)';combos(:,2)'; h; p; abs([stat(:).tstat])];

%%
contS = [1:8 9 10]
subj = [1  3:10]
figure, plot([1],[(conITC_roi1(contS,1)-conITC_roi1(contS,2))], 'ob'), ...
%     (PDITC_roi1(subj,1)-PDITC_roi1(subj,2))],'ob')
xlim([0.5 2.5])
hold on, plot([2],(PDITC_roi1(subj,1)-PDITC_roi1(subj,2)),'or')
title('1vs2')



figure, plot([1],[(conITC_roi1(contS,1)-conITC_roi1(contS,6))], 'ob'), ...
%     (PDITC_roi1(subj,1)-PDITC_roi1(subj,2))],'ob')
xlim([0.5 2.5])
hold on, plot([2],(PDITC_roi1(subj,1)-PDITC_roi1(subj,6)),'or')
title('1vs6')

%%
save conITC60_100 conITC60_100
save conITC40_100 conITC40_100

%% Plot ITC TFR averaged across SEF N20 and neighbours

% cfg           = [];
% cfg.channel   = {'MEG*2','MEG*3'};
% cfg.method    = 'distance';
% cfg.template = 'neuromag306planar_neighb.mat';
% cfg.neighbourdist = 6;
% cfg.grad      = SEFITC_cmb{1,1}.grad;
% cfg.feedback  = 'yes';
% neighbours    = ft_prepare_neighbours(cfg);

cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relative';
cfg.parameter               = 'powspctrm';

for subj = 1:11
    for cond = 1:6
        ITC_bsl{subj, cond} = ft_freqbaseline(cfg,SEFITC_cmb{subj, cond});
    end
end
%%
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0 0.2];
cfg.zlim = [0 1];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
% cfg.channel =  {neighbours(1,60).neighblabel{:}, neighbours(1,60).label};
% for cond = 1:6
% subplot(2,3,cond)
% figure,ft_multiplotTFR(cfg, SEFITC_cmb{1,1})
figure,ft_multiplotTFR(cfg, ITC_bsl)


%%
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	

for subj = 1:row
    for cond = 1:col
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                    neighbours(1,maxPosN20(subj,cond)).label};
        cfg.channel =  chan_sel;
        hFig = figure(subj);
        set(hFig, 'Position', [20 60 1824 1050]);
        subplot(2,3,cond),
        ft_singleplotTFR(cfg, SEFITC_cmb{subj,cond})
        title([SEFITC_cmb{subj,1}.filename(1:3), '-', num2str(cond)])
    end
       set(gcf, 'Color', 'white'); % white bckgr
    export_fig( gcf, ...      % figure handle
    ['Corr-DBS OFF(120) vs MED ON'],...
    '-painters', ...      % renderer
    '-jpg', ...           % file format
    '-r250' );             % resolution in dpi
end


%% AVG of neigh per subj/cond
% roi_1 = {'MEG0222+0223', 'MEG0412+0413', 'MEG0422+0423', 'MEG0442+0443', ...
%             'MEG0632+0633', 'MEG0642+0643', 'MEG0712+0713', 'MEG0222+0223'};
% roi_2 = {'MEG0432+0433', 'MEG0742+0743', ...
%     'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

cfg_itc = [];
cfg_itc.combinemethod = 'itc';

cfg                     = [];
cfg.parameter           = 'powspctrm';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.ylim              = [55 95];
cfg.colorbar = 'yes';

cfg.xlim             = [0.02 0.08];
% cfg.zlim = [0 0.25];

for subj = [1 2 3:10]
    hFig = figure(subj);
    set(hFig, 'Position', [20 60 1824 1050]);
    
    for cond = 1:6
        subplot(2,3,cond)
%         ft_topoplotTFR(cfg, ft_combineplanar_itc(cfg_itc, GA_ITC_control{1, cond}))
        ft_topoplotTFR(cfg, SEF_TFR_cmb{subj,cond})
        title(['Cond-', num2str(cond)])
    end
    suptitle(['Control-', num2str(subj)])
    
%     set(gcf, 'Color', 'w'); % white bckgr
%     export_fig( gcf, ...      % figure handle
%         ['GA-Control-ITC'],... % name of output file without extension
%         '-painters', ...      % renderer
%         '-png', ...           % file format
%         '-r250' );             % resolution in dpi
end

% figure,boxplot(SEF_ITC_summary_1)
% figure,boxplot(SEF_ITC_summary_2)

%%
cfgact                     = [];
cfgact.parameter           = 'powspctrm';
cfgact.foilim              = [10 100];
cfgact.latency             = [0.02 0.08];

cfgbsl                     = [];
cfgbsl.parameter           = 'powspctrm';
cfgbsl.foilim              = [10 100];
cfgbsl.latency             = [-0.08 -0.02];

for subj = 1:120
    disp(['###### ',num2str(subj)])
%     for cond = 1:6
%         disp(['###### ',num2str(cond)])
        
        SEF_ITC_act{subj} = ft_selectdata(cfgact, SEF{subj});
        SEF_ITC_bsl{subj} = ft_selectdata(cfgbsl, SEF{subj});
        SEF_ITC_bsl{subj}.time = SEF_ITC_act{subj}.time;
        
    end
% end

save SEF_ITC_act SEF_ITC_act
save SEF_ITC_bsl SEF_ITC_bsl

%% Avg across conditions

cfg = [];
% cfg.keepindividual = 'yes';
cfg.parameter = 'powspctrm';

for cond = 1:7
   
    GA_gamma_Contorl_ind{cond} = ft_freqgrandaverage(cfg, allGamma_Control{:,cond});
    
end

%% SEF ITC cond

cfg                     = [];
cfg.parameter = 'powspctrm';
% cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.zlim = [0.6 1.2]

for cond = 1:col
        subplot(2,3,cond)
%         hFig = figure(cond);
%         set(hFig, 'Position', [20 60 1824 1050]);
        ft_topoplotER(cfg, GA_gamma_PD{cond})
%         title(['Cond-', num2str(cond)])
%         set(gcf, 'Color', 'white'); % white bckgr
        export_fig( gcf, ...      % figure handle
        ['GA-GammaPD-1-3to5-7to10-12'],... 
        '-painters', ...      % renderer
        '-jpg', ...           % file format
        '-r250' );             % resolution in dpi

end

%% Select max from Evoked power

cfg = []
cfg.parameter           = 'powspctrm';


for subj =1:10
    for cond = 1:6
        
        cfg.channel              = TFRcmb_ev{subj,cond}.label(maxPosN20(subj,cond));
        TFRcmb_ev_max{subj,cond} = ft_selectdata(cfg, TFRcmb_ev{subj,cond});
        TFRcmb_ev_max{subj,cond}.label = {'maxPosN20'};
    end
end

%% Cond-wise grand average

cfg = [];
cfg.parameter = 'powspctrm';

for cond = 1:6
    
    grandCond{cond} = ft_freqgrandaverage(cfg, TFRcmb_ev_max{:,cond});
    
end

cfg                     = [];
cfg.parameter           = 'powspctrm';
cfg.colorbar            = 'yes';
cfg.xlim                = [0 0.2];
cfg.ylim                = [10 40];
cfg.zlim                = [0 1000];
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	

for loop = 1:6
    
    subplot(2,3,loop)
    ft_singleplotTFR(cfg, grandCond{loop})
    title(['Cond: ',num2str( loop)])
end

%%



