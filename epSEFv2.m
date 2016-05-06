clc;  clear all; close all;

%%

proc =[];
proc.dataFolder = 'j:/MEG_Research/SEF/SEFVisClean/temp';
proc.saveFolder =  'j:/MEG_Research/SEF/SEFepresults/';
% mkdir(proc.saveFolder);
cd(proc.dataFolder)

filenames      = dir('*.mat');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub = reshape(file_sub, 6, [])';
[row col] = size(file_sub);

%%

cfgtlck = [];
cfgtlck.removemean = 'yes';
% cfgbsl.baseline = [-0. -0.01];


for subj = 1:row
    disp(['######   ', num2str(subj)])
    for cond = 1:col
            disp(['******** ', char(file_sub(subj, cond))])
            load(char(file_sub(subj, cond)));
            temp = ft_timelockanalysis(cfgtlck,visClean);
            
            tlck{subj,cond}  = temp;
            tlck{subj,cond}.name = char(file_sub(subj, cond));
            
            tlckCmb{subj,cond} = ft_combineplanar([],tlck{subj,cond});
            tlckCmb{subj,cond}.name = char(file_sub(subj, cond));
        
    end
end

%% Max 20-30 ms

cfg               = [];
cfg.latency       = [0.02 0.03];

cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 1:row
    disp(['######   ', num2str(subj)])

    for cond = 1:col
        topo_2030         = ft_selectdata(cfg, tlckVisClean{subj,cond});
        [Y I] =  max(topo_2030.avg,[],2);
        max_2030(subj,cond,:)    = Y;
        max_2030_lat(subj,cond,:) = (I./topo_2030.fsample)+ topo_2030.time(1);
        
        max_2030Struct{subj,cond} = ...
                       ft_selectdata(cfg_sing, tlckVisClean{subj,cond});

        max_2030Struct{subj,cond}.avg = max(topo_2030.avg,[],2);
        max_2030Struct{subj,cond}.name= tlckVisClean{subj,cond}.name;
        
    end
end

blankForPlot = max_2030Struct{1,1};
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
cfg.parameter = 'avg';
% cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
% cfg.xlim                = [0.02 0.06];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';

for subj = 1:7
    [maxVal I] = max(max_2030(subj, :,chanPos), [], 3);
    maxPos(subj,: )= chanPos(I);
    
    for cond = 1:6
        cfg.highlightchannel =  {max_2030Struct{subj,cond}.label{maxPos(subj,cond)}};
        cfg.highlight        =  'numbers';
        figure(subj),
        subplot(2,3,cond),ft_topoplotER(cfg, max_2030Struct{subj,cond})
        title({max_2030Struct{subj,cond}.label{maxPos(subj,cond)}})
   end
end
 
for subj = 1:7
    for cond = 1:6
        max2030ChExtract(subj, cond) = max(max_2030(subj,cond,maxPos(subj,cond)));
        max2030LatExtract(subj, cond)= max(max_2030_lat(subj,cond,maxPos(subj,cond)));
    end
end

maxPosN20 = maxPos;

%% Max 30-40 ms

cfg               = [];
cfg.latency       = [0.03 0.045];

cfg_sing               = [];
cfg_sing.latency       = [0.03];

for subj = 1:7
    disp(['######   ', num2str(subj)])

    for cond = 1:6
        topo_3040         = ft_selectdata(cfg, tlckVisClean{subj,cond});
        [Y I] =  max(topo_3040.avg,[],2);
        max_3040(subj,cond,:)    = Y;
        max_3040_lat(subj,cond,:) = (I./topo_3040.fsample)+ topo_3040.time(1);
        
        max_3040Struct{subj,cond} = ...
                       ft_selectdata(cfg_sing, tlckVisClean{subj,cond});

        max_3040Struct{subj,cond}.avg = max(topo_3040.avg,[],2);
        max_3040Struct{subj,cond}.name= tlckVisClean{subj,cond}.name;
        
    end
end

blankForPlot = max_3040Struct{1,1};
blankForPlot.avg = zeros(204,1);

chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
  chanSel_2                 = {'MEG0432+0433', 'MEG0442+0443'};

for loop = 1:length(chanSel_1) 
    chanPos(loop) =  find(strcmp(blankForPlot.label,chanSel_1(loop)));
end

cfg                     = [];
cfg.parameter = 'avg';
% cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
% cfg.xlim                = [0.02 0.06];
cfg.shading                 = 'flat';
cfg.maskstyle               = 'saturation';	
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';

for subj = 1:7
    [maxVal I] = max(max_3040(subj, :,chanPos), [], 3);
    maxPos(subj,: )= chanPos(I);
    for cond = 1:6
        cfg.highlightchannel =  {max_3040Struct{subj,cond}.label{maxPos(subj,cond)}};
        cfg.highlight        =  'numbers';
        figure(subj),
        subplot(2,3,cond),ft_topoplotER(cfg, max_2030Struct{subj,cond})
        title({max_2030Struct{subj,cond}.label{maxPos(subj,cond)}})
   end
end
 
for subj = 1:7
    for cond = 1:6
        max3040ChExtract(subj, cond) = max(max_3040(subj,cond,maxPos(subj,cond)));
        max3040LatExtract(subj, cond)= max(max_3040_lat(subj,cond,maxPos(subj,cond)));
    end
end


%% Max 75-115 ms

cfg               = [];
cfg.latency       = [0.07 0.110];

cfg_sing               = [];
cfg_sing.latency       = [0.07];

for subj = 1:7
    disp(['######   ', num2str(subj)])

    for cond = 1:6
        topo_77115         = ft_selectdata(cfg, tlckVisClean{subj,cond});
        [Y I] =  max(topo_77115.avg,[],2);
        max_77115(subj,cond,:)    = Y;
        max_77115_lat(subj,cond,:) = (I./topo_77115.fsample)+ topo_77115.time(1);
        
        max_77115Struct{subj,cond} = ...
                       ft_selectdata(cfg_sing, tlckVisClean{subj,cond});

        max_77115Struct{subj,cond}.avg = max(topo_77115.avg,[],2);
        max_77115Struct{subj,cond}.name= tlckVisClean{subj,cond}.name;
        
    end
end

blankForPlot = max_77115Struct{1,1};
blankForPlot.avg = zeros(204,1);

chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
  chanSel_2                 = {'MEG0432+0433', 'MEG0442+0443'};

for loop = 1:length(chanSel_1) 
    chanPos(loop) =  find(strcmp(blankForPlot.label,chanSel_1(loop)));
end

cfg                     = [];
cfg.parameter = 'avg';
% cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
% cfg.xlim                = [0.02 0.06];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';

for subj = 1:7
    [maxVal I] = max(max_77115(subj, :,chanPos), [], 3);
    maxPos(subj,: )= chanPos(I);
    for cond = 1:6
        cfg.highlightchannel =  {max_77115Struct{subj,cond}.label{maxPos(subj,cond)}};
        cfg.highlight        =  'numbers';
        figure(subj),
        subplot(2,3,cond),ft_topoplotER(cfg, max_77115Struct{subj,cond})
        title({max_77115Struct{subj,cond}.label{maxPos(subj,cond)}})
    end
end
 
for subj = 1:7
    for cond = 1:6
        max77115ChExtract(subj, cond) = max_77115(subj,cond,maxPos(subj,cond));
        max77115LatExtract(subj, cond)= max_77115_lat(subj,cond,maxPos(subj,cond));
    end
end

%%
save max2030ChExtract max2030ChExtract
save max2030LatExtract max2030LatExtract

save max3040ChExtract max3040ChExtract
save max3040LatExtract max3040LatExtract

save max77115ChExtract max77115ChExtract
save max77115LatExtract max77115LatExtract

%%
figure
select = [1 3:7];
subplot 221, boxplot(max2030ChExtract(select,:)), title('20-30 ms Amp'),  ylabel('Amplitude T')
subplot 222, boxplot(max3040ChExtract(select,:)), title('30-40 ms Amp')
% subplot 233, boxplot(max77115ChExtract(select,:)), title('77-115 ms Amp')

subplot 223, boxplot(max2030LatExtract(select,:)), title('20-30 ms Latency'), 
xlabel('Conditions'), ylabel('Latency in ms')
subplot 224, boxplot(max3040LatExtract(select,:)), title('30-40 ms Latency'), xlabel('Conditions')
% subplot 236, boxplot(max77115LatExtract(select,:)), title('77-115 ms latency'),xlabel('Conditions')

%% Grand average

cfg = [];
cfg.channel   = 'MEGGRAD';
cfg.latency   = 'all';
cfg.parameter = 'avg';
% cfg.keepindividual = 'yes';
for loop=1:6
    grandAvg{loop}      = ft_timelockgrandaverage(cfg,tlckVisClean{[1 3:7],loop});  
end

%%

cfg=[];
tlck = ft_timelockanalysis(cfg, visClean)

%%
cfg                     = [];
cfg.parameter = 'avg';
cfg.channel = 'MEGGRAD';
cfg.layout              = 'neuromag306all.lay';
cfg.hlim                = [0 0.2];
cfg.vlim    = [-5e-12 5e-12];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
ft_multiplotER(cfg,tlck)
%%
cfg.layout              = 'neuromag306cmb.lay';
figure,
ft_multiplotER(cfg,ft_combineplanar([],tlck))

%%
cfg                     = [];
cfg.parameter = 'avg';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
% cfg.xlim                = [-0.01 0.2];
cfg.shading                 = 'interp';
% cfg.maskstyle               = 'saturation';	
cfg.graphcolor = 'k';
cfg.renderer = 'zbuffer';
cfg.interpolation = 'v4' 
cfg.ylim = [0 4.2e-12]
% cfg.zlim = [0.5e-12 4e-12]
cfg.xlim = [-0.01 0.20];
% ft_topoplotER(cfg,grandAvg{1})
cfg.linewidth = 1.5;

cfg.channel = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0412+0413',...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', 'MEG0712+0713', ...
    'MEG0742+0743', 'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833',...
    'MEG1842+1843'};
cfg.graphcolor = 'rbgkcm'
ft_singleplotER(cfg,grandAvg{:})

%

hleg = legend('DBS ON', 'DBS OFF(0min)', 'DBS OFF(60min)', ...
    'DBS OFF(90min)', 'DBS OFF(120min)', 'MED ON')
set(hleg,'Fontsize',10);
set(hleg,'Fontweight','bold');
set(hleg,'Fontname','Georgia');
legend boxoff

xlab = xlabel('Time in s');
set(xlab,'Fontsize',12);
set(xlab,'Fontweight','bold');
set(xlab,'Fontname','Georgia');

ylab = ylabel('Amplitude in a.u.');
set(ylab,'Fontsize',12);
set(ylab,'Fontweight','bold');
set(ylab,'Fontname','Georgia');

title_h = title('Mean of selection')
set(title_h,'Fontsize',14);
set(title_h,'Fontweight','bold');
set(title_h,'Fontname','Georgia');

%%
subj = [1 3:5 6 7 ];
baseline = repmat(updrs(subj,1),1,3);
updrsNorm = (updrs(subj,[1 2 7])-baseline)./baseline;
updrs_res = reshape(updrsNorm', 1,[]);

baseline = repmat(max3040LatExtract(subj,1),1,3);
gammaNorm = (max3040LatExtract(subj,[1 2 6])-baseline)./baseline;
gamma_res = reshape(gammaNorm', 1,[]);

baseline = repmat(currents(subj,1),1,3);
currNorm = (currents(subj,[1 2 7])-baseline)./baseline;
curr_res = reshape(currNorm', 1,[]);

% baseline = repmat(max77115LatExtract(subj,1),1,3);
% max77115N = (max77115LatExtract(subj,[1 2 6])-baseline);%./baseline;
% max77115N_res = reshape(max3040LatExtract(subj, [1 2 6]), 1,[]);

figure(1),scatter(updrs_res, gamma_res), lsline;
legend('Data',['Rho: ', num2str(corr2(updrs_res, gamma_res))])
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
figure,plot(maxGammaChExtract_now([1 2 3:6],1:6)','o')
% axis([1 4 2.5 5 ])
% set(gca,'YTickLabel',...
%     {'Stim OFF 0', 'Stim OFF 30', 'Stim OFF 60','Stim ON'},...
%     'fontsize', 9, 'fontname', 'Georgia') ;
lsline

%%

[p, table, stats] = anova_rm(maxGammaChExtract([1 2:6],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.05,'ctype','lsd');


%%  Evoked TFR

% subj = 1;  cond = 1;
for subj = 1:6
    disp('####################')
    for cond=1 :6
        cfg              = [];
        cfg.paramter     = 'avg';
        cfg.output       = 'pow';
        cfg.channel      = procLabel{subj,cond};
        cfg.method       = 'wavelet';
        cfg.pad          = 2;
        cfg.taper        = 'dpss';
        cfg.foi          = 1:2:100;                          
        cfg.t_ftimwin    = 2./cfg.foi;   
        cfg.width        = 7; 

        cfg.toi          = -0.5:0.010:0.5;            
        cfg.tapsmofrq    = 4;
        TFRmtm           = ft_freqanalysis(cfg, tlck{subj,cond});

        cfg = [];
        cfg.baseline                = [-0.090 -0.01];
        cfg.baselinetype            = 'absolute';
        cfg.parameter               = 'powspctrm';
        TFRmtm_bsl       = ft_freqbaseline(cfg, TFRmtm);
        
        TFRcmb{subj,cond} = ft_combineplanar([],TFRmtm_bsl);
    end
end

cfg         = [];
cfg.xlim    = [0 0.3];
cfg.ylim    = [10 40];
% cfg.showlabels   = 'yes';	        
cfg.maskstyle  = 'saturation';
cfg.layout          = 'neuromag306cmb.lay';

for subj=1:6
    for cond = 1:6

        figure(subj), subplot(2,3,cond)
        ft_singleplotTFR(cfg, TFRcmb{subj,cond})
    end
end

%%

cfg_lateB                         = [];
cfg_lateB.maskstyle               = 'saturation';
cfg_lateB.zlim                    = 'maxmin';
cfg_lateB.xlim                    = [0.02 0.13]; 
cfg_lateB.ylim                    = [12 24];

cfg_lowG                         = [];
cfg_lowG.maskstyle               = 'saturation';
cfg_lowG.zlim                    = 'maxmin';
cfg_lowG.xlim                    = [0.02 0.1]; 
cfg_lowG.ylim                    = [30 40];

cfg_highG                         = [];
cfg_highG.maskstyle               = 'saturation';
cfg_highG.zlim                    = 'maxmin';
cfg_highG.xlim                    = [0.02 0.06]; 
cfg_highG.ylim                    = [50 90];

%% Induced
% cd(proc.dataFolder);
% maxPosres = reshape(maxPos',[],1);

for subj = 1:6
    for cond = 1:6
   
        label = char(tlckVisClean{1,1}.label(maxPos(subj,cond)));
        procLabel{subj,cond} = {label(1:7), ['MEG', label(9:12)]};
    end
end
maxPosEvoked = reshape(procLabel',[],1);

%% Calc Induced power
% filenames      = dir('*.mat');
TFRwave = {};   inducedBslWave = {};  TFRWaveProc  = {};  inducedBslWave = {};
% for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
    
    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.channel      = maxPosRes{loop};
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
    cfg.foi          = 65;                          
    cfg.t_ftimwin    = 5./cfg.foi;   
    cfg.toi          = -0.500:0.020:0.500;           
    cfg.tapsmofrq    = 25;
    TFRwave{loop}        = ft_freqanalysis(cfg, visClean);
%     cfg              = [];
%     cfg.paramter     = 'avg';
% %     cfg.keeptrials   = 'yes';
%     cfg.output       = 'pow';
%     cfg.channel      = maxPosRes{loop};
%     cfg.method       = 'wavelet';
%     cfg.pad          = 2;
%     cfg.foi          = 1:2:100;                          
%     cfg.width        = 5; 
%     cfg.toi          = -0.5:0.010:0.5;        
%     TFRfourier{loop}    = ft_freqanalysis(cfg, visClean);
    TFRwave{loop}.name = b;
    
    cfg                         = [];
    cfg.baseline                = [-0.090 -0.010];
    cfg.baselinetype            = 'relative';
    cfg.parameter               = 'powspctrm';
    TFRwave_bsl                 = ft_freqbaseline(cfg, TFRwave{loop});

    cfg                         = [];
    inducedBslWave{loop}        = ft_combineplanar(cfg, TFRwave_bsl);
    inducedBslWave{loop}.name   = b;

% end
% inducedBslWaveProc = reshape(inducedBslWave,6,6)';
% TFRWaveProc = reshape(TFRwave,6,6)';

%%

% cfg.layout          = 'neuromag306cmb.lay';
% zlims = [    -0.2e-22  3e-22;  ...
%            -.06e-22 2.5e-22; ...
%           -0.2e-22 35e-22  ;... 
%              -2e-22 8e-22; ...
%              -0.5e-22 6e-22; ...
%           -1e-22 10e-22];
for subj=1:6
    subj
    for cond = 1:6
        cond
        cfg                         = [];
    cfg.baseline                = [-0.090 -0.010];
    cfg.baselinetype            = 'relchange';
    cfg.parameter               = 'powspctrm';
    TFRwave_bsl                 = ft_freqbaseline(cfg, TFRWaveProc{subj,cond});

    cfg                         = [];
    TFRwaverelch{subj,cond}        = ft_combineplanar(cfg, TFRwave_bsl);
%     inducedBslWave{loop}.name   = b;


% inducedBslWaveLog{subj,cond}    = ft_math(cfg, TFRWaveProc{subj,cond});
    end
end

cfg = [];
cfg.xlim         = [0 0.250];
% cfg.zlim         = 'maxmin';%[-2e-22 35e-22];
% cfg.maskstyle    = 'saturation';

for subj= 1:6
%     cfg.zlim = zlims(subj,:);

    for cond =1:6


        h3=figure(subj);

        subplot(2,3,cond),
        ft_singleplotER(cfg, inducedBslWaveProc{subj,cond});

    end
end

% figure(1)
% suptitle('Low Gamma')
% 
% figure(2)
% suptitle('High Gamma')
% 
% figure(3)
% suptitle('Late Beta')
% 

for cond = 1:6
    cfg = [];
    cfg.parameter = 'powspctrm';
%     cfg.channel = {'maxN20'};

    TFRgrandAvg{cond} =  ft_freqgrandaverage(cfg,TFRwaverelch{:,cond});
end

cfg = [];
cfg.xlim         = [0 0.250];
cfg.zlim         = 'maxmin';%[-2e-22 35e-22];
cfg.maskstyle  = 'saturation';
    for cond =1:6
        h3=figure(cond);
%         subplot(2,3,cond),
        ft_singleplotTFR(cfg, TFRgrandAvg{cond});

    end

    % statistics on induced data
    cfg = [];
cfg.channel   = 'MEGGRAD';
cfg.statistic = 'depsamplesT';
cfg.ivar      = 1;
cfg.ivar      = 2;
cond1 = 1; cond2=6;
design = zeros(2,size(inducedBslWaveProc{1,cond1}.powspctrm,1) ...
    + size(inducedBslWaveProc{2,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{3,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{4,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{5,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{6,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{1,cond2}.powspctrm,1) ...
    + size(inducedBslWaveProc{2,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{3,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{4,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{5,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{6,cond2}.powspctrm,1));

design(1:size(inducedBslWaveProc{1,cond1}.powspctrm,1))  = 1;
design(1:size(inducedBslWaveProc{1,cond1}.powspctrm,1)...
    +size(inducedBslWaveProc{2,cond1}.powspctrm,1))  = 2;

for subj = 1:6
    
    design(1:size(inducedBslWaveProc{1,cond1}.powspctrm,1))  = subj;
    
end

design(1,1:size(inducedBslWaveProc{1,cond1}.powspctrm,1) ...
    + size(inducedBslWaveProc{2,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{3,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{4,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{5,cond1}.powspctrm,1)...
    + size(inducedBslWaveProc{6,cond1}.powspctrm,1)) = 1;

design(1,1:size(inducedBslWaveProc{1,cond2}.powspctrm,1) ...
    + size(inducedBslWaveProc{2,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{3,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{4,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{5,cond2}.powspctrm,1)...
    + size(inducedBslWaveProc{6,cond2}.powspctrm,1)) = 2;

% design(1,1:size(freqFIC_planar_cmb.powspctrm,1)) = 1;

cfg.design = design;
cfg.method    = 'analytic';
cfg.correctm  = 'no';
TFR_stat1     = ft_freqstatistics(cfg, inducedBslWaveProc{:,1},inducedBslWaveProc{:,6});
    
    
%% Connectivity between S1 and other areas

chan_L_S1S2 = {'MEG0112+0113', 'MEG0122+0123', 'MEG0132+0133', ...
    'MEG0142+0143', 'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', ...
    'MEG0242+0243', 'MEG0312+0313', 'MEG0322+0323', 'MEG0332+0333', ...
    'MEG0342+0343', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', ...
    'MEG0442+0443', 'MEG0512+0513', 'MEG0522+0523', 'MEG0532+0533', ...
    'MEG0542+0543', 'MEG0612+0613', 'MEG0642+0643', 'MEG1512+1513',...
    'MEG1812+1813', 'MEG1822+1823'};

chan_R_S1S2 = {'MEG0912+0913', 'MEG0922+0923', 'MEG0932+0933', 'MEG0942+0943', ...
    'MEG1022+1023', 'MEG1032+1033', 'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133',...
    'MEG1142+1143', 'MEG1212+1213', 'MEG1222+1223', 'MEG1232+1233', 'MEG1242+1243',...
    'MEG1312+1313', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1412+1413',...
    'MEG1422+1423', 'MEG1432+1433', 'MEG1442+1443', 'MEG2212+2213', 'MEG2222+2223',...
    'MEG2412+2413', 'MEG2612+2613', 'MEG2622+2623'};



% cfg           = [];
% cfg.method    = 'mtmfft';
% cfg.taper     = 'dpss';
% cfg.output    = 'fourier';
% cfg.tapsmofrq = 2;
% freq          = ft_freqanalysis(cfg, visClean);


%% ITC
TFRfourier = {};
ITC = {};
for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)

    cfg              = [];
    cfg.paramter     = 'avg';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'fourier';
    cfg.channel      = maxPosEvoked{loop};
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
    cfg.foi          = 1:2:100;                          
    cfg.t_ftimwin    = 5./cfg.foi;
    cfg.tapsmofrq    = 0.4 *cfg.foi;

    cfg.toi          = -0.5:0.010:0.5;        
    TFRfourier{loop}    = ft_freqanalysis(cfg, visClean);
    TFRfourier{loop}.name = b;
    
    tmpdat = TFRfourier{loop}.fourierspctrm;
    tmpdat = tmpdat./(abs(tmpdat)); % this will normalise each trial for its
    ITC{loop}= abs(mean(tmpdat)); % this will give the itc
    
end

ITCProc = reshape(ITC,6,6)';

for subj= 1:6
%     cfg.zlim = zlims(subj,:);

    for cond =1:6
        
        figure(subj)
        subplot(2,3,cond)
        
        imagesc(TFRfourier{1}.time, TFRfourier{1}.freq,squeeze(sum(ITCProc{subj,cond},2)));
        axis xy
        axis([0 0.1 0 100])
%         caxis([0 1.2])
        colorbar
        
    end
  
end

%%
cfg             = [];
% cfg.variance    = 'yes';
cfg.jackknife   = 'yes';
% cfg.foilim      = [30 90];
cfg.toilim      = [0.01 0.4];


for subj = 1:6
    for cond = 1:6
        inducedFreqDesc{subj, cond} = ft_freqdescriptives(cfg, inducedBslWaveProc{subj,cond});
        
    end
end

for subj =1:6
    for cond = 1:6
          figure(subj)
        subplot(2,3,cond)
                    errorbar(inducedFreqDesc{subj,cond}.time,squeeze(inducedFreqDesc{subj,cond}.powspctrm), ...
                squeeze(inducedFreqDesc{subj,cond}.powspctrmsem));
            axis([0 0.4 2 4.5])
    end
end

%%

cfg           = [];
cfg.channel   = 'MEG';
cfg.method    = 'triangulation';
cfg.grad      = gammaFreq.grad;
cfg.feedback  = 'yes';
neighbours    = ft_prepare_neighbours(cfg);

%%
for subj = 1:7
for cond = 1:6
    cond
    allGamma_copy{subj,cond}.dimord = 'chan_time';
    allGamma_copy{subj,cond}.avg = squeeze(allGamma_copy{subj,cond}.powspctrm);
    allGamma_copy{subj,cond} = rmfield(allGamma_copy{subj,cond}, {'powspctrm','powspctrmsem', 'freq' });
end
end
%% 
subj = [1 3:7];
cfg           = [];
cfg.statistic = 'depsamplesF';
design(1,:) = [ones(1,6), 2*ones(1,6) 3*ones(1,6) 4*ones(1,6) 5*ones(1,6) 6*ones(1,6)] ;...
%           1:6, 1:6, 1:6, 1:6, 1:6, 1:6];
design(2,:) = [1:6 1:6 1:6 1:6 1:6 1:6 ];
% design =
cfg.design = design;


% channel_sel = {'MEG0132+0133', 'MEG0212+0213', 'MEG0222+0223', ...
%     'MEG0232+0233', 'MEG0242+0243', 'MEG0412+0413', 'MEG0422+0423',...
%     'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', 'MEG0712+0713', ...
%     'MEG0742+0743', 'MEG1522+1523', 'MEG1612+1613', 'MEG1622+1623', ...
%     'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
% cfg.channel = chanSel_1;
lat_gl   = [0.02 0.1];
cfg.latency     = lat_gl;
% cfg.channel = 'MEGGRAD';
cfg.alpha       = 0.01;
cfg.clusteralpha = 0.01;
cfg.ivar      = 1;
cfg.uvar      = 2;
cfg.tail        = 1; % two-sided test
% cfg.method    = 'analytic';
cfg.correcttail = 'prob';
% cfg.correctm  = 'bonferoni';

cfg.method            = 'montecarlo';
cfg.correctm          = 'cluster';

cfg.numrandomization  = 500; % 1000 is recommended, but that takes longer
cfg.neighbours        = neighbours;
% subj
ERF_stat1     = ft_timelockstatistics(cfg, tlckVisClean{subj,1}, tlckVisClean{subj,2}, ...
                    tlckVisClean{subj,3}, tlckVisClean{subj,4}, tlckVisClean{subj,5},...
                    tlckVisClean{subj,6});

figure(1),imagesc(ERF_stat1.mask)

cfg = [];
cfg.latency     = lat_gl;
% cfg.channel = chanSel_1;
avg = ft_selectdata(cfg, tlckVisClean{1,1});

cfg = [];
cfg.layout        = 'neuromag306cmb.lay';
cfg.maskparameter = 'mask';
% cfg.xlim     = [0.02 0.1];
% cfg.channel = 'MEGGRAD';
cfg.maskstyle     = 'box';
avg.mask = ERF_stat1.mask;  % copy the significance mask into the ERF
figure(2); ft_multiplotER(cfg, avg);

%%










