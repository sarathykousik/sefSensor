clc;  clear all; close all;

%%

proc =[];
proc.dataFolder = 'j:/MEG_Research/SEF/SEFVisClean/';
proc.saveFolder =  'j:/MEG_Research/SEF/SEFep/';
mkdir(proc.saveFolder);
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
cfgbsl.baseline = [-0.09 -0.01];


for subj = 1:6
    disp(['######   ', num2str(subj)])
    for cond = 1:6
            disp(['******** ', char(file_sub(subj, cond))])
            load(char(file_sub(subj, cond)));
            temp = ft_timelockanalysis(cfgtlck,visClean);
            tlck{subj,cond}  = ft_timelockbaseline(cfgbsl, temp);
            tlckVisClean{subj,cond} = ft_combineplanar([],tlck{subj,cond});
            tlckVisClean{subj,cond}.name = char(file_sub(subj, cond));
        
    end
end

%%
cfg                 = [];
cfg.parameter       = 'avg';
cfg.layout          = 'neuromag306planar.lay';
% cfg.layout          = 'neuromag306planar.lay';
% cfg.colorbar = 'yes';
cfg.xlim            = [-.10 0.120];
cfg.ylim            = 'maxmin';
cfg.interplimits    = 'electrodes';
% cfg.baselinetype    = 'relative';
% cfg.baseline        = [-0.09 -0.01]; 
cfg.shading         = 'interp';
% cfg.channel = {'MEG0122+0123', 'MEG0132+0133', 'MEG0212+0213', ...
% cfg.channel = setdiff(tlckVisClean{1,1}.label,{'MEG0122+0123', 'MEG0132+0133', 'MEG0212+0213', ...
% 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', 'MEG0322+0323', ...
%     'MEG0332+0333', 'MEG0342+0343', 'MEG0412+0413', 'MEG0422+0423',...
%     'MEG0432+0433', 'MEG0442+0443', 'MEG1522+1523', 'MEG1612+1613', ...
%     'MEG1622+1623', 'MEG1632+1633', 'MEG1642+1643', 'MEG1812+1813', ...
%     'MEG1822+1823', 'MEG1842+1843', 'MEG0112+0113', 'MEG0142+0143', ...
%     'MEG1512+1513', 'MEG1542+1543', 'MEG0312+0313', 'MEG0512+0513', ...
%     'MEG0522+0523', 'MEG0532+0533', 'MEG0542+0543', 'MEG0612+0613', ...
%     'MEG0632+0633', 'MEG0642+0643', 'MEG0712+0713', 'MEG0742+0743', ...
%     'MEG1532+1533', 'MEG1712+1713', 'MEG1722+1723', 'MEG1732+1733', ...
%     'MEG1742+1743', 'MEG1832+1833', 'MEG1912+1913', 'MEG1922+1923', ...
%     'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2042+2043', ...
%     'MEG2142+2143'});


% for subj = 1:6
subj = 4;        
% cond = 1;
%     disp(['######   ', num2str(subj)])
%             figure(subj)
% clf(gcf)
%     for cond = 1:6
%         disp(['*** ', num2str(cond)])

%         subplot(2,3,cond)
ft_multiplotER(cfg, proc.tlck)
%         ft_multiplotER(cfg, tlckVisClean{subj,1}, tlckVisClean{subj,2},tlckVisClean{subj,3},...
%         tlckVisClean{subj,4}, tlckVisClean{subj,5},tlckVisClean{subj,6});
    
        
%     end 
% end

%%

cfg = [];
cfg.channel     = 'MEGGRAD';
cfg.latency     = [0.02 0.03];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesT'
cfg.alpha       = 0.05;
cfg.correctm    = 'yes';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.correcttail = 0;
cfg.numrandomization = 1000;
 
Nsub = 4;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% subj = 6;
% design = zeros(2,2*subj);
% for i = 1:subj
%   design(1,i) = i;
% end
% for i = 1:subj
%   design(1,subj+i) = i;
% end
% design(2,1:subj)        = 1;
% design(2,subj+1:2*subj) = 2;

cfg.ivar                = 1; % cfg.design contains the independent variable
cfg.uvar                = 2; % cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,tlckVisClean{1,[1 4:6]},tlckVisClean{4,[1 4:6]});
figure(1),plot(stat.prob);figure(gcf);

cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'neuromag306cmb.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure(2); ft_topoplotER(cfg, tlckVisClean{1,1})
title('Nonparametric: significant without multiple comparison correction')

% figure(3), imagesc(stat.time, [1:102],stat.prob)

%% Max 20-30 ms

cfg               = [];
cfg.latency       = [0.02 0.03];

cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 1:6
    disp(['######   ', num2str(subj)])

    for cond = 1:6
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

for subj = 1:6
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
 
for subj = 1:6
    for cond = 1:6
        max2030ChExtract(subj, cond) = max(max_2030(subj,cond,maxPos(subj,cond)));
        max2030LatExtract(subj, cond)= max(max_2030_lat(subj,cond,maxPos(subj,cond)));
    end
end

maxPosN2 = maxPos;
%% Max 30-40 ms

cfg               = [];
cfg.latency       = [0.03 0.045];

cfg_sing               = [];
cfg_sing.latency       = [0.03];

for subj = 1:6
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

for subj = 1:6
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
 
for subj = 1:6
    for cond = 1:6
        max3040ChExtract(subj, cond) = max(max_3040(subj,cond,maxPos(subj,cond)));
        max3040LatExtract(subj, cond)= max(max_3040_lat(subj,cond,maxPos(subj,cond)));
    end
end


%% Max 75-115 ms

cfg               = [];
cfg.latency       = [0.075 0.115];

cfg_sing               = [];
cfg_sing.latency       = [0.09];

for subj = 1:6
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

for subj = 1:6
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
 
for subj = 1:6
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
select = [1   2 3:6];
subplot 231, boxplot(max2030ChExtract(select,:)), title('20-30 ms Amp'),  ylabel('Amplitude T')
subplot 232, boxplot(max3040ChExtract(select,:)), title('30-40 ms Amp')
subplot 233, boxplot(max77115ChExtract(select,:)), title('77-115 ms Amp')

subplot 234, boxplot(max2030LatExtract(select,:)), title('20-30 ms Latency'), 
xlabel('Conditions'), ylabel('Latency in ms')
subplot 235, boxplot(max3040LatExtract(select,:)), title('30-40 ms Latency'), xlabel('Conditions')
subplot 236, boxplot(max77115LatExtract(select,:)), title('77-115 ms latency'),xlabel('Conditions')


%%

[p, table, stats] = anova1(max77115LatExtract([1 3 5:6],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.05,'ctype','lsd');

%%  Evoked TFR

% subj = 1;  cond = 1;
for subj = 1:6
    for cond=1 :6
        cfg              = [];
        cfg.paramter     = 'avg';
        cfg.output       = 'pow';
        cfg.channel      = 'MEGGRAD';
        cfg.method       = 'wavelet';
        cfg.pad          = 2;
        % cfg.taper        = 'dpss';
        cfg.foi          = 1:2:100;                          
        % cfg.t_ftimwin    = 2./cfg.foi;   
        cfg.width        = 7; 

        cfg.toi          = -0.5:0.010:0.5;            
        % cfg.tapsmofrq    = 4;
        TFRmtm           = ft_freqanalysis(cfg, tlck{subj,cond});

        cfg = [];
        cfg.baseline                = [-0.090 -0.01];
        cfg.baselinetype            = 'absolute';
        cfg.parameter               = 'powspctrm';
        TFRmtm_bsl       = ft_freqbaseline(cfg, TFRmtm);
        
        TFRcmb{subj,cond} = ft_combineplanar([],TFRmtm_bsl);
    end
end

for loop=1:6
cfg = [];
% cfg.zlim         = [-3e-25 3e-25];
% cfg.showlabels   = 'yes';	        
  cfg.maskstyle  = 'saturation';
cfg.layout          = 'neuromag306cmb.lay';
figure(loop)
ft_multiplotTFR(cfg, TFRcmb{4,loop})

end


cfg_lateB                         = [];
cfg_lateB.maskstyle               = 'saturation';
cfg_lateB.zlim                    = 'maxmin';
cfg_lateB.xlim                    = [0.2 0.3]; 
cfg_lateB.ylim                    = [22 30];

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

%%

cfg = [];
cfg.xlim         = [0 0.250];
% cfg.zlim         = [-2e-22 35e-22];
cfg.maskstyle  = 'saturation';
cfg.layout          = 'neuromag306cmb.lay';
zlims = [    -0.2e-22  3e-22;  ...
           -.06e-22 2.5e-22; ...
          -0.2e-22 35e-22  ;... 
             -2e-22 8e-22; ...
             -0.5e-22 6e-22; ...
          -1e-22 10e-22];

for subj= 1:6
    cfg.zlim = zlims(subj,:);
    for cond =1:6

    %     h1=figure(1)
    %     subplot(2,3,cond),
    %     ft_topoplotTFR(cfg_lowG, TFRcmb{subj,cond});
    %     
    %     h2=figure(2)
    %     subplot(2,3,cond),
    %     ft_topoplotTFR(cfg_highG, TFRcmb{subj,cond});
        cfg.channel = TFRcmb{subj,cond}.label{maxPos(subj,cond)};
        h3=figure(subj);
        subplot(2,3,cond),
        ft_singleplotTFR(cfg, TFRcmb{subj,cond});

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










