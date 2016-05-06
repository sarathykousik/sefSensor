clc;  clear all; close all;

%%
proc =[];
proc.dataFolder = 'J:\MEG_Research\SEF\controlNew\xscan30\others';
% proc.saveFolder =  'J:\MEG_Research\SEF\SEFGamma\2tap-75Hz-10ms\freq\temp\';
% mkdir(proc.saveFolder);
cd(proc.dataFolder)

filenames      = dir('*SaveData.mat');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub = reshape(file_sub, 3, [])';
[row col] = size(file_sub);

cfgtlck = [];
cfgtlck.channel = 'MEGGRAD';
cfgtlck.removemean = 'yes';

% cfgbsl = [];
% cfgbsl.removemean = 'yes';
% cfgbsl.baseline = [-0.05 -0.01];

for subj = 1:row
    disp(['######   ', num2str(subj)])
    for cond = 1:col
            disp(['******** ', char(file_sub(subj, cond))])
            load(char(file_sub(subj, cond)));
                        
            tlck{subj,cond}           = ft_timelockanalysis(cfgtlck,saveData.dataArtRej);
            tlck{subj,cond}.name      = char(file_sub(subj, cond));
            
            tlckCmb{subj,cond}        = ft_combineplanar([],tlck{subj,cond});
            tlckCmb{subj,cond}.name   = char(file_sub(subj, cond));
        
    end
end

% save tlckControl tlckControl
% save tlckCmbControl tlckCmbControl

%% Plot gamma
figure
cfg                     = [];
% cfg.parameter           = 'powspctrm';
% cfg.colorbar            = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.shading             = 'interp';
cfg.maskstyle           = 'saturation';	
cfg.hlim = [0 0.2];
cfg.channel = {'all','-MEG1442+1443'};
% for cond=1:10
%     subplot(3,4,cond)
%     ft_topoplotER(cfg, tlckCmb_new{cond})
% end
tlckCmb{1,1}=ft_combineplanar([],tlck)
    ft_multiplotER(cfg, ft_combineplanar([],tlck))

%% Max 20-30 ms

cfg               = [];
cfg.latency       = [0.02 0.03];

cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 1:row
    disp(['######   ', num2str(subj)])

    for cond = 1:col
        topo_2030         = ft_selectdata(cfg, tlckCmb{subj,cond});
        [Y I] =  max(topo_2030.avg,[],2);
        max_2030(subj,cond,:)    = Y;
        max_2030_lat(subj,cond,:) = (I./topo_2030.fsample)+ topo_2030.time(1);
        
        max_2030Struct{subj,cond} = ...
                       ft_selectdata(cfg_sing, tlckCmb{subj,cond});

        max_2030Struct{subj,cond}.avg = max(topo_2030.avg,[],2);
        max_2030Struct{subj,cond}.name= tlckCmb{subj,cond}.name;
        
    end
end

blankForPlot = max_2030Struct{1,1};
blankForPlot.powspctrm = zeros(204,1);

chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
%   chanSel_2                 = {'MEG0432+0433', 'MEG0442+0443'};

for loop = 1:length(chanSel_1) 
    chanPos(loop) =  find(strcmp(blankForPlot.label,chanSel_1(loop)));
end

%
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
% cfg.channel = {'all','-MEG2422+2423','-MEG1342+1343'};

for subj = 1:row
    [maxVal I] = max(max_2030(subj, :,chanPos), [], 3);
    maxPos(subj,: )= chanPos(I);
    
    for cond = 1:col
        cfg.highlightchannel =  {max_2030Struct{subj,cond}.label{maxPos(subj,cond)}};
        cfg.highlight        =  'numbers';
%         cfg.channel = {'all','-MEG2422+2433','-MEG1342+1343'};

        hFig = figure(subj);
        set(hFig, 'Position', [10 80 1824 968]);
        subplot(3,3,cond),ft_topoplotER(cfg, max_2030Struct{subj,cond})
        title({max_2030Struct{subj,cond}.label{maxPos(subj,cond)}})
    end
   suptitle(tlckCmb{1,1}.name(1:7))
%   set(gcf, 'Color', 'white'); % white bckgr
%     export_fig( gcf, ...      % figure handle
%     ['N20m-', num2str(subj)],... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpeg', ...           % file format
%     '-r250' );             % resolution in dpi
  
end
 
for subj = 1:row
    for cond = 1:col
        max2030ChExtract(subj, cond) = max(max_2030(subj,cond,maxPos(subj,cond)))
        max2030LatExtract(subj, cond)= max(max_2030_lat(subj,cond,maxPos(subj,cond)));
    end
end

% maxPosN20PD = maxPos;
% save maxPosN20PD maxPosN20PD;
% 
% maxPosN20Amp = max2030ChExtract;
% save maxPosN20Amp maxPosN20Amp;
% 
% maxPosN20Lat = max2030LatExtract;
% save maxPosN20Lat maxPosN20Lat;


%%

hold on
plot(mean(N20mAmpControl), '-ob')
plot(mean(N20mAmpPD), '-vr')

plot(mean(N20mLatControl), '--ob')
plot(mean(N20mLatPD), '--vr')

xlim([0.5 6.5])
legend('Control - Amp', 'PD - Amp', 'Control - lat', 'PD lat')
title('N20m maxima amplitudes')

%% Grand average
% allTlck = [tlckCmbPD; tlckCmbControl];

cfg = [];
% cfg.channel   = 'MEGGRAD';
cfg.latency   = 'all';
cfg.parameter = 'avg';
for subj = 1:20
%     for cond = 1:7
%         cfg.channel   = tlckCmb_bsl{subj,cond}.label(maxPosN20(subj, cond));
%         maxChanSelect{subj,cond} = ft_selectdata(cfg, tlckCmb_bsl{subj,cond});
%         maxChanSelect{subj,cond}.label = {'max'};
%         maxChanSelect{subj,cond}.dimord = 'chan_time';
%         grandAvgPD{cond}      = ft_timelockgrandaverage(cfg,tlckCmbPD{:,cond});  
%         grandAvgControl{cond}      = ft_timelockgrandaverage(cfg,tlckCmbControl{:,cond});  
%         grandAvgCondsControl{cond}      = ft_timelockgrandaverage(cfg,allTlck{11:20,cond});  
        grandAvgSubj{subj}      = ft_timelockgrandaverage(cfg,allTlck{subj,:});  

end
% end
% GGA      = ft_timelockgrandaverage(cfg,grandAvgall{:});  


%%
cfg                     = [];
cfg.parameter = 'avg';
cfg.colorbar = 'no';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0 0.12];
cfg.ylim                = [0 5e-12];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.graphcolor = 'brgkmcbrgkmcbrgkmc';
% cfg.colorbar = 'yes';
% cfg.channel = {'MEG0232+0233', 'MEG0432+0433', 'MEG0442+0443', ...
%                'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823'}

figure,
loop=1
for subjLoop = 6:10
%     for condLoop = 1:7

        subplot(2,3,loop)
        ft_singleplotER(cfg, grandAvgSubj{subjLoop})
        title(['PD-',num2str(subjLoop)])
        
%     end
loop=loop+1;
end

% subplot(2,1,1)
% ft_singleplotER(cfg, grandAvgSubj{[1:10]})
% subplot(2,1,2)
% ft_singleplotER(cfg, grandAvgSubj{[11:20]})


suptitle('brgkmcb')

%%
title_str='ERF-C-6-10';
% title(title_str)
export_fig( gcf, title_str,'-transparent', ...
        '-painters','-pdf', '-r250' ); 


%%
cfg = [];
% cfg.channel   = 'MEGGRAD';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';

for loop=1:6
    loop
    grandAvg_ind{loop}      = ft_timelockgrandaverage(cfg, maxChanSelect{[1:10],loop});  %tlckCmb{1,:});
end

grandAvgPDind      = ft_timelockgrandaverage(cfg,grandAvg_ind{:});  

grandAvgPDind.SEM(1,:,:) = ((mean(grandAvgPDind.individual,1))...
         -(std(grandAvgPDind.individual,1)/size(grandAvgPDind.individual,1)));
     
grandAvgPDind.SEM(2,:,:) = ((mean(grandAvgPDind.individual,1))...
         +(std(grandAvgPDind.individual,1)/size(grandAvgPDind.individual,1)));
grandAvgPDind.avg = squeeze(mean(grandAvgPDind.individual,1));


% grandAvgCind      = grandAvg_ind{:,:};  

% for cond = 1:6
%     cond
%     grandAvgCind{cond}   = grandAvg_ind{cond};
%     grandAvgCind{cond}.SEM(1,:) = ((mean(grandAvgCind{cond}.individual,1))...
%              -(std(grandAvgCind{cond}.individual,1)/size(grandAvgCind{cond}.individual,1)));
% 
%     grandAvgCind{cond}.SEM(2,:) = ((mean(grandAvgCind{cond}.individual,1))...
%              +(std(grandAvgCind{cond}.individual,1)/size(grandAvgCind{cond}.individual,1)));
%     grandAvgCind{cond}.avg = squeeze(mean(grandAvgCind{cond}.individual,1))';
% end
for cond = 1:6
    cond
    grandAvgPDind{cond}   = grandAvg_ind{cond};
    grandAvgPDind{cond}.SEM(1,:) = ((mean(grandAvgPDind{cond}.individual,1))...
             -(std(grandAvgPDind{cond}.individual,1)/size(grandAvgPDind{cond}.individual,1)));

    grandAvgPDind{cond}.SEM(2,:) = ((mean(grandAvgPDind{cond}.individual,1))...
             +(std(grandAvgPDind{cond}.individual,1)/size(grandAvgPDind{cond}.individual,1)));
    grandAvgPDind{cond}.avg = squeeze(mean(grandAvgPDind{cond}.individual,1))';
end


%%
% for cond =1:6
cond=6    
h1 = boundedline(grandAvgPDind{cond}.time, grandAvgPDind{cond}.avg, ...
    [grandAvgPDind{cond}.SEM(2,  :)'-...
     grandAvgPDind{cond}.SEM(1,  :)'], '-b', 'transparency', 0.09)
%%  
% end
% boundedline(grandAvgPDind.time, squeeze(mean(grandAvgPDind.avg)), ...
%     [squeeze(mean(grandAvgPDind.SEM(2, :),2))'-...
%       squeeze(mean(grandAvgPDind.SEM(1, :),2))'], '-r')
  axis([0 0.2 0 1.2e-11])
%   legend('Control-SEM','Control','PD-SEM','PD' )
  xlabel('Time in s'), ylabel('Amplitude in a.u.   (1e-12 a.u. = 10 fT/cm)')
  title('ERF Grand averages of PD cohort over all conditions')
% 
%   plot(grandAvgCind{cond}.time, grandAvgCind{cond}.avg, 'k'), hold on
% plot(grandAvgCind{cond}.time, grandAvgCind{cond}.SEM(1, :), 'b')
% plot(grandAvgCind{cond}.time, grandAvgCind{cond}.SEM(2, :), 'r')


%%
cfg = [];
cfg.channel = {'MEG0232+0233', 'MEG0432+0433', 'MEG0442+0443', 'MEG1622+1623', 'MEG1812+1813'};
cfg.avgoverchan = 'yes';
avgControl = ft_selectdata(cfg, grandAvgCind)

chan_nos = match_str(grandAvgPDind.label,channel);
figure
plot(grandAvgPDind.time, squeeze(mean(grandAvgPDind.avg(chan_nos,:))), 'k'), hold on
plot(grandAvgPDind.time, squeeze(mean(grandAvgPDind.SEM(1, chan_nos, :),2)), 'b')
plot(grandAvgPDind.time, squeeze(mean(grandAvgPDind.SEM(2, chan_nos, :),2)), 'r')

boundedline(grandAvgCind.time, squeeze(mean(grandAvgCind.avg(chan_nos,:))), ...
    [squeeze(mean(grandAvgCind.SEM(2, chan_nos, :),2))'-...
      squeeze(mean(grandAvgCind.SEM(1, chan_nos, :),2))'], '-b')
boundedline(grandAvgPDind.time, squeeze(mean(grandAvgPDind.avg(chan_nos,:))), ...
    [squeeze(mean(grandAvgPDind.SEM(2, chan_nos, :),2))'-...
      squeeze(mean(grandAvgPDind.SEM(1, chan_nos, :),2))'], '-r')
  axis([0 0.2 0 6e-12])
  legend('Control-SEM','Control','PD-SEM','PD' )
  xlabel('Time in s'), ylabel('Amplitude in a.u.   (1e-12 a.u. = 10 fT/cm)')
  title('ERF Grand averages of PD and control cohorts')
  
  %%
% for subj = 1:10  
% figure(subj);
    set(gcf, 'Color', 'w'); % white bckgr
    export_fig( gcf, ...      % figure handle
    'ITCall',... % name of output file without extension
    '-painters', ...      % renderer
    '-png', ...           % file format
    '-r250' );             % resolution in dpi
% end



%%  Evoked TFR

% subj = 1;  cond = 1;
for subj = 1:10
    disp(['#################### ', num2str(subj)])
    for cond= 1:6
        disp(['#################### ', num2str(cond)])
        cfg              = [];
        cfg.paramter     = 'avg';
        cfg.output       = 'pow';
        cfg.channel      = 'MEGGRAD';%procLabel{subj,cond};
        cfg.method       = 'wavelet';
        cfg.pad          = 2;
%         cfg.taper        = 'dpss';
        cfg.foi          = 6:2:100;                          
%         cfg.t_ftimwin    = 2./cfg.foi;   
        cfg.width        = 5; 
        cfg.toi          = -0.2:0.020:0.3;            
%         cfg.tapsmofrq    = 4;
        TFRwave_ev{subj, cond}       = ft_freqanalysis(cfg, tlck{subj,cond});

        cfg = [];
        cfg.baseline                = [-0.15 -0.10];
        cfg.baselinetype            = 'relative';
        cfg.parameter               = 'powspctrm';
        TFRwave_ev_bsl{subj, cond}       = ft_freqbaseline(cfg, TFRwave_ev{subj, cond});
        
        TFRcmb_ev{subj, cond} = ft_combineplanar([],TFRwave_ev_bsl{subj, cond});
    end
end




%% 
cfg                     = [];
cfg.statistic           = 'ft_statfun_depsamplesFunivariate';
% cfg.statistic           = 'ft_statfun_depsamplesT';
lat_gl                  = [0 0.2];
cfg.latency             = lat_gl;
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.ivar                = 1;
cfg.uvar                = 2;
cfg.clustertail         = 1;
cfg.tail                = 1; 
cfg.correcttail         = 'prob';
cfg.clusterstatistic    = 'maxsum';
cfg.method              = 'montecarlo';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 500; 
cfg_neighb.method       = 'distance';
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb, tlckCmb{1,1}.grad);
    
subj                    = [1 2:6];
nsubj                   = length(subj);
design                  = [];
design(1,:) = [ones(1,nsubj), 2*ones(1,nsubj) 3*ones(1,nsubj)...
    4*ones(1,nsubj) 5*ones(1,nsubj) 6*ones(1,nsubj)];
design(2,:) = [1:nsubj 1:nsubj 1:nsubj 1:nsubj 1:nsubj 1:nsubj ];
cfg.design              = design;
cfg.ivar                = 1;
cfg.uvar                = 2;
ERF_stat1               = ft_timelockstatistics(cfg,  tlckCmb{subj,:});%, tlckCmb{subj,6});




