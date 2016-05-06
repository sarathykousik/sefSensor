clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFVisClean';
proc.save_folder                 = 'J:\MEG_Research\SEF\hilbERD-ord8';
mkdir(proc.save_folder)

%%
cd(proc.data_folder)
filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)

    cfg = []; 
    cfg.hilbert = 'abs';
    cfg.bpfilter  = 'yes';
    cfg.detrend   = 'yes';
    cfg.bpfreq   = [40 90];    
    cfg.bpfiltdir = 'twopass';
    cfg.bpfilttype = 'but';
    cfg.bpfiltord     = 8;
    cfg.bpinstabilityfix = 'reduce';
    % cfg.trials = 1:100;
    hilbData = ft_preprocessing(cfg, visClean);

    %
    cfg                     = [];
    hilbComb                = ft_combineplanar(cfg, hilbData);

    cfg = [];
    tlckHilb = ft_timelockanalysis(cfg, hilbComb);

    cfg = [];
    cfg.baseline = [-0.09 -0.01];
    cfg.baselinetype = 'absolute';
    % cfg.parameter = 'trial';
    hilbBase = ft_timelockbaseline(cfg, tlckHilb);

    save ([proc.save_folder,'\',b, '-hilbBase'], 'hilbBase', '-v7.3');
end

%% Result matrix - loading and saving
cd(proc.save_folder)
filenames_ERD      = dir('*.mat');

for loop = 1: length(filenames_ERD)
    
    file_sub(loop) = {filenames_ERD(loop).name};

end

file_sub = reshape(file_sub, 6, [])';
[row col] = size(file_sub);

hilbERD=[];
for subjLoop = 1:row
    disp(['######   ', num2str(subjLoop)])
    for condLoop = 1:col
        disp(['******** ', char(file_sub(subjLoop, condLoop))])
        data = load(char(file_sub(subjLoop, condLoop)));
        hilbERD{subjLoop, condLoop} =  data.hilbBase;
        hilbERD{subjLoop, condLoop}.name = char(file_sub(subjLoop, condLoop));
    end
end

hilbERD = hilbERD';         
%%

% toi = find(hilbERD{1}.hilbBase.time>0.02 & hilbERD{1}.hilbBase.time<0.061);

% for subj = 1:6
subj = 1;    
cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.xlim = [-0.1 0.1];
cfg.ylim = 'maxmin';
ft_multiplotER(cfg, hilbERD{1,subj},hilbERD{2,subj},hilbERD{3,subj},hilbERD{4,subj},...
                hilbERD{5,subj},hilbERD{subj,6})
title(num2str(subj));

% end

%% topoplots
for subj= 1:6

    for cond = 1:6
        cfg = [];
        cfg.parameter = 'avg';
        cfg.layout    = 'neuromag306cmb.lay';
    %     cfg.channel  = {'MEG0432+0433', 'MEG0442+0443'};
        cfg.xlim = [0.02 0.06];
        cfg.zlim = [-4.43e-13 14.7e-13];
        cfg.colorbar = 'yes';
        figure
        ft_topoplotER(cfg, hilbERD{subj,cond})
        title(['Subj: ' num2str(subj), '  Cond: ' num2str(cond)']);
        saveas(gcf,['J:\MEG_Research\SEF\hilbTopo','\',num2str(subj),'-',num2str(cond),'-gammaHilb'],'jpg')
    end
end


%% Extract only 20-60 ms - mean of all trials and max - topos
% t-test on topos
% close all
subj = 6;
for cond = 1:6

cfg               = [];
cfg.latency       = [0.02 0.06];
topoGamma         = ft_selectdata(cfg, hilbERD{cond,subj});

maxGamma_ch = mean(topoGamma.avg,2);
% 
% [h p ci stats] = ttest(maxG_Z,mean(maxG_Z),0.05, 'right')


cfg               = [];
cfg.latency       = [0.02];
topoGamma_thresh  = ft_selectdata(cfg, hilbERD{cond,subj});
topoGamma_thresh.avg = maxGamma_ch.* (maxGamma_ch>(mean(maxGamma_ch)+2*std(maxGamma_ch))); 
% figure,plot(topoGamma_thresh.avg)

cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.colorbar = 'yes';
figure(cond),
subplot 121, ft_topoplotER(cfg, topoGamma_thresh), title(['Thresh', num2str(cond)])
subplot 122, ft_topoplotER(cfg, topoGamma), title('True')
suptitle([num2str(cond), ' ** ', subjID{subj}])
end

%%
cfg               = [];
cfg.latency       = [0.02 0.06];

cfg_sing               = [];
cfg_sing.latency       = [0.02];


for subj = 1:6
    for cond = 1:6
    toi        = find(hilbERD{cond,subj}.time>0.02 ...
        & hilbERD{cond,subj}.time<0.061)   ;     
    
    topoGamma         = ft_selectdata(cfg, hilbERD{cond,subj});
    maxGamma(subj,cond,:)     = max(topoGamma.avg,[],2);
    
    maxGammaStruct{subj,cond} = ...
                       ft_selectdata(cfg_sing, hilbERD{cond,subj});
   
    maxGammaStruct{subj,cond}.avg = max(topoGamma.avg,[],2);
    maxGammaStruct{subj,cond}.name= hilbERD{cond,subj}.name;
        
    end
end

maxGammaSubjAvg = squeeze(mean(maxGamma,2));


% %% Plotting subj means
% 
% for subj = 1:6
%     
%    maxGammaStruct{1,1}.avg = maxGammaSubjAvg(subj,:)';
%    cfg = [];
%    cfg.parameter = 'avg';
%    cfg.layout    = 'neuromag306cmb.lay';
%    cfg.colorbar = 'yes';
%    cfg.shading = 'interp';
%    figure(subj),
%    ft_topoplotER(cfg, maxGammaStruct{1,1}), title([num2str(subj)])
%     
% end

% Are maxima in each subject/condition ok?
blankForPlot = maxGammaStruct{1,1};
blankForPlot.avg = zeros(204,1);


chanSel = {'MEG0222+0223', 'MEG0232+0233', 'MEG0322+0323', 'MEG0332+0333', 'MEG0342+0343',...
     'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0642+0643',...
     'MEG1622+1623', 'MEG1632+1633', 'MEG1642+1643', 'MEG1812+1813', 'MEG1822+1823', ...
     'MEG1842+1843', 'MEG1912+1913'};
chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
     'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
     'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
     'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
 
 chanSel_2                 = {'MEG0432+0433', 'MEG0442+0443'};

 for loop = 1:length(chanSel_2) 
    chanPos(loop) =  find(strcmp(blankForPlot.label,chanSel_2(loop)));
 end
 
%
cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.colorbar = 'no';
cfg.shading = 'interp';
maxVal = [];  maxPos = [];
% subj=6
for subj = 1:6
 
    maxVal = max(maxGamma(subj, :,chanPos), [], 3)
    
    for cond = 1:6
        maxPos(subj,cond) = find(maxGamma(subj, cond,:)==maxVal(cond));
    end
    cfg.highlightsize    = 8;
    cfg.highlightfontsize    = 8;
    cfg.highlightsymbol  = 'o';
    cfg.highlightchannel =  {maxGammaStruct{1,1}.label{maxPos(subj,:)}}
    cfg.highlight        =  'numbers';
%     figure(cond),ft_topoplotER(cfg, maxGammaStruct{subj,:})
%     title({maxGammaStruct{1,1}.label{maxPos(subj,:)}})

end
 
% Extract max Gamma

for subj = 1:6
    for cond = 1:6
        maxGammaChExtract(subj, cond) = max(maxGamma(subj,cond,[16]));
    end
end
%%
figure
plot(maxGammaChExtract([1 2 3 4 5 6],:)')
title('all')

%% t-test between groups

for loop = 2:6
    [~, p(loop-1), ~, ~] = ...
        ttest2(maxGammaChExtract(:,1), maxGammaChExtract(:,loop),0.05);
end

disp(p)

%%

[p, table, stats] = anova1(maxGammaChExtract([1 4 5 6],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.01,'ctype','lsd');

%%

baseline = repmat(maxGammaChExtract(:,1),1,6);
maxGammaValNorm = (maxGammaChExtract-baseline);
plot(maxGammaValNorm([1 4 5 6],:))
title('all-norm abs 1 4 5 6')

% anova_rm(maxGammaValNorm([1,4,5,6],:))

%%

cfg = [];
cfg.method  = 'superimposed';
cfg.fsample = 1000;
cfg.numtrl  = 1;
cfg.trllen  = 1;
cfg.s1.freq = 65;
cfg.s1.ampl = 1;
cfg.s1.phase = 0;
cfg.noise.ampl = 0;
data = ft_freqsimulation(cfg);
figure
plot(data.time{1}, data.trial{1}(1,:))

%%
cfg                 = [];
cfg.method          = 'mtmconvol';
cfg.output          = 'pow';
cfg.pad             = 2;
cfg.foi             = 65     ;
cfg.t_ftimwin       = 5./cfg.foi;   
cfg.toi             = -0.500:0.020:0.500;           
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 12;
freq                = ft_freqanalysis(cfg, data);
plot(freq.time,squeeze(freq.powspctrm(1,:,:))./max(squeeze(freq.powspctrm(1,:,:))), 'b-');

cfg.tapsmofrq = 15;
freq          = ft_freqanalysis(cfg, data);
hold on
plot(freq.time,squeeze(freq.powspctrm(1,:,:))./max(squeeze(freq.powspctrm(1,:,:))), 'g-');


cfg.tapsmofrq = 25;
% cfg.taper  = 2;
freq          = ft_freqanalysis(cfg, data);
hold on
plot(freq.time,squeeze(freq.powspctrm(1,:,:))./max(squeeze(freq.powspctrm(1,:,:))), 'r-');


legend({'12 Hz', '15 Hz', '25 Hz'});

%%



cfg = []; 
cfg.method = 'hilbert';
% cfg.hilbert = 'abs';
% cfg.bpfilter  = 'yes';
% cfg.detrend   = 'yes';
% cfg.bpfreq   = [40 90];    
cfg.pad             = 2;
cfg.foi             = [6:2:100];
% cfg.t_ftimwin       = 5./cfg.foi;   
cfg.toi             = -0.500:0.020:0.500;       
cfg.bpfiltdir = 'twopass';
cfg.bpfilttype = 'fir';
cfg.bpfiltord     = 500;
cfg.bpinstabilityfix = 'reduce';
% cfg.trials = 1:100;
hilbData = ft_freqanalysis(cfg, visClean);

%%