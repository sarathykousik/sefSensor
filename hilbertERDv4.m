clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFVisClean\temp';
proc.save_folder                 = 'J:\MEG_Research\SEF\hilbERD-fir-500';
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
    cfg.hilbert = 'angle';
    cfg.bpfilter  = 'yes';
%     cfg.detrend   = 'yes';
    cfg.bpfreq   = [40 90];    
    cfg.bpfiltdir = 'twopass';
    cfg.bpfilttype = 'but';
    cfg.bpfiltord     = 4;
%     cfg.bpinstabilityfix = 'reduce';
    % cfg.trials = 1:100;
    hilbData = ft_preprocessing(cfg, visClean);

    cfg = [];
    tlckHilb = ft_timelockanalysis(cfg, hilbData);

    
    cfg = [];
    cfg.baseline = [-0.09 -0.05];
    cfg.baselinetype = 'relative';
    % cfg.parameter = 'trial';
    hilbBase = ft_timelockbaseline(cfg, tlckHilb);

    %
    cfg                     = [];
    hilbComb                = ft_combineplanar(cfg, hilbBase);


    
    save ([proc.save_folder,'\',b, '-hilbBase'], 'hilbBase', '-v7.3');
end

%% Result matrix - loading and saving
% cd(proc.save_folder)
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
cfg.xlim = [0.02 0.04];
cfg.ylim = 'maxmin';
ft_topoplotER(cfg, hilbBase)%{subj,1},hilbERD{subj,2},hilbERD{subj,3},hilbERD{subj,4},...
%                 hilbERD{subj,5});%,hilbERD{subj,6})
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
%         cfg.zlim = [-4.43e-13 14.7e-13];
        cfg.colorbar = 'yes';
        figure
        ft_topoplotER(cfg, hilbERD{subj,cond})
        title(['Subj: ' num2str(subj), '  Cond: ' num2str(cond)']);
%         saveas(gcf,['J:\MEG_Research\SEF\hilbTopo','\',num2str(subj),'-',num2str(cond),'-gammaHilb'],'jpg')
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

 for loop = 1:length(chanSel_1) 
    chanPos(loop) =  find(strcmp(blankForPlot.label,chanSel_1(loop)));
 end
 
%
cfg = [];
cfg.parameter = 'avg';
cfg.layout    = 'neuromag306cmb.lay';
cfg.colorbar = 'yes';
cfg.shading = 'interp';
maxVal = [];  maxPos = [];
% subj=6
%
for subj = 1:6
     maxVal = max(maxGamma(subj, :,chanPos), [], 3)
    for cond = 1:6
        maxPos(subj,cond) = find(maxGamma(subj, cond,:)==maxVal(cond));
%     end
    cfg.highlightsize    = 8;
    cfg.highlightfontsize    = 8;
    cfg.highlightsymbol  = 'o';
    cfg.highlightchannel =  {maxGammaStruct{1,1}.label{maxPos(subj,cond)}}
    cfg.highlight        =  'numbers';
%     cfg.zlim             = []
%     figure(subj),
%     subplot(2,3,cond),ft_topoplotER(cfg, maxGammaStruct{subj,cond})
%     title({maxGammaStruct{1,1}.label{maxPos(subj,cond)}})
   end
end
 
% Extract max Gamma

for subj = 1:6
    for cond = 1:6
        maxGammaChExtract(subj, cond) = max(maxGamma(subj,cond,maxPos(subj,cond)));
    end
end
%
figure
boxplot(maxGammaChExtract([1  3 4 5 6],:))
title('Abs val')

baseline = repmat(maxGammaChExtract(:,1),1,6);
maxGammaValNorm = (maxGammaChExtract-baseline)./baseline;
figure, boxplot(maxGammaValNorm([1:6],:))
title('norm')

% anova_rm(maxGammaValNorm([1,3,4],:))
%% t-test between groups

for loop = 2:6
    [~, p(loop-1), ~, ~] = ...
        ttest2(maxGammaChExtract([1 3 ],1), maxGammaChExtract(:,loop),0.05);
end

disp(p)

%%

[p, table, stats] = anova1(maxGammaValNorm([1 3 4],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.1,'ctype','dunn-sidak');



%%
% Sim 10, 15, 45, 70, 90 Hz =>  append   => TFR with same settings as dpss
loop =1;
for freqloop = [91:3:120]
    disp(['################## ', num2str(loop)])
cfg = [];
cfg.method  = 'superimposed';
cfg.fsample = 1000;
cfg.numtrl  = 50;
cfg.trllen  = 0.2;

cfg.s1.freq  = [freqloop]; 
cfg.s1.ampl = 1;
cfg.s2.phase = 0;

% cfg.s2.freq = [60];
% cfg.s2.ampl = 1;
% cfg.s2.phase = 0;

% cfg.s3.freq = [80];
% cfg.s3.ampl = 1;
% cfg.s3.phase = 0;

data_1 = ft_freqsimulation(cfg);

cfg = [];
cfg.method  = 'superimposed';
cfg.fsample = 1000;
cfg.numtrl  = 50;
cfg.trllen  = 0.2;

cfg.s1.freq  = [6]; 
cfg.s1.ampl = 1;
cfg.s1.phase = 0;

cfg.s2.freq = [12];
cfg.s2.ampl = 1;
cfg.s2.phase = 0;

cfg.s2.freq = [25];
cfg.s2.ampl = 1;
cfg.s2.phase = 0;

data_2 = ft_freqsimulation(cfg);

% data = ft_appenddata([], data_1, data_2);

data = [];
for tr  = 1:numel(data_2.trial)
    data.trial{tr} = [data_2.trial{tr} data_2.trial{tr} ...
        data_1.trial{tr} data_2.trial{tr} data_2.trial{tr}];
    data.time{tr} = 1e-3:1e-3:1.1;
    data.fsample = 1000;
    data.label = data_1.label;
    data.cfg = data_1.cfg;
    data.dimord = 'rpt_time';
end

%
% figure(1)
% plot(data.time{1}, data.trial{1}(1,:))

   cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'mix';
    cfg.method       = 'mtmconvol';
%     cfg.pad          = 2;
%     cfg.foi          = 21;                          
    cfg.foi          = 65;                          
    cfg.t_ftimwin    = 5./cfg.foi;   
    cfg.toi          = 0:0.010:1.200;           
    cfg.tapsmofrq    = 25;
%     cfg.tapsmofrq    = 8;
    gammaFreq        = ft_freqanalysis(cfg, data);
    gammaDesc{loop} = ft_freqdescriptives([],gammaFreq);
    loop=loop+1;
end
freq_str = num2str(91:3:120);
leg_freq = strsplit(freq_str);

    cfg                         = [];
    cfg.zlim                    = 'maxmin';
    cfg.channel = 'mix'
    cfg.ylim                    = [0.0 0.2]; 
%     cfg.xlim                    = [0.12 0.5]
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.axes                    = 'no';
    cfg.graphcolor              = 'brgkycmbrgkycm';
    ft_singleplotER(cfg, gammaDesc{1:10});
    legend(leg_freq{1:10})
    
title('Simulation(Noise:6,12,25) - dpss,Fc=65Hz BW=25Hz')
%     figure,plot(gammaFreq.freq,squeeze(mean(gammaFreq.powspctrm(:,1,:))))
    
%%
cfg = []; 
% cfg.hilbert = 'abs';
cfg.bpfilter  = 'yes';
% cfg.detrend   = 'yes';
cfg.channel = 'mix';
cfg.bpfreq   = [40 90];    
cfg.padding   =  2;
cfg.padtype   = 'zero';
cfg.bpfiltdir = 'twopass';
cfg.bpfiltwintype  = 'kaiser';
cfg.bpfilttype = 'fir';
cfg.bpfiltord     = 500;
% cfg.bpinstabilityfix = 'reduce';
% cfg.trials = 1:100;
gammaData = ft_preprocessing(cfg, data);

cfg = []; 
cfg.hilbert = 'abs';
cfg.bpfilter  = 'yes';
% cfg.detrend   = 'yes';
cfg.channel = 'mix';
cfg.bpfreq   = [40 90];    
cfg.padding   =  2;
cfg.padtype   = 'zero';
cfg.bpfiltdir = 'twopass';
cfg.bpfiltwintype  = 'kaiser';
cfg.bpfilttype = 'fir';
cfg.bpfiltord     = 500;
% cfg.bpinstabilityfix = 'reduce';
% cfg.trials = 1:100;
hilbData_1 = ft_preprocessing(cfg, data);

cfg.bpfiltord     = 1000;
hilbData_2 = ft_preprocessing(cfg, data);

cfg= [];
hilbDatatlk_1 = ft_timelockanalysis(cfg,hilbData_1);
hilbDatatlk_2 = ft_timelockanalysis(cfg,hilbData_2);

% plot(hilbDatatlk.time, hilbDatatlk.avg(1,:)),hold on, 
% plot(data.time{1,1},data.trial{1,1}(1,:), 'r')

%
% clf(gcf)
% cfg = []; 
% cfg.method = 'hilbert';
% cfg.hilbert = 'abs';
% cfg.bpfilter  = 'yes';
% cfg.detrend   = 'yes';
% cfg.bpfreq   = [40 90];    
% cfg.pad             = 4;
% cfg.foi             = [40:5:90];
% cfg.t_ftimwin       = 0.01;   
% cfg.toi             = 0:0.02:1;       
% % cfg.correctt_ftimwin  = 'yes';
% cfg.bpfiltdir = 'twopass';
% cfg.bpfilttype = 'fir';
% cfg.bpfiltord     = 1000;
% cfg.bpfiltwintype  = 'kaiser';
% % cfg.bpinstabilityfix = 'reduce';
% % cfg.trials = 150:250;
% hilbDatamtm = ft_freqanalysis(cfg, data);
% 
% subplot 211,imagesc(hilbDatamtm.time, hilbDatamtm.freq, squeeze(hilbDatamtm.powspctrm(1,:,:)))
% subplot 212, plot(hilbDatamtm.time,  squeeze((hilbDatamtm.powspctrm(1,:,:)))')

% TFR

cfg              = [];
cfg.paramter     = 'trial';
% cfg.keeptrials   = 'yes';
cfg.output       = 'pow';
cfg.channel = 'mix';
cfg.method       = 'mtmconvol';
cfg.pad          = 2;
cfg.foi          = 65;                          
cfg.t_ftimwin    = 3./cfg.foi;   
cfg.toi          = 0:0.020:1;           
cfg.tapsmofrq    = 17;
gammaFreq_1        = ft_freqanalysis(cfg, data);

cfg.tapsmofrq    = 25;
gammaFreq_2        = ft_freqanalysis(cfg, data);

cfg.t_ftimwin    = 4./cfg.foi;   
gammaFreq_3        = ft_freqanalysis(cfg, data);

cfg.t_ftimwin    = 5./cfg.foi;   
gammaFreq_4        = ft_freqanalysis(cfg, data);

%
figure
clf(gcf)
subplot 311,plot(data.time{1,1},data.trial{1,1}(1,:), 'b'), title('Simulated data')

subplot 312, hold on,
% plot(gammaFreq_1.time,squeeze(gammaFreq_1.powspctrm), '-db'),
plot(gammaFreq_2.time,squeeze(gammaFreq_2.powspctrm), '-or'),
% plot(gammaFreq_3.time,squeeze(gammaFreq_3.powspctrm), '-vk'),
plot(gammaFreq_4.time,squeeze(gammaFreq_4.powspctrm), '-xm'),
legend('smo:25+3c','sm:25+5c')
title('Power evolution TFR')

subplot 313,
plot(hilbDatatlk_1.time, hilbDatatlk_1.avg(1,:), '-r'),hold on, 
plot(hilbDatatlk_2.time, hilbDatatlk_2.avg(1,:), '-k'),hold on, 
% plot(data.time{1,1},data.trial{1,1}(1,:), 'b'),
plot(gammaData.time{1,1}, gammaData.trial{1,1}, '-b')
title('Sig + hilb pow envelope')
legend('ord:500', 'ord:1000', 'signal')

%
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.pad    = 'maxperlen';
cfg.foilim = [0 200];
cfg.taper  = 'hanning';
freq       = ft_freqanalysis(cfg, data);
% figure
% semilogy(freq.freq, freq.powspctrm(1,:), 'b-');
% 
% cfg.taper     = 'dpss';
% cfg.tapsmofrq = 2;   % i.e. no real smoothing
% freq          = ft_freqanalysis(cfg, data);
% hold on
% semilogy(freq.freq, freq.powspctrm(1,:), 'r-');
% legend({'hanning', 'dpss'});


