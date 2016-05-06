clc;    clear all;  close all;
ft_defaults

%%

[temp.num,temp.txt,temp.raw] = xlsread('J:\data-tsss-sef\SEF-TSSS\subjectDataList.xls');
temp.num_cell = num2cell(temp.num);
temp.data_folder                 = 'J:\procFtResults_trials';
subjName = {temp.raw{1,6:10}};
condName = {temp.raw{4:10,1}};

%% Initial processing
cfg_t= [];
cfg_t.baselinetype      = 'absolute';
cfg_t.baseline          = [-0.09 -0.01];

cfg_f = [];
cfg_f.baselinetype      = 'absolute';
cfg_f.baseline          = [-0.09 -0.01];

cfg = [];
for subjLoop = 1:length(subjName)
    %Choose filenames
    [row, col]=find(strcmp(temp.raw, subjName{subjLoop})==1);
    disp('==========================================')
    disp(['   *******      ',subjName{subjLoop}, '    *********'])
    disp('==========================================')
    
    %     fileList = dir('*.fif');
    fileList = {temp.raw{4:10,col}}; 
%     cd([temp.data_folder,'\', subjName{subjLoop}])
    
        for condLoop  = 1:length(fileList)
            load_file = [temp.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'-SaveData','.mat'];
            disp(['   ****    ',num2str(condName{condLoop}), ' ***  ',load_file,'     ****    '])
            workData{condLoop,subjLoop}     = load(load_file);
            workData{condLoop,subjLoop}.loadFile = load_file; 
            workData{condLoop,subjLoop}.origFile = workData{condLoop,subjLoop}.saveData.file; ; 
      
            cmbData_tlck{condLoop,subjLoop}     = ft_combineplanar(cfg_t, workData{condLoop,subjLoop}.saveData.tlck_avg);
            cmbData_tlck{condLoop,subjLoop}.filename = load_file;
            cmbData_tlck{condLoop,subjLoop}.origFile = workData{condLoop,subjLoop}.saveData.file;
            
            cfg = [];
            cmbData_TFR{condLoop,subjLoop}       = ft_combineplanar(cfg,...
                ft_freqbaseline(cfg_f, workData{condLoop,subjLoop}.saveData.TFR_hann));
            cmbData_TFR{condLoop,subjLoop}.filename  = load_file;
            cmbData_TFR{condLoop,subjLoop}.origFile = workData{condLoop,subjLoop}.saveData.file;
    
            
        end
end

%% Save data
% save cmbData_TFR cmbData_TFR
% save cmbData_tlck cmbData_tlck
% save workData workData

%% Grand average timelock
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_30         = ft_timelockgrandaverage(cfg,cmbData_tlck{1,:});  
GA_60        = ft_timelockgrandaverage(cfg,cmbData_tlck{2,:});
GA_90        = ft_timelockgrandaverage(cfg,cmbData_tlck{3,:});
GA_120        = ft_timelockgrandaverage(cfg,cmbData_tlck{4,:});
GA_150        = ft_timelockgrandaverage(cfg,cmbData_tlck{5,:});
GA_180        = ft_timelockgrandaverage(cfg,cmbData_tlck{6,:});
GA_210        = ft_timelockgrandaverage(cfg,cmbData_tlck{7,:});

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GGA           = ft_timelockgrandaverage(cfg,GA_30, GA_60, GA_90,GA_120,...
                GA_150, GA_180, GA_210);


%%
toi_20 = find(GA_30.time>0.020 & GA_30.time<0.030);
toi_30 = find(GA_30.time>0.082 & GA_30.time<0.121);
ep20 = zeros(size(cmbData_tlck,2), size(cmbData_tlck,1));
ep30 = zeros(size(cmbData_tlck,2), size(cmbData_tlck,1));

MEGChanSelect = {'MEG0442+0443', 'MEG1812+1813'};

chanNos = [ find(strcmp(GA_30.label, MEGChanSelect{1})) ...
    find(strcmp(GA_30.label, MEGChanSelect{2}))];

for subjLoop = 1:size(cmbData_tlck,2)
        for condLoop  = 1:size(cmbData_tlck,1)
 
            ep30(subjLoop, condLoop) ...
                = max(mean(cmbData_tlck{condLoop, subjLoop}.avg(chanNos,toi_30),1));
            ep20(subjLoop, condLoop) ...
                = max(mean(cmbData_tlck{condLoop, subjLoop}.avg(chanNos,toi_20),1));
            
%             ep20(subjLoop, condLoop) ...
%                 = max(max(cmbData_tlck{condLoop, subjLoop}.avg(chanNos,toi_20)));
%             ep30(subjLoop, condLoop) ...
%                 = max(max(cmbData_tlck{condLoop, subjLoop}.avg(chanNos,toi_30)));
        
        end
end

figure, 
subplot 221,boxplot(ep20), title('EP 20 ms'), xlabel('Conditions'), ylabel('Amplitude')
subplot 222,boxplot(ep30), title('EP 30 ms'), xlabel('Conditions'), ylabel('Amplitude')
subplot 223,plot(mean(ep20,1)), title('EP 20 ms'), xlabel('Conditions'), ylabel('Amplitude')
subplot 224,plot(mean(ep30,1)), title('EP 30 ms'), xlabel('Conditions'), ylabel('Amplitude')

%%
% values = mean(ep30,3);
[h,p,ci,stats] = ttest2(out.avg, out.avg, 0.05,'right' );
% H0: mean = 0, alpha 0.05

%%

% cfg = [];
% for loop = 1:5
    cfg              = [];
    cfg.parameter    = 'avg';
    cmbChk           = ft_combineplanar(cfg, GA_30)
    
    GA = GA_30;
    GA.avg = GA_60.avg-GA_210.avg;
    cfg = []
    cfg.showlabels   = 'no';
    cfg.baseline     = [-0.09 -0.01];
    cfg.baselinetype = 'relative';
    cfg.layout    	 = 'neuromag306cmb.lay';
    figure, ft_multiplotER(cfg,GA_30, GA_120, GA_180);
%     title(subjName{loop})
% end

%% Statistics

%

cfg = [];
cfg.method = 'template';         
% cfg.method      = 'template'; % try 'distance' as well
cfg.template    = 'neuromag306cmb_neighb.mat';               % specify type of template
cfg.layout    	= 'neuromag306cmb.lay';
cfg.feedback    = 'yes';                             % show a neighbour plot 
neighbours      = ft_prepare_neighbours(cfg, GA_30); % define neighbouring channels

%%
cfg = [];
% cfg.channel     = MEGChanSelect;
cfg.latency             = [0.02 0.03];
cfg.avgovertime         = 'yes';
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';
cfg.alpha               = 0.05;
cfg.correctm            = 'cluster';
cfg.correcttail         = 'prob';
cfg.numrandomization    = 1000;
cfg.neighbours          = neighbours;
cfg.tail                = 1;

Nsub = 5;  Ncond = 4;
% cfg.design = zeros(2*Ncond,Nsub);
% cfg.design(1:2:Nsub+2,:) = [ones(1,Nsub); 2*ones(1,Nsub); 3*ones(1,Nsub);...
%                                 4*ones(1,Nsub)];
% cfg.design(2:2:Nsub+3,:) = repmat([1 2 3 4 5],4,1);

cfg.design(1:2:Nsub,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = [1 2 3 4]; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = [1 2 3 4 5]; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,cmbData_tlck{1,:},cmbData_tlck{7,:});   % don't forget the {:}!

% plot uncorrected "significant" channels
cfg = [];
cfg.style     = 'blank';
cfg.layout    	 = 'neuromag306mag.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_30)

% figure
% plot(stat.time,stat.prob(16,:)'*10), hold on
% plot(stat.time,-stat.stat(16,:)', 'r')
% title('significant without multiple comparison correction')

%%

cfg = [];
% cfg.channel     = 'MEG0222+0223';
cfg.latency     = [0.02 0.03];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'depsamplesFmultivariate';
% cfg.alpha       = 0.025;
cfg.correctm    = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.neighbours = neighbours;
cfg.correcttail = 'none';
cfg.tail                = 1;
cfg.numrandomization = 500;
 
Nsub = 5;  Ncond = 4;
% cfg.design = zeros(2*Ncond,Nsub);
% cfg.design(1:2:Nsub+2,:) = [ones(1,Nsub); 2*ones(1,Nsub); 3*ones(1,Nsub);...
%                                 4*ones(1,Nsub)];
% cfg.design(2:2:Nsub+3,:) = repmat([1 2 3 4 5],4,1);
% cfg.ivar                = [1 2 3 4]; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = [1 2 3 4 5]; % the 2nd row in cfg.design contains the subject number

cfg.design(1,1:4*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub) 3*ones(1,Nsub)...
                            4*ones(1,Nsub)];% 5*ones(1,Nsub) 6*ones(1,Nsub) 7*ones(1,Nsub)];
cfg.design(2,1:4*Nsub)  = [1:Nsub 1:Nsub 1:Nsub 1:Nsub];% 1:Nsub 1:Nsub 1:Nsub];
cfg.uvar = 2;
cfg.ivar = 1;

stat = ft_timelockstatistics(cfg,cmbData_tlck{1,:},cmbData_tlck{3,:},...
        cmbData_tlck{5,:},cmbData_tlck{7,:});%,cmbData_tlck{5,:},...
%         cmbData_tlck{6,:},cmbData_tlck{7,:});   % don't forget the {:}!
 
% make the plot
cfg = [];
cfg.style     = 'blank';
cfg.layout    	 = 'neuromag306cmb.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GGA)
% 
% figure
% plot(stat.time,stat.prob(16,:)'), hold on
% plot(stat.time,stat.stat(16,:)', 'r')



%% Grand average TFR
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
% cfg.parameter = 'avg';
TFR_30         = ft_freqgrandaverage(cfg,cmbData_TFR{1,:});  
TFR_60         = ft_freqgrandaverage(cfg,cmbData_TFR{2,:});
TFR_90         = ft_freqgrandaverage(cfg,cmbData_TFR{3,:});
TFR_120        = ft_freqgrandaverage(cfg,cmbData_TFR{4,:});
TFR_150        = ft_freqgrandaverage(cfg,cmbData_TFR{5,:});
TFR_180        = ft_freqgrandaverage(cfg,cmbData_TFR{6,:});
TFR_210        = ft_freqgrandaverage(cfg,cmbData_TFR{7,:});

%%

cfg = [];
% cfg.baselinetype      = 'absolute';
% cfg.channel           = [16];
% cfg.baseline          = [-0.09 -0.01];
cfg.layout            = 'neuromag306cmb.lay';
% cfg.zlim              = [-1.5e-23 1.5e-23];
cfg.maskstyle         = 'saturation';
% for loop =1:7
figure,ft_multiplotTFR(cfg,cmbData_TFR{1,1});
% end

cfg            = [];
cfg.channel     = {'MEG***2', 'MEG***3'}
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 15;
cfg.foi        = 55;
freq_cmb       = ft_freqanalysis(cfg, proc.data_epoched);

[h,p,ci,stats] = ttest2(squeeze(TFR_30.powspctrm(16,:,:)),...
    squeeze(TFR_180.powspctrm(16,:,:)), 0.05,'right' );


for x=1:24
%   for y=1:46
    [h(x,:),p(x,:), ~, ~ ] =ttest2(squeeze(TFR_30.powspctrm(16,x,:)),...
        squeeze(TFR_30.powspctrm(67,x,:)),0.5,'right');
%   end
end
imagesc(h), figure, imagesc(p)

%% Gamma
foi_gamma = find(TFR_30.freq>32 & TFR_30.freq<96);
toi_gamma = find(TFR_30.time>0.010 & TFR_30.time<0.040);

foi_beta = find(TFR_30.freq>20 & TFR_30.freq<36);
toi_beta = find(TFR_30.time>0.110 & TFR_30.time<0.200);

foi_alpha = find(TFR_30.freq>8 & TFR_30.freq<15);
toi_alpha = find(TFR_30.time>0.050 & TFR_30.time<0.150);

for subjLoop = 1:size(cmbData_TFR,2)
        for condLoop  = 1:size(cmbData_TFR,1)
 
            TFR_Gamma(subjLoop, condLoop) ...
                = mean(mean(mean(cmbData_TFR{condLoop, subjLoop}.powspctrm(chanNos,foi_gamma,toi_gamma))));
            TFR_beta(subjLoop, condLoop) ...
                = mean(mean(mean(cmbData_TFR{condLoop, subjLoop}.powspctrm(chanNos,foi_beta,toi_beta))));
            TFR_alpha(subjLoop, condLoop) ...
                = nanmean(nanmean(nanmean(cmbData_TFR{condLoop, subjLoop}.powspctrm(chanNos,foi_alpha,toi_alpha),1)));
            
        end
end

%% Frequency statistics

% t_values = zscore(cmbData_TFR{1,1}.powspctrm);
% cmbData_TFR{1,1}.powspctrm = t_values;

sample_data = cmbData_TFR{1,1};
sample_data.powspctrm = zscore(sample_data.powspctrm);

cfg                   = [];
cfg.baselinetype      = 'absolute';
cfg.baseline          = [-0.09 -0.01];
freq                  = ft_freqbaseline(cfg, sample_data);

cfg = [];
% cfg.baselinetype      = 'absolute';
% cfg.channel           = [16];
% cfg.baseline          = [-0.09 -0.01];
cfg.layout            = 'neuromag306cmb.lay';
cfg.zlim              =  [-2.5 2];
cfg.maskstyle         = 'saturation';
% for loop =1:7
figure,ft_multiplotTFR(cfg,freq);







