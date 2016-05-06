
clc; clear all; close all;
ft_defaults

addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')
proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEF-TSSS';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              =  10;
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFArtData';
mkdir(proc.save_folder)

%% Files

[temp.num,temp.txt,temp.raw] = xlsread('J:\MEG_Research\SEF\SEF-TSSS\subjectDataList.xls');
temp.num_cell = num2cell(temp.num);

subjName = {temp.raw{1,5:10}};
condName = {temp.raw{4:10,1}};

%% Import file

for subjLoop = 1:6
%   Choose filenames
    [row, col]=find(strcmp(temp.raw, subjName{subjLoop})==1);
    disp('==========================================')
    disp('==========================================')
    disp(['   *******      ',subjName{subjLoop}, '    *********'])
    disp('==========================================')
    disp('==========================================')
    mkdir([proc.save_folder,'\',subjName{subjLoop}]);

    fileList = {temp.raw{4:10,col}}; 
%     cd([proc.data_folder,'\', subjName{subjLoop}])
    cd(proc.data_folder)
        for condLoop  = [1,2,4:length(fileList)]
            disp('#######################################')
            disp(['   ****    ',num2str(condName{condLoop}),'     ****    '])
            disp('#######################################')

%             disp(subjName{subjLoop})
%             disp(fileList{condLoop})
%             disp(num2str(condName{condLoop}))
%             checkFlag(subjLoop, condLoop) = ...
%                 exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'],'file');

            file_name = [proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'];
            cfg                         = [];
            cfg.dataset                 = file_name ;
            proc.data_import            = ft_preprocessing(cfg);

            % Find epochs
            cfg = [];
            cfg.trialdef.prestim        = 0.5;
            cfg.trialdef.poststim       = 0.5;
            cfg.trialdef.eventtype      = proc.stim_chan ; 
            cfg.trialdef.eventvalue     = 5;
            cfg.dataset                 = file_name;
            proc.data_trial             = ft_definetrial(cfg);

            cfg = [];
            cfg.channel                 =  {'MEG'}; 
            proc.data_epoched           = ft_redefinetrial(proc.data_trial, proc.data_import);

            cfg = [];
            cfg.removemean     = 'yes';
            tlck306{subjLoop, condLoop} = ft_timelockanalysis(cfg, proc.data_epoched);
            
        
        end
end

%%

for subj=1:6
    for cond = 1:6
        
    cfg = [];
    tlk204{subj,cond} = ft_combineplanar(cfg, tlck306cond{subj,cond});

    end
end

%%

cfg                 = [];
cfg.parameter       = 'avg';
cfg.layout          = 'neuromag306cmb.lay';
cfg.xlim            = [-.10 0.120];
cfg.ylim            = 'maxmin';
ft_multiplotER(cfg, tlk204{1,1})

%%
stimChan = {'MEG0222+0223', 'MEG0232+0233', 'MEG0332+0333', 'MEG0412+0413', ...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0622+0623', 'MEG0632+0633',...
    'MEG0642+0643', 'MEG0712+0713', 'MEG0722+0723', 'MEG0732+0733', 'MEG0742+0743',...
    'MEG1032+1033', 'MEG1042+1043', 'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', ...
    'MEG1142+1143', 'MEG1242+1243', 'MEG1312+1313', 'MEG1342+1343', 'MEG1622+1623',...
    'MEG1632+1633', 'MEG1642+1643', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833',...
    'MEG1842+1843', 'MEG1912+1913', 'MEG1942+1943', 'MEG2012+2013', 'MEG2022+2023',...
    'MEG2032+2033', 'MEG2042+2043', 'MEG2212+2213', 'MEG2222+2223', 'MEG2232+2233',...
    'MEG2242+2243', 'MEG2312+2313', 'MEG2322+2323', 'MEG2412+2413', 'MEG2432+2433',...
    'MEG2442+2443'};

cfg               = [];
cfg.latency       = [-0.005 0.015];

cfg_sing               = [];
cfg_sing.latency       = [0.02];

for subj = 1:6
    disp(['######   ', num2str(subj)])

    for cond = 1:6
        topo_2030         = ft_selectdata(cfg, tlk204{subj,cond});
        [Y I] =  max(topo_2030.avg(1:102,:),[],2);
        max_2030(subj,cond,:)    = Y';
%         max_2030_lat(subj,cond,:) = (I./topo_2030.fsample)+ topo_2030.time(1);
        
%         max_2030Struct{subj,cond} = ...
%                        ft_selectdata(cfg_sing, tlckVisClean{subj,cond});

%         max_2030Struct{subj,cond}.avg = max(topo_2030.avg,[],2);
%         max_2030Struct{subj,cond}.name= tlckVisClean{subj,cond}.name;
        
    end
end
% maxOfMax =  max(max_2030(:,:,chanPos),[],3);
% blankForPlot = tlk204{1,1};
% blankForPlot.powspctrm = zeros(204,1);
% 
for loop = 1:length(stimChan) 
    chanPos(loop) =  find(strcmp(tlk204{1,1}.label,stimChan(loop)));
end

% cfg                     = [];
% cfg.parameter = 'avg';
% % cfg.colorbar = 'yes';
% cfg.layout              = 'neuromag306cmb.lay';
% % cfg.xlim                = [0.02 0.06];
% cfg.shading                 = 'interp';
% cfg.maskstyle               = 'saturation';	
% cfg.highlightsize    = 8;
% cfg.highlightfontsize    = 8;
% cfg.highlightsymbol  = 'o';

for subj = 1:6
    [maxVal I] = max(max_2030(subj, :,chanPos), [], 3);
    maxPos(subj,: )= chanPos(I);
    
%     for cond = 1:6
%         cfg.highlightchannel =  {max_2030Struct{subj,cond}.label{maxPos(subj,cond)}};
%         cfg.highlight        =  'numbers';
%         figure(subj),
%         subplot(2,3,cond),ft_topoplotER(cfg, max_2030Struct{subj,cond})
%         title({max_2030Struct{subj,cond}.label{maxPos(subj,cond)}})
%    end
end
 
for subj = 1:6
    for cond = 1:6
        maxArtChExtract(subj, cond) = max(max_2030(subj,cond,maxPos(subj,cond)));
%         max2030LatExtract(subj, cond)= max(max_2030_lat(subj,cond,maxPos(subj,cond)));
    end
end

%%

maxArt = reshape(maxArtChExtract,1,[]);
maxGamma = reshape(maxGammaVal,[],1);

figure,scatter(maxArt,maxGamma)
lsline

%%

[p, table, stats] = anova1(maxGammaVal([1 2 3 4 5:6],:));

figure,
[c, m, h, nms] = multcompare(stats,'alpha',0.05,'ctype','dunn-sidak');

%%