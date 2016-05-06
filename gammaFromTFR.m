clc;        clear all;      close all

%%

proc.data_folder = 'J:\MEG_Research\SEF\SEFTFRwave\freq';
cd(proc.data_folder)

filenames = dir('*.mat');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub = reshape(file_sub, 6, [])';
[row col] = size(file_sub);


for subj = 1:row
    for cond = 1:col
        
        disp(['******** ', char(file_sub(subj, cond))])
        load(char(file_sub(subj, cond)));
        

        cfg = [];
        cfg.baseline                = [-0.15 -0.10];
        cfg.baselinetype            = 'relative';
        cfg.parameter               = 'powspctrm';
        TFRwave_bsl               = ft_freqbaseline(cfg, TFRwave);
        
        cfg                             = [];
        cfg.parameter                   = 'powspctrm';
        cfg.avgoverrpt                  = 'yes';
        cfg.avgoverfreq                 = 'yes';
        cfg.latency                     = [0 0.1];
        cfg.foilim                      = [40 100];
        TFRgamma_bsl_40_100{subj, cond} = ft_selectdata(cfg, TFRwave_bsl);

        cfg.foilim                      = [60 100];
        TFRgamma_bsl_60_100{subj, cond} = ft_selectdata(cfg, TFRwave_bsl);
        
        TFRgamma_bsl_40_100_cmb{subj, cond} = ...
            ft_combineplanar([],TFRgamma_bsl_40_100{subj,cond});
        
        TFRgamma_bsl_60_100_cmb{subj, cond} = ...
            ft_combineplanar([],TFRgamma_bsl_60_100{subj,cond});
        
    end
end

% for subj = 1:row
%     for cond = 1:col
%    
%         TFRgamma_bsl_40_100_cmb{subj, cond} = ft_combineplanar([],TFRgamma_bsl_40_100{subj,cond});
%         TFRgamma_bsl_60_100_cmb{subj, cond} = ft_combineplanar([],TFRgamma_bsl_60_100{subj,cond});
%         
%     end
% end

%%
TFRgamma = TFRgamma_bsl_40_100_cmb;
blankForPlot = TFRgamma{1,1};
blankForPlot.powspctrm = zeros(204,1);

for subj = 1:row
    for cond = 1:col
    
        maxGamma(subj, cond, :)  =  max(TFRgamma{subj,cond}.powspctrm, [], 3)';

    end
end

% chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
%     'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633',...
%     'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623', 'MEG1812+1813',...
%     'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

chanSel_1 = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
    'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG1612+1613',...
    'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823'};
chan_CMC_paper = {'MEG0222+0223', 'MEG0232+0233', 'MEG0412+0413', ...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', ...
    'MEG0712+0713', 'MEG0742+0743', 'MEG1622+1623', 'MEG1632+1633', ...
    'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};

% for loop = 1:length(chan_CMC_paper) 
% %     chanPos(loop) =  find(strcmp(blankForPlot.label,chan_CMC_paper(loop)));
%     chanPos(loop) ={neighbours(1,maxPosN20).neighblabel{:}, neighbours(1,60).label};
% end


% gammaCmb = maxGammaStruct{1,1};
cfg                     = [];
cfg.parameter = 'powspctrm';
cfg.colorbar = 'yes';
cfg.layout              = 'neuromag306cmb.lay';
cfg.xlim                = [0.02 0.08];
cfg.shading                 = 'interp';
cfg.maskstyle               = 'saturation';	
cfg.highlightsize    = 8;
cfg.highlightfontsize    = 8;
cfg.highlightsymbol  = 'o';

for subj = 1:row
    for cond = 1:col
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(TFRgamma{subj,cond}.label, chan_sel);
        [maxVal maxPoschk] = max(maxGamma(subj, cond,chanPos), [], 3);
%         [maxVal] = mean(maxGamma(subj, cond,chanPos), 3);
        maxPos(subj,cond) = chanPos(maxPoschk);
        cfg.highlightchannel = {TFRgamma{subj,cond}.label{maxPos(subj,cond)}};
        cfg.highlight        =  'numbers';
%         cfg.channel = {'all', '-MEG1612+1613', '-MEG2342+2343'};
%         cfg.zlim                = [0 1.5];
%         hFig = figure(subj);
%         set(hFig, 'Position', [10 80 1824 968]);
%         subplot(2,3,cond),ft_topoplotER(cfg, maxGammaStruct{subj,cond})
%         title({maxGammaStruct{subj,cond}.label{maxPos(subj,cond)}})
    end
%    suptitle(TFRgamma{subj,1}.name(1:3))
end
 
for subj = 1:row
    for cond = 1:col
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                neighbours(1,maxPosN20(subj,cond)).label};
        chanPos = match_str(TFRgamma{subj,cond}.label, chan_sel);
        
        maxGammaChExtract(subj, cond) = mean(maxGamma(subj,cond,chanPos),3);
%         maxGammaChExtract(subj, cond) = mean(maxGamma(subj,cond,maxPosN20(subj,cond)));
    end
end
figure,boxplot(maxGammaChExtract([1 2:10],:))


