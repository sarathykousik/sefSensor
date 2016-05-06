
clc;    clear all;      close all;

%%
filenames = dir('*.mat');

for loop = 1: length(filenames)
    
    file_sub(loop) = {filenames(loop).name};

end

file_sub = reshape(file_sub, 6, []);
file_sub = file_sub';
[row col] = size(file_sub);

%% AVG of neigh per subj/cond; TFR plot; 

% cfg_plot                     = [];
% cfg_plot.parameter           = 'powspctrm';
% cfg_plot.colorbar            = 'yes';
% cfg_plot.layout              = 'neuromag306cmb.lay';
% cfg_plot.xlim                = [0 0.2];
% cfg_plot.shading             = 'interp';
% cfg_plot.maskstyle           = 'saturation';	
% 
cfg_sel                      = [];
cfg_sel.parameter            = 'powspctrm';
cfg_sel.avgoverchan          = 'yes';

for subj = 1:row
    for cond = 1:col
        
        % load data
%         [a b c]     =   fileparts(char(file_sub(subj, cond)));
%         disp(['   ****    ',b,'     ****    '])
%         load(char(file_sub(subj, cond)));  
%         [subj,cond]
%         
        cfg = [];
        cfg.baseline                = [-0.15 -0.10];
        cfg.baselinetype            = 'relative';
        cfg.parameter               = 'powspctrm';
        TFRwave_bsl               = ft_freqbaseline(cfg, TFRwave);
        
%         cfg_desc = [];
%         cfg_desc.jackknife   = 'yes';
%         TFRInd = ft_freqdescriptives(cfg_desc, TFRwave_bsl);
%         TFRInd_cmb{subj,cond} = ft_combineplanar([], TFRInd);
%         TFRInd_cmb{subj,cond}.filename = b;
        
        % Select neigh data
        chan_sel = {neighbours(1,maxPosN20(subj,cond)).neighblabel{:}, ...
                    neighbours(1,maxPosN20(subj,cond)).label};
%         chan_sel = {neighbours(1,maxPosN20(subj,cond)).label};
        cfg_sel.channel =  chan_sel;
%         SEF_TFR_rpt{subj,cond} = ...
%             ft_selectdata(cfg_sel, ft_combineplanar([],TFRcmb_ev{subj,cond}));
%         SEF_TFR_rpt{subj,cond} = ft_combineplanar([],TFRwave_bsl);
%         SEF_TFR_rpt{subj,cond}.label = {'meanNeigh'};
%         SEF_TFR_rpt{subj,cond}.filename   = b;
        SEF_TFR_sel{subj,cond} = ft_selectdata(cfg_sel, TFRcmb_ev{subj,cond});        
%         SEF_TFR_sel{subj,cond} = ft_selectdata(cfg_sel, ft_combineplanar([],TFRwave_desc));
%         SEF_TFR_cmb{subj,cond} = ft_combineplanar([],TFRwave_desc);
        SEF_TFR_sel{subj,cond}.label = {'meanNeigh'};
%         SEF_TFR_sel{subj,cond}.filename   = b;
        
        
        % TFR plot - mean of neighbours
%         hFig = figure(subj);
%         set(hFig, 'Position', [20 60 1824 1050]);
%         subplot(2,3,cond),
%         ft_singleplotTFR(cfg_plot, SEF_TFR_sel{subj,cond})
%         title([SEF_TFR_sel{subj,cond}.filename(1:3), '-', num2str(cond)])
        
        
    end
end

%%
%     set(gcf, 'Color', 'white'); % white bckgr
%     export_fig( gcf, ...      % figure handle
%     ['Subj-', '009', '-Ind'],... % name of output file without extension
%     '-painters', ...      % renderer
%     '-jpg', ...           % file format
%     '-r250' );             % resolution in dpi

%% Grand average TFR

cfg = [];
cfg.parameter = 'powspctrm';

cfg_ind = [];
cfg_ind.parameter = 'powspctrm';
cfg_ind.keepindividual = 'yes';

for cond = 1:6
   
%     SEF_TFR_cond_avg{cond} = ft_freqgrandaverage(cfg, SEF_TFR_sel{:,cond});
    SEF_TFR_cond_ind{cond} = ft_freqgrandaverage(cfg_ind, SEF_TFR_sel{[1 3:10],cond});
    
end

% TFR condition-wise

cfg_plot                     = [];
cfg_plot.parameter           = 'powspctrm';
cfg_plot.colorbar            = 'yes';
% cfg_plot.layout              = 'neuromag306cmb.lay';
cfg_plot.xlim                = [0 0.2];
cfg_plot.ylim                = [10 100];
cfg_plot.zlim                = [0 600];
cfg_plot.shading             = 'interp';
cfg_plot.maskstyle           = 'saturation';

% % for subj = 1:row
    for cond = 1:col
        hFig = figure(cond);
%         subplot(2,3,cond)
%         set(hFig, 'Position', [20 60 1824 1050]);
        ft_singleplotTFR(cfg_plot, SEF_TFR_cond_ind{cond});
%         title(num2str(cond))

    %       set(gcf, 'Color', 'white'); % white bckgr
    %     export_fig( gcf, ...      % figure handle
    %     ['Cond-', num2str(cond), '-Ind'],... % name of output file without extension
    %     '-painters', ...      % renderer
    %     '-jpg', ...           % file format
    %     '-r250' );             % resolution in dpi
%     end    
%     suptitle(num2str(subj))
end
