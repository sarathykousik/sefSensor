%% TFRcalc
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFVisClean';
proc.save_folder                 = 'J:\MEG_Research\SEF\betaTopos';
mkdir(proc.save_folder)

%%
cd(proc.data_folder)
filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)

    % Beta mtp: Fc = 21 Hz, BW: 18 Hz
    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'MEG';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
    cfg.foi          = 21;                          
    cfg.t_ftimwin    = 4./cfg.foi;   
    cfg.toi          = -0.500:0.010:0.500;           
    cfg.tapsmofrq    = 9;
    betaFreq        = ft_freqanalysis(cfg, visClean);

    cfg = [];
    cfg.baseline                = [-0.090 -0.01];
    cfg.baselinetype            = 'absolute';
    cfg.parameter               = 'powspctrm';
    betaFreq_bsl               = ft_freqbaseline(cfg, betaFreq);

    cfg                         = [];
    betaFreq_bsl                = ft_combineplanar(cfg, betaFreq_bsl);

    cfg                         = [];
    cfg.zlim                    = 'maxmin';
    cfg.xlim                    = [0.0 0.2]; 
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.axes                    = 'no';
    ft_multiplotER(cfg, betaFreq_bsl);
    saveas(gcf,[proc.save_folder,'\',b,'-betaER'],'tif')
    
    %
    cfg                         = [];
    cfg.zlim                    = [-10e-25 30e-25];
    cfg.xlim                    = [0.03 0.09]; 
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.axes                    = 'no';
    cfg.colorbar                = 'yes';
    cfg.shading                 = 'interp';
    figure,ft_topoplotER(cfg, betaFreq_bsl);
    saveas(gcf,[proc.save_folder,'\',b,'-betaTopoEarly'],'png')
    
    cfg                         = [];
    cfg.zlim                    = [-10e-25 30e-25];
    cfg.xlim                    = [0.1 0.18]; 
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.axes                    = 'no';
    cfg.colorbar                = 'yes';
    cfg.shading                 = 'interp';
    figure,ft_topoplotER(cfg, betaFreq_bsl);
    saveas(gcf,[proc.save_folder,'\',b,'-betaTopoLate'],'png')

    
    cfg                         = [];
    cfg.zlim                    = 'maxmin';
    cfg.channel                 = {'MEG0432+0433', 'MEG0442+0443'}
    cfg.xlim                    = [-0.09 0.2]; 
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.axes                    = 'no';
    cfg.colorbar                = 'yes';
    cfg.shading                 = 'interp';
    figure,ft_singleplotER(cfg, betaFreq_bsl);
    saveas(gcf,[proc.save_folder,'\', b,'-betaEvolution'],'tif')

    % select data
    cfg = [];
    cfg.channel                 = {'MEG0432+0433', 'MEG0442+0443'}
    cfg.avgoverchan             = 'no';
    cfg.latency                 = [0.03 0.18];
    data_beta_trials           = ft_selectdata(cfg, betaFreq_bsl);
    save ([proc.save_folder,'\',b, '-betaRptChan'], 'data_beta_trials', '-v7.3');
    %     
    cfg = [];
    cfg.avgoverrpt              = 'yes';
    data_beta                  = ft_selectdata(cfg, data_beta_trials);
    beta_pow_max(loop)         = max(max(data_beta.powspctrm));
    
    close all;
    disp('#######################################')
    
end   
%%
% Add function to save images per condition


%        Check for files
%                 disp(subjName{subjLoop})
%                 disp(fileList{condLoop})
%                 disp(num2str(condName{condLoop}))
%                 exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'-SaveData.mat'],'file')