%% TFRcalc
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFVisClean\temp';
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFGammanewtSSS\';
mkdir(proc.save_folder)

cd(proc.data_folder)
filenames      = dir('*visClean.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
%
    % Gamma mtp: Fc = 75 Hz, BW: 20 Hz - uses 2 tapers
    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'MEGGRAD';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
%     cfg.foi          = 21;                          
    cfg.foi          = 75;                          
    cfg.t_ftimwin    = 5./cfg.foi;   
    cfg.toi          = -0.500:0.010:0.500;           
    cfg.tapsmofrq    = 20;
%     cfg.tapsmofrq    = 8;
    gammaFreq        = ft_freqanalysis(cfg, visClean);
    save([proc.save_folder , b, '-gammaFreq'], 'gammaFreq', '-v7.3');
        disp('#######################################')
    
end

%% ITPC
    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'fourier';
    cfg.channel      = 'MEGGRAD';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
%     cfg.foi          = 21;                          
    cfg.foi          = 75;                          
    cfg.t_ftimwin    = 5./cfg.foi;   
    cfg.toi          = -0.500:0.010:0.500;           
    cfg.tapsmofrq    = 20;
%     cfg.tapsmofrq    = 8;


for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
%
    % Gamma mtp: Fc = 75 Hz, BW: 20 Hz - uses 2 tapers

    TFRfourier        = ft_freqanalysis(cfg, visClean);
    
    tmpdat = TFRfourier.fourierspctrm;
    tmpdat = tmpdat./(abs(tmpdat)); % this will normalise each trial for its
    ITC = abs(mean((tmpdat))); % this will give the itc
        
    TFRfourier = rmfield(TFRfourier, {'fourierspctrm','trialinfo', 'cumtapcnt'});
    TFRfourier.dimord = 'chan_freq_time';

    TFRfourier.powspctrm = permute(squeeze(ITC), [1 3 2])
    save (['J:\MEG_Research\SEF\SEF_ITC_Control_mtm\',b, '-ITC'], 'TFRfourier', '-v7.3');

end

%%
        
cfg = [];
cfg.baseline                = [-0.15 -0.1];
cfg.baselinetype            = 'relchange';
cfg.parameter               = 'powspctrm';
gammaFreq_bsl               = ft_freqbaseline(cfg,TFRwave);

cfg             = [];
cfg.jackknife   = 'yes';
gammaFreq_bsl               = ft_freqdescriptives(cfg,gammaFreq_bsl);

cfg                         = [];
gammaCmb               = ft_combineplanar(cfg, gammaFreq_bsl);
        
%%
figure
cfg                         = [];
% cfg.zlim                    = [0 7];%'maxmin';
cfg.xlim                    = [-0.1 0.2]; 
cfg.ylim                    = [20 100];
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.shading                 = 'interp';
cfg.colorbar                = 'yes';
cfg.axes                    = 'no';
% cfg.channel                 = {'MEG'};
ft_multiplotTFR(cfg, gammaCmb);
    
%%    
%     saveas(gcf,[proc.save_folder,'\',b,'-gammaER'],'tif')
    
    %
%     cfg                         = [];
%     cfg.zlim                    = [-10e-25 30e-25];
%     cfg.xlim                    = [0.02 0.06]; 
%     cfg.maskstyle               = 'saturation';	
%     cfg.masknans                = 'yes';
%     cfg.layout                  = 'neuromag306cmb.lay';
%     cfg.axes                    = 'no';
%     cfg.colorbar                = 'yes';
%     cfg.shading                 = 'interp';
%     figure,ft_topoplotER(cfg, gammaFreq_bsl);
%     saveas(gcf,[proc.save_folder,'\',b,'-gammaTopo'],'png')

    %
%     cfg                         = [];
%     cfg.zlim                    = 'maxmin';
%     cfg.channel                 = {'MEG0432+0433', 'MEG0442+0443'}
%     cfg.xlim                    = [-0.09 0.1]; 
%     cfg.maskstyle               = 'saturation';	
%     cfg.masknans                = 'yes';
%     cfg.layout                  = 'neuromag306cmb.lay';
%     cfg.axes                    = 'no';
%     cfg.colorbar                = 'yes';
%     cfg.shading                 = 'interp';
%     figure,ft_singleplotER(cfg, gammaFreq_bsl);
%     saveas(gcf,[proc.save_folder,'\', b,'-gammaEvolution'],'tif')

%     % select data
%     cfg = [];
%     cfg.channel                 = {'MEG0432+0433', 'MEG0442+0443'}
%     cfg.avgoverchan             = 'no';
%     cfg.latency                 = [0.02 0.06];
%     data_gamma_trials           = ft_selectdata(cfg, gammaFreq_bsl);
%     save ([proc.save_folder,'\',b, '-GammaRptChan'], 'data_gamma_trials', '-v7.3');
% %     
%     cfg = [];
%     cfg.avgoverrpt              = 'yes';
%     data_gamma                  = ft_selectdata(cfg, data_gamma_trials);
%     gamma_pow_max(loop)         = max(max(data_gamma.powspctrm));
    
%     close all;

%% Add function to save images per condition


%        Check for files
%                 disp(subjName{subjLoop})
%                 disp(fileList{condLoop})
%                 disp(num2str(condName{condLoop}))
%                 exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'-SaveData.mat'],'file')