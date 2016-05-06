%% Gamma evolution calc
% Use cleanData-VisClean

clc;    clear all;     close all;

ft_defaults;

proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEFVisClean';  
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFGammanewtSSS';
mkdir(proc.save_folder)

%
cd(proc.data_folder)
filenames      = dir('*.mat');

for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)

    % Gamma mtp: Fc = 65 Hz, BW: 50 Hz
    cfg              = [];
    cfg.paramter     = 'trial';
    cfg.taper        = 'dpss';
    cfg.keeptrials   = 'yes';
    cfg.output       = 'pow';
    cfg.channel      = 'MEGGRAD';
    cfg.method       = 'mtmconvol';
    cfg.pad          = 2;
    cfg.foi          = 65;                          
    cfg.t_ftimwin    = 7./cfg.foi;   
    cfg.toi          = -0.500:0.010:0.500;           
    cfg.tapsmofrq    = 25;
    gammaFreq        = ft_freqanalysis(cfg, visClean);

    save([proc.save_folder,'\freq\' , b, '-gammaFreq'], 'gammaFreq', '-v7.3');

%     cfg = [];
%     cfg.baseline                = [-0.15 -0.10];
%     cfg.baselinetype            = 'relchange';
%     cfg.parameter               = 'powspctrm';
%     gammaFreq_bsl               = ft_freqbaseline(cfg, gammaFreq);
% 
%     cfg                         = [];
%     gammaCmb               = ft_combineplanar(cfg, gammaDesc);
% %     save ([proc.save_folder, '\bsl\',b, '-gammaBsl'], 'gammaFreq_bsl', '-v7.3');
% %     
%     cfg             = [];
%     cfg.jackknife   = 'yes';
% %     cfg.biascorrect = 'yes';
% %     cfg.variance = 'yes';
%     gammaDesc       = ft_freqdescriptives(cfg, gammaFreq_bsl);
%     save ([proc.save_folder,'\desc\', b, '-gammaDesc'], 'gammaDesc', '-v7.3');

        disp('#######################################')
    
end

%%
    cfg                         = [];
    cfg.xlim                    = [0.03 0.08]; 
    cfg.channel                 = {'all', '-MEG1642+1643', '-MEG1612+1613',...
        '-MEG0742+0743', '-MEG1132+1133','-MEG0442+0443', '-MEG1522+1523', '-MEG2042+2043'};
    cfg.ylim                    = 'maxmin';
    cfg.maskstyle               = 'saturation';	
    cfg.masknans                = 'yes';
    cfg.layout                  = 'neuromag306cmb.lay';
    cfg.shading                 = 'interp';
    cfg.axes                    = 'no';
    ft_topoplotER(cfg, gammaCmb);
    
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

%%
figure,
% plot(proc.data_import.time{1}, proc.data_import.trial{1}(324,:),'b'), hold on
plot(proc.data_import.time{1}, proc.blank_seq' .* proc.data_import.trial{1}(324,:),'r')


%% Add function to save images per condition


%        Check for files
%                 disp(subjName{subjLoop})
%                 disp(fileList{condLoop})
%                 disp(num2str(condName{condLoop}))
%                 exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'-SaveData.mat'],'file')