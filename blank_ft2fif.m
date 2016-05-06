function blank_ft2fif(data_folder, stim_chan, pre_stim_time, post_stim_time, save_folder)

% blank_ft2fif - Blanking artefacts in an MEG dataset(.fif)
%
% Takes the following arguments: 
% data_folder     = Folder in which data which is to be blanked is stored
% stim_chan       = Name of the trigger channel, for detecting events (default: 'STI001')
% pre_stim_time   = Time before stimulation to be blanked
%                      (must be a negative integer), in ms (default: 2 ms)
% post_stim_time  = Time after stimulation to be blanked, in ms (default: 7 ms)
% save_folder     = Folder to which blanked files are written  (default: [data_folder\blanked]) 
%
% Written on 02-12-2013

if nargin<1
    error('Too few arguments')
elseif nargin<2
    stim_chan = 'STI001';
elseif nargin<3
    pre_stim_time = -2;
elseif nargin<4
    post_stim_time = 7;  
elseif nargin<5
    cd(data_folder);
    mkdir('blanked');
    save_folder = [data_folder, '\blanked'];
end
% 
cd(data_folder);
% ft_defaults;
% 
% % Extract filenames from data folder
% filename = dir('*.fif');
% for loop = 1:length(filename);
%     fprintf('\n'); disp(['##########   ', num2str(loop), '    ###############'])
%     fprintf('\n'); disp(filename(loop).name)
% 
%     proc = [];
%     % Import file
%     cfg = [];
%     cfg.datafile            = filename(loop).name;
%     cfg.headerfile          = filename(loop).name;
%     cfg.continuous          = 'yes';
%     [proc.data]             = ft_preprocessing(cfg);
% 
%     % Extract events
%     cfg.trialdef.eventtype  = stim_chan; 
%     cfg.trialdef.eventvalue = 5;
%     cfg.dataset             = filename(loop).name;
%     proc.data_trial         = ft_definetrial(cfg);
% 
%     % Build Blanking sequence
%     proc.blank_seq          = ones(size(proc.data.trial{1},2),1);
%     for ind = pre_stim_time:post_stim_time
%         proc.blank_seq(proc.data_trial.trl(:,1)+ind) = 0;
%     end
% 
%     proc.MEG_channel = find(strncmp(proc.data.label, 'MEG',3)==1);                    % Pick only MEG channels
%     proc.data.trial{1}(proc.MEG_channel,:) = ...
%             bsxfun(@times, proc.data.trial{1}(proc.MEG_channel,:), proc.blank_seq');  % Multiply blanking sequence with MEG channels

    % Write data back to .fif format
    fieldtrip2fiff([save_folder,'\', '015_06_try'],proc.preproc_data_MEG);
    disp(['File written:',save_folder,'\', filename(loop).name]);
% end

return