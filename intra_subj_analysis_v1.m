%% Code for inter-subject analysis
% Written on 06-02-2014

%%
clc;       clear all;       close all;

ft_defaults;
% Get data blanked
proc.data_folder                 = ...
    'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -2;
proc.post_stim_time              = 7;
proc.save_folder                 = 'J:\results_FT';
subjList = {'014_NJO','013_7UF', '012_KBU', '011_EHI', '010_XKG', '009_3KE'...
            '008_ONB', '007_FUE', '006_TYB'};
%% 

% proc.blanked_folder         = proc.save_folder;                                  %'C:\Program Files\MATLAB\R2012b\toolbox_add_on\data-tsss-sef\014_NJO\blanked';
cd(proc.data_folder)
filename                    = []; 
filename                    = dir('*.fif');
% loop                        = 2;

% parameters to pass: file_name, proc 

for loop = 1:length(filename)
   
    fif_file_preproc(filename(loop).name, proc);
    
end


