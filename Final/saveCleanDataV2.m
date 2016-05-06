%% saveCleanData
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc;    clear all;     close all;
ft_defaults;
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')
proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEF-TSSS\Control\028_D4Q';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              =  10;
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFControlcleanData';
mkdir(proc.save_folder)

%% Files

% [temp.num,temp.txt,temp.raw] = xlsread('J:\MEG_Research\SEF\SEF-TSSS\subjectDataList.xls');
% temp.num_cell = num2cell(temp.num);

% subjName = {temp.raw{1,5:10}};
% condName = {temp.raw{4:10,1}};
% subjName = '009_3KE';
% condName = {'30', '60', '120','1','150', '180', '210' }

% Import, epoch, filter, artefact data
cd(proc.data_folder)
filenames      = dir('*.fif');


for condLoop  = 1:length(filenames)%[1,2,4:length(fileList)]
            [a b c]     =   fileparts(filenames(condLoop).name);
            disp(['   ****    ',b,'     ****    '])

            
%             fif_file_preproc([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'],...
%                         subjName{subjLoop}, num2str(condName{condLoop}),proc);
            tok_name = tokenize(b, '_');
            subjName= [tok_name{1},'_',tok_name{2}];
            fif_file_preproc(filenames(condLoop).name,...
                        subjName, num2str(condLoop),proc);
%        Check for files
%                 checkFlag(subjLoop, condLoop) = ...
%                     exist(filenames(condLoop).name, 'file')

end





