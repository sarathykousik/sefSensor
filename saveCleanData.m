%% saveCleanData
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc;    clear all;     close all;
ft_defaults;
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')
proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\SEF\SEF-TSSS\013_7UF';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              =  10;
proc.save_folder                 = 'J:\MEG_Research\SEF\SEFcleanData';
mkdir(proc.save_folder)

%% Files

[temp.num,temp.txt,temp.raw] = xlsread('J:\MEG_Research\SEF\SEF-TSSS\subjectDataList.xls');
temp.num_cell = num2cell(temp.num);

subjName = {temp.raw{1,5:10}};
condName = {temp.raw{4:10,1}};
% subjName = '009_3KE';
% condName = {'30', '60', '120','1','150', '180', '210' }

%% Import, epoch, filter, artefact data

% for subjLoop = 1:length(subjName)
subjLoop = 5
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
        for condLoop  = 7%[1,2,4:length(fileList)]
            disp('#######################################')
            disp(['   ****    ',num2str(condName{condLoop}),'     ****    '])
            disp('#######################################')
            
%             fif_file_preproc([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'],...
%                         subjName{subjLoop}, num2str(condName{condLoop}),proc);
                    
            fif_file_preproc([fileList{condLoop},'.fif'],...
                        subjName, num2str(condName{condLoop}),proc);
%        Check for files
%                 disp(subjName{subjLoop})
%                 disp(fileList{condLoop})
%                 disp(num2str(condName{condLoop}))
%                 checkFlag(subjLoop, condLoop) = ...
%                     exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'],'file');
        end
% end




