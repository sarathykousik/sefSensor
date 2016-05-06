clc;       clear all;       close all;

ft_defaults;
proc                             = [];
proc.data_folder                 = 'J:\data-tsss-sef\SEF-TSSS';
proc.stim_chan                   = 'STI001';
proc.pre_stim_time               = -5;
proc.post_stim_time              =  10;
proc.save_folder                 = 'J:\procFtResults_trials';
mkdir(proc.save_folder)

%%

[temp.num,temp.txt,temp.raw] = xlsread('J:\data-tsss-sef\SEF-TSSS\subjectDataList.xls');
temp.num_cell = num2cell(temp.num);

subjName = {temp.raw{1,6:10}};
condName = {temp.raw{4:10,1}};

%% Initial processing
for subjLoop = 1:length(subjName)
    %Choose filenames
    [row, col]=find(strcmp(temp.raw, subjName{subjLoop})==1);
    disp('==========================================')
    disp('==========================================')
    disp(['   *******      ',subjName{subjLoop}, '    *********'])
    disp('==========================================')
    disp('==========================================')
    mkdir([proc.save_folder,'\',subjName{subjLoop}]);

%     fileList = dir('*.fif');
    fileList = {temp.raw{4:10,col}}; 
    cd([proc.data_folder,'\', subjName{subjLoop}])
    
        for condLoop  = 1:length(fileList)
            disp('#######################################')
            disp(['   ****    ',num2str(condName{condLoop}),'     ****    '])
            disp('#######################################')
            
            fif_file_preproc([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'],...
                        subjName{subjLoop}, num2str(condName{condLoop}),proc);

%                 subjName{subjLoop}
%                 fileList{condLoop}
%                 num2str(condName{condLoop})
%                 exist([proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif'],'file')
        end
end

%% Check

% 
% file_name = [proc.data_folder,'\',subjName{subjLoop},'\',fileList{condLoop},'.fif']
%  subj=                       subjName{subjLoop}
% condN = num2str(condName{condLoop})