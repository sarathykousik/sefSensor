%% Calculate Induced, evoked gamma
% Use cleanData.m
% Import, epoch, filter, artefact data
% Save data in patients folders/conditions

clc; clear all; close all;

%%
filePath    = 'J:\MEG_Research\SEF\tsssNEW';
save_folder = 'J:\MEG_Research\SEF\tsssNEW';
mkdir(save_folder)
cd(filePath)
filenames      = dir('*.mat');
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')

for loop=4:length(filenames)
    [a b c]     =   fileparts(filenames(loop).name);    
    disp(['#########  ',b])
    load(filenames(loop).name)
   % Calc ERF
%     cfg = []; 
%     tlck = ft_timelockanalysis(cfg, visClean)
%     ind_dat=[];
   % Calc Ind data (visClean-ERF)
%     for trialLoop = 1: length(visClean.trial)
%         ind_dat{1,trialLoop} = ...
%             visClean.trial{1,trialLoop} - tlck.avg;
%     end

%   visCleanInd=visClean;
%   visCleanInd.trial = ind_dat;
   
  % Calc Ind gamma
   gammaFreq = calcGamma(saveData.dataArtRej); %calcGamma(visClean);%
  
   % Calc Ind gamma
%    indGamma = calcGamma(visCleanInd);
   
   % Calc Evoked gamma
%    evokedGamma = calcGamma(tlck);
   
   save ([save_folder, '/',b, '-gammaFreq'], 'gammaFreq');
%    save ([save_folder, '/',b, '-indGamma'], 'indGamma');
%    save ([save_folder, '/',b, '-evokedGamma'], 'evokedGamma');
end







