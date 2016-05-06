
clc;     clear all;     close all
filenames      = dir('*.mat');


for loop  = 1:length(filenames)
    disp('#######################################')
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['   ****    ',b,'     ****    '])
    load(filenames(loop).name)
    
%     tap(loop)    = gammaFreq.cfg.taper;
    padlen(loop) = gammaFreq.cfg.pad;
    foicfg(loop) = gammaFreq.cfg.foi;
    tapsmo(loop) = gammaFreq.cfg.tapsmofrq;
    tfwin(loop)  = gammaFreq.cfg.t_ftimwin;
    
end