clc;     clear all;    close all;


%%


cfg = [];
cfg.latency = [0.045 0.055];  % specify latency window around M50 peak
cfg.numdipoles = 1;
cfg.hdmfile = 'bauer_m.hdm';
cfg.feedback = 'textbar';
cfg.grid.resolution = 2;
cfg.grid.unit = 'cm';
dipM50 = ft_dipolefitting(cfg, avg);
cfg.latency = [0.100 0.120]; % specify latency window around M100 peak
dipM100 = ft_dipolefitting(cfg, avg);

