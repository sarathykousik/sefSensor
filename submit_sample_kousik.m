%% Pseudo-generic sample script for submitting jobs to the cluster

% [2014-01-15: setting it up as a generic sample script]

% clusterconfig('wait',1);                            % waiting for the job to finish before moving on to the next line in the code
% % NB! locks the cmd-line, so don't do it

% clusterconfig('scheduler','none'); % Run everything without any clusterizing code at all
% clusterconfig('scheduler','local'); % Run everything clusterized but only to the machine you're working on (NB! NOT Hyades01!)
% clusterconfig('scheduler','cluster'); % Run everything truly clusterized ('default')

% clusterconfig('long_running',1); % for jobs with a duration > 2 hrs


%% INPUT for MY_FUNCTION

clear all

% input-structure with keys (i.e. numbers) for conditions and/or groups -
% these should thus refer to some sort of key in my_function, so that
% input(1) selects the relevant group to analyze and input(2) selects the
% relevant condition to analyze

conds = 1:6;
groups = 1:2;
for i = length(groups):-1:1
    for j = length(conds):-1:1
        input{j+(i-1)*length(conds)} = [groups(i) conds(j)];    % has to be cell array for the job2cluster submission to accept it
    end
end


%% ALTERNATIVE INPUT structure tailored at Kousik's ITC cluster code

filePath = '/projects/yours&mine/raw/';     % whatever dir you're sitting in when running the cluster-code you sent me

cd(filePath)

% Pick all filenames
filenames = dir('*.mat');

for i = length(filenames):-1:1
    input{i} = i;
end

% now instead of looping over 
% for loop  = 1:length(filenames)

% you just call
% [a b c] = fileparts(filenames(input).name);

% in your function, and you have to specify and cd to the filePath in that
% function, as well as call the dir-command as you're already doing
% Hope this makes sense


%% Actually submit the job

display(input(:))

jobid = job2cluster(@my_function, input);


% jobstatus(jobid)


% % example function
% function out = my_own_function(patient)
%   out = [];
%   filedir = fullfile('/my/file/location', patient);
%  spm_something(filedir);
 
% clusterconfig('scheduler', [])     % [] = no clusterization; 'local' = clusterized on the local machine; 'cluster' = truly clusterized
% clusterconfig('wait', 1);          % in order for results = jobresults(job_id) not to be submitted until the jub has actually finished
%  
%  
%  
% e=joberrors(job_id);
% o=joboutput(job_id);

% % kill a job
% jobdestroy(job_id);

% % who's running what?
% clusterjobs();

