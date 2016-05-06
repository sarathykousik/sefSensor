
%% Gamma
total_data = [gammaC_trials; gammaPD_trials];
% total_trial = 0;

fid = fopen('gammaVal-newtSSS.txt','w');
fprintf(fid,['Value, Trial, Condition, Group, Subj, Area \n']);
total_trial = 1;
% gammaVal = zeros(79267,6);

for subjLoop = 1:20
    
    if(subjLoop<=10)
        group = 1;  % Control
    else
        group = 2;  % PD
    end
    
    for condLoop = 1:7
%        if(condLoop==1)
%            area = 1;
%        elseif(condLoop<6 & condLoop>1)
%            area = 2;
%        elseif(condLoop==7)
%            area = 3;
%        end
       
       trial_no = numel(total_data{subjLoop,condLoop}.powspctrm);
       
       for trialLoop = 1:trial_no
            
            gammaVal(total_trial+trialLoop,:,:,:,:) = ...
                [total_data{subjLoop, condLoop}.powspctrm(trialLoop) trialLoop condLoop group subjLoop];
            
            fprintf(fid,'%f,%d,%d,%d,%d \n',...
                [total_data{subjLoop, condLoop}.powspctrm(trialLoop) trialLoop condLoop group subjLoop]);
            
        end
        total_trial = total_trial+trial_no;
        
        
    end
end

fclose(fid);

%% ERF Amp N20m
ControlERF  = ControlN20mAmp./1e-12.*10;
PDERF       = maxPosN20Amp./1e-12.*10;
total_data = [ControlERF; PDERF];

% total_trial = 0;
ERF = [];
fid = fopen('ERF.txt','w');
fprintf(fid,['Value, Condition, Group, Subj, Area \n']);
total_trial = 1;
% gammaVal = zeros(79267,6);

for subjLoop = 1:20
    
    if(subjLoop<=10)
        group = 1;  % Control
    else
        group = 2;  % PD
    end
    
    for condLoop = 1:7
       if(condLoop==1)
           area = 1;
       elseif(condLoop<6 & condLoop>1)
           area = 2;
       elseif(condLoop==7)
           area = 3;
       end
       
%        trial_no = numel(total_data{subjLoop,condLoop}.powspctrm);
       
%        for trialLoop = 1:trial_no
            
            ERF(total_trial,:,:,:,:) = ...
                [total_data(subjLoop, condLoop) condLoop group subjLoop area];
            
            fprintf(fid,'%f,%d,%d,%d,%d \n',...
                 [total_data(subjLoop, condLoop) condLoop group subjLoop area]);
            
%         end
%         total_trial = total_trial+condLoop;
            total_trial = size(ERF,1)+1;
%             disp([total_trial,subjLoop ,condLoop])
        
        
    end
end

fclose(fid);

%% UPDRS
total_data = [updrs(:,[1 2 3 7])];

% total_trial = 0;
UPDRS = [];
fid = fopen('updrs.txt','w');
fprintf(fid,['Value, Condition, Subj, Area \n']);
total_trial = 1;
% gammaVal = zeros(79267,6);

for subjLoop = 1:10
    
       for condLoop = 1:4
%            if(condLoop==1)
%                area = 1;
%            elseif(condLoop<6 & condLoop>1)
%                area = 2;
%            elseif(condLoop==7)
%                area = 3;
%        end
       
%        trial_no = numel(total_data{subjLoop,condLoop}.powspctrm);
       
%        for trialLoop = 1:trial_no
            
            UPDRS(total_trial,:,:,:,:) = ...
                [total_data(subjLoop, condLoop) condLoop subjLoop area];
            
            fprintf(fid,'%f,%d,%d,%d \n',...
                 [total_data(subjLoop, condLoop) condLoop subjLoop area]);
            
%         end
%         total_trial = total_trial+condLoop;
            total_trial = size(UPDRS,1)+1
%             disp([total_trial,subjLoop ,condLoop])
        
        
    end
end

fclose(fid);

