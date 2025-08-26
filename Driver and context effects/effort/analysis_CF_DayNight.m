%% Effort  – Day  vs Night   (rescaled (within participant) + trimmed)  ────────────────────────
clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ─────────────── trimming settings (same as your “FAST” script) ──────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central 50 %  (≈ ±0.674 σ)
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ─────────────── collectors ──────────────────────────────────────────────
all_eff_Day   = [];   all_pid_Day   = [];
all_eff_Night = [];   all_pid_Night = [];

% (optional) list of indices to skip
excluded = [];                    

for k = 1:numel(participants_CF)
    if ismember(k,excluded),  continue,  end
    
    % raw signals ---------------------------------------------------------
    THW   = participants_CF(k).THW(:);
    Eff   = participants_CF(k).effort_scaled(:);  % ← still 0-100 in the .mat
    
    % 1) rescale to 0-1  (driver-wise)
    Eff_n = local_normalize01(Eff);

    % 2) build trimming mask (same helper as in your scatter code)
    keep  = local_buildKeepMask(THW,Eff_n,trim);

    % 3) keep only the central 50 %, *then* store by condition
    Eff_keep = Eff_n(keep);

    if contains(participants_CF(k).filename,'i4Driving_D')      % Day
        all_eff_Day   = [all_eff_Day ;   Eff_keep];
        all_pid_Day   = [all_pid_Day ;   repmat(k,numel(Eff_keep),1)];
        
    elseif contains(participants_CF(k).filename,'i4Driving_N')  % Night
        all_eff_Night = [all_eff_Night ; Eff_keep];
        all_pid_Night = [all_pid_Night ; repmat(k,numel(Eff_keep),1)];
    end
end

%% Participant-level mean Effort (already 0-1 & trimmed) ───────────────────
meanEff_Day   = splitapply(@mean, all_eff_Day  , findgroups(all_pid_Day  ));
meanEff_Night = splitapply(@mean, all_eff_Night, findgroups(all_pid_Night));

%% Box-plot + statistics ──────────────────────────────────────────────────
figure
boxplot([meanEff_Day; meanEff_Night], ...
        [repelem({'Day'},numel(meanEff_Day)) ...
         repelem({'Night'},numel(meanEff_Night))]', ...
        'Symbol','')
ylabel('Efforttest')
title('Effort  – Day vs Night')
set(gca,'XLim',[0.5 2.5])

% two-sample t-test  (driver means)
[H,p,CI,stats] = ttest2(meanEff_Day,meanEff_Night);
std_Day   = std(meanEff_Day);
std_Night = std(meanEff_Night);

fprintf('\nDay  mean = %.3f (SD = %.4f)   |  Night mean = %.3f (SD = %.4f)\n', ...
        mean(meanEff_Day), std_Day, mean(meanEff_Night), std_Night)

fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

%% save
if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
exportgraphics(gcf,'./Figure/Effort_Day_vs_Night_trimmed.png','Resolution',300)

% ──────────────────────────  helper functions  ───────────────────────────
function keep = local_buildKeepMask(x,y,t)
    % identical copy of your FAST version (shortened header)
    keep = false(size(x));
    if strcmpi(t.method,'sd'), k = sqrt(2)*erfinv(t.centralPct/100); end
    for b = 1:numel(t.edgesTHW)-1
        in = x>=t.edgesTHW(b) & x<t.edgesTHW(b+1);
        if nnz(in)<t.minPts, keep(in)=true; continue, end
        switch lower(t.method)
            case 'percentile'
                tail=(100-t.centralPct)/2;
                lo=prctile(y(in),tail); hi=prctile(y(in),100-tail);
                keep(in)=y(in)>=lo & y(in)<=hi;
            case 'sd'
                mu=mean(y(in)); sd=std(y(in),0);
                keep(in)=abs(y(in)-mu)<=k*sd;
        end
    end
end

function y_n = local_normalize01(y)
    d = max(y)-min(y);
    if d>0,   y_n=(y-min(y))./d;      % map to [0,1]
    else      y_n=0.5*ones(size(y));  % constant vector→mid-level
    end
end
