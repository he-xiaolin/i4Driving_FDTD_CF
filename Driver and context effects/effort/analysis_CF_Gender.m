%% Effort  – Male vs Female (rescaled (within participant) + trimmed) ───
clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ─────────────── trimming settings (same as your “FAST” script) ──────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central band (与前面脚本保持一致)
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ─────────────── collectors ──────────────────────────────────────────────
all_eff_M = [];  all_pid_M = [];
all_eff_F = [];  all_pid_F = [];

% (optional) list of indices to skip
excluded = [];

for k = 1:numel(participants_CF)
    if ismember(k,excluded),  continue,  end
    
    % raw signals ---------------------------------------------------------
    THW = participants_CF(k).THW(:);
    Eff = participants_CF(k).effort_scaled(:);   % ← still 0-100 in the .mat
    
    % 1) rescale to 0-1  (driver-wise)
    Eff_n = local_normalize01(Eff);

    % 2) build trimming mask (same helper as in your scatter code)
    keep  = local_buildKeepMask(THW,Eff_n,trim);

    % 3) keep only the central band, *then* store by condition
    Eff_keep = Eff_n(keep);

    if strcmpi(participants_CF(k).genderGroup,'M')      % Male
        all_eff_M = [all_eff_M ; Eff_keep];
        all_pid_M = [all_pid_M ; repmat(k,numel(Eff_keep),1)];

    elseif strcmpi(participants_CF(k).genderGroup,'F')  % Female
        all_eff_F = [all_eff_F ; Eff_keep];
        all_pid_F = [all_pid_F ; repmat(k,numel(Eff_keep),1)];
    end
end

%% Participant-level mean Effort (already 0-1 & trimmed) ───────────────────
meanEff_M = splitapply(@mean, all_eff_M, findgroups(all_pid_M));
meanEff_F = splitapply(@mean, all_eff_F, findgroups(all_pid_F));

%% Box-plot + statistics ──────────────────────────────────────────────────
figure
boxplot([meanEff_M; meanEff_F], ...
        [repelem({'Male'},numel(meanEff_M)) ...
         repelem({'Female'},numel(meanEff_F))]', ...
        'Symbol','')
ylabel('Effort')
title('Effort  – Male vs Female')
set(gca,'XLim',[0.5 2.5])

% two-sample t-test  (driver means)
[H,p,CI,stats] = ttest2(meanEff_M,meanEff_F);
std_M = std(meanEff_M); 
std_F = std(meanEff_F);

fprintf('\nMale mean = %.3f (SD = %.4f)   |  Female mean = %.3f (SD = %.4f)\n', ...
        mean(meanEff_M), std_M, mean(meanEff_F), std_F)
fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

%% save
if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
exportgraphics(gcf,'./Figure/Effort_Gender_Male_vs_Female_trimmed.png','Resolution',300)

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
