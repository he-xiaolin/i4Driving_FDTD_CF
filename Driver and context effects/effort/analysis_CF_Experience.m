%% Effort  – Low vs High Driving Experience (rescaled (within participant) + trimmed) ───
clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ─────────────── trimming settings (same as your “FAST” script) ──────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central 50 %  (≈ ±0.674 σ)  ← 保持与上面脚本一致
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ─────────────── collectors ──────────────────────────────────────────────
all_eff_L   = [];   all_pid_L   = [];
all_eff_H   = [];   all_pid_H   = [];

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

    % 3) keep only the central band, *then* store by condition
    Eff_keep = Eff_n(keep);

    if strcmpi(participants_CF(k).drivingExperience,'L')      % Low experience
        all_eff_L = [all_eff_L ; Eff_keep];
        all_pid_L = [all_pid_L ; repmat(k,numel(Eff_keep),1)];

    elseif strcmpi(participants_CF(k).drivingExperience,'H')  % High experience
        all_eff_H = [all_eff_H ; Eff_keep];
        all_pid_H = [all_pid_H ; repmat(k,numel(Eff_keep),1)];
    end
end

%% Participant-level mean Effort (already 0-1 & trimmed) ───────────────────
meanEff_L = splitapply(@mean, all_eff_L, findgroups(all_pid_L));
meanEff_H = splitapply(@mean, all_eff_H, findgroups(all_pid_H));

%% Box-plot + statistics ──────────────────────────────────────────────────
figure
boxplot([meanEff_L; meanEff_H], ...
        [repelem({'Low'},numel(meanEff_L)) ...
         repelem({'High'},numel(meanEff_H))]', ...
        'Symbol','')
ylabel('Effort')
title('Effort  – Low vs High Driving Experience')
set(gca,'XLim',[0.5 2.5])

% two-sample t-test  (driver means)
[H,p,CI,stats] = ttest2(meanEff_L,meanEff_H);
std_L = std(meanEff_L);
std_H = std(meanEff_H);

fprintf('\nLow  mean = %.3f (SD = %.4f)   |  High mean = %.3f (SD = %.4f)\n', ...
        mean(meanEff_L), std_L, mean(meanEff_H), std_H)

fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

%% save
if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
exportgraphics(gcf,'./Figure/Effort_Experience_Low_vs_High_trimmed.png','Resolution',300)

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
