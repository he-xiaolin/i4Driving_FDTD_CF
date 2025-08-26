%% Effort – NDRT vs NoNDRT (rescaled (within participant) + trimmed) ───
clc; clear; close all
load participants_CF.mat

% ─────────────── trimming settings ──────────
trim.method      = 'sd';
trim.centralPct  = 95;
trim.edgesTHW    = 0 : 0.25 : 10;
trim.minPts      = 5;

% ─────────────── collectors ──────────
all_eff_NDRT   = [];   all_pid_NDRT   = [];
all_eff_noNDRT = [];   all_pid_noNDRT = [];

% 判断是否 NDRT 的函数
isNDRTfile = @(fname) (numel(fname) >= 13 && fname(13) == '1');

% 可选剔除名单
excluded = [];

for k = 1:numel(participants_CF)
    if ismember(k,excluded), continue, end

    THW = participants_CF(k).THW(:);
    Eff = participants_CF(k).effort_scaled(:);

    % 归一化到 0–1
    Eff_n = local_normalize01(Eff);

    % 修剪
    keep = local_buildKeepMask(THW, Eff_n, trim);
    Eff_keep = Eff_n(keep);

    % 分组
    if isNDRTfile(participants_CF(k).filename)
        all_eff_NDRT   = [all_eff_NDRT; Eff_keep];
        all_pid_NDRT   = [all_pid_NDRT; repmat(k,numel(Eff_keep),1)];
    else
        all_eff_noNDRT = [all_eff_noNDRT; Eff_keep];
        all_pid_noNDRT = [all_pid_noNDRT; repmat(k,numel(Eff_keep),1)];
    end
end

%% Participant-level mean Effort
meanEff_NDRT   = splitapply(@mean, all_eff_NDRT,   findgroups(all_pid_NDRT));
meanEff_noNDRT = splitapply(@mean, all_eff_noNDRT, findgroups(all_pid_noNDRT));

%% Box-plot + statistics
figure
boxplot([meanEff_NDRT; meanEff_noNDRT], ...
        [repelem({'NDRT'},numel(meanEff_NDRT)) ...
         repelem({'No NDRT'},numel(meanEff_noNDRT))]', ...
        'Symbol','')
ylabel('Effort')
title('Effort – NDRT vs NoNDRT')
set(gca,'XLim',[0.5 2.5])

% t-test
std_NDRT   = std(meanEff_NDRT);
std_noNDRT = std(meanEff_noNDRT);
[H,p,CI,stats] = ttest2(meanEff_NDRT, meanEff_noNDRT);

fprintf('\nNDRT mean = %.3f (SD = %.4f)   |  NoNDRT mean = %.3f (SD = %.4f)\n', ...
        mean(meanEff_NDRT), std_NDRT, mean(meanEff_noNDRT), std_noNDRT)
fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

%% Save
if ~exist('./Figure','dir'), mkdir('./Figure'); end
exportgraphics(gcf,'./Figure/Effort_NDRT_vs_NoNDRT_trimmed.png','Resolution',300)

% ─────────────── helpers ──────────
function keep = local_buildKeepMask(x,y,t)
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
    if d>0,   y_n=(y-min(y))./d;
    else      y_n=0.5*ones(size(y));
    end
end
