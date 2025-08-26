%% Pupil diameter – Day vs Night (rescaled (within participant) + trimmed)
clc; clear; close all
load participants_CF.mat

% ─ trimming settings ─
trim.method      = 'sd';
trim.centralPct  = 95;
trim.edgesTHW    = 0 : 0.25 : 10;
trim.minPts      = 5;

% ─ collectors ─
all_work_Day   = []; all_pid_Day   = [];
all_work_Night = []; all_pid_Night = [];

excluded = [];

for k = 1:numel(participants_CF)
    if ismember(k,excluded), continue, end

    THW      = participants_CF(k).THW(:);
    thisWork = participants_CF(k).pupil_diameter_scaled(:);

    Work_n = local_normalize01(thisWork);
    keep   = local_buildKeepMask(THW, Work_n, trim);
    Work_keep = Work_n(keep);

    if contains(participants_CF(k).filename,'i4Driving_D')      % Day
        all_work_Day   = [all_work_Day ; Work_keep];
        all_pid_Day    = [all_pid_Day  ; repmat(k,numel(Work_keep),1)];
    elseif contains(participants_CF(k).filename,'i4Driving_N')  % Night
        all_work_Night = [all_work_Night ; Work_keep];
        all_pid_Night  = [all_pid_Night  ; repmat(k,numel(Work_keep),1)];
    end
end

% participant-level mean
meanPD_Day   = splitapply(@mean, all_work_Day  , findgroups(all_pid_Day));
meanPD_Night = splitapply(@mean, all_work_Night, findgroups(all_pid_Night));

% boxplot + stats
figure
boxplot([meanPD_Day; meanPD_Night], ...
        [repelem({'Day'},numel(meanPD_Day))' ; ...
         repelem({'Night'},numel(meanPD_Night))' ], ...
        'Symbol','')
ylabel('Pupil diameter')
title('Pupil diameter – Day vs Night')
set(gca,'XLim',[0.5 2.5])

sd_Day   = std(meanPD_Day);
sd_Night = std(meanPD_Night);
[H,p,CI,stats] = ttest2(meanPD_Day, meanPD_Night);

fprintf('\nDay   mean = %.3f (SD = %.4f)  |  Night mean = %.3f (SD = %.4f)\n', ...
        mean(meanPD_Day), sd_Day, mean(meanPD_Night), sd_Night)
fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

% save
if ~exist('./Figure','dir'), mkdir('./Figure'); end
exportgraphics(gcf,'./Figure/Pupil_Day_vs_Night_trimmed.png','Resolution',300)

% ──────────────────────────  helper functions  ───────────────────────────
function keep = local_buildKeepMask(x,y,t)
    keep = false(size(x));
    if strcmpi(t.method,'sd'), k = sqrt(2)*erfinv(t.centralPct/100); end
    for b = 1:numel(t.edgesTHW)-1
        in = x>=t.edgesTHW(b) & x<t.edgesTHW(b+1);
        if nnz(in)<t.minPts, keep(in)=true; continue, end
        switch lower(t.method)
            case 'percentile'
                tail=(100-t.centralPct)/2; 
                lo=prctile(y(in),tail); 
                hi=prctile(y(in),100-tail);
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

