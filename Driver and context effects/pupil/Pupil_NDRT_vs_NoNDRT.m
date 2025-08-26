%% Pupil diameter – NDRT vs NoNDRT (rescaled + trimmed)
clc; clear; close all
load participants_CF.mat

trim.method='sd'; trim.centralPct=95; trim.edgesTHW=0:0.25:10; trim.minPts=5;

all_work_NDRT=[];   all_pid_NDRT=[];
all_work_noNDRT=[]; all_pid_noNDRT=[];

isNDRTfile = @(fname) (numel(fname)>=13 && fname(13)=='1');
excluded=[];

for k=1:numel(participants_CF)
    if ismember(k,excluded), continue, end
    THW=participants_CF(k).THW(:);
    thisWork=participants_CF(k).pupil_diameter_scaled(:);

    Work_n=local_normalize01(thisWork);
    keep=local_buildKeepMask(THW,Work_n,trim);
    Work_keep=Work_n(keep);

    if isNDRTfile(participants_CF(k).filename)
        all_work_NDRT=[all_work_NDRT;Work_keep];
        all_pid_NDRT=[all_pid_NDRT;repmat(k,numel(Work_keep),1)];
    else
        all_work_noNDRT=[all_work_noNDRT;Work_keep];
        all_pid_noNDRT=[all_pid_noNDRT;repmat(k,numel(Work_keep),1)];
    end
end

meanPD_NDRT  = splitapply(@mean,all_work_NDRT ,findgroups(all_pid_NDRT));
meanPD_noNDRT= splitapply(@mean,all_work_noNDRT,findgroups(all_pid_noNDRT));

figure
boxplot([meanPD_NDRT; meanPD_noNDRT], ...
        [repelem({'NDRT'},numel(meanPD_NDRT))' ; ...
         repelem({'No NDRT'},numel(meanPD_noNDRT))' ], ...
        'Symbol','')
ylabel('Pupil diameter')
title('Pupil diameter – NDRT vs NoNDRT')
set(gca,'XLim',[0.5 2.5])

sd_NDRT=std(meanPD_NDRT); sd_noNDRT=std(meanPD_noNDRT);
[H,p,CI,stats]=ttest2(meanPD_NDRT,meanPD_noNDRT);
fprintf('\nNDRT   mean = %.3f (SD = %.4f)  |  NoNDRT mean = %.3f (SD = %.4f)\n', ...
        mean(meanPD_NDRT),sd_NDRT,mean(meanPD_noNDRT),sd_noNDRT)
fprintf('t(%d) = %.2f,  p = %.4f\n\n',stats.df,stats.tstat,p)

if ~exist('./Figure','dir'), mkdir('./Figure'); end
exportgraphics(gcf,'./Figure/Pupil_NDRT_vs_NoNDRT_trimmed.png','Resolution',300)

% helpers
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
