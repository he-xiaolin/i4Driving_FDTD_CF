%% Pupil diameter – Low vs High Driving Experience (rescaled + trimmed)
clc; clear; close all
load participants_CF.mat

trim.method='sd'; trim.centralPct=95; trim.edgesTHW=0:0.25:10; trim.minPts=5;

all_work_L=[]; all_pid_L=[];
all_work_H=[]; all_pid_H=[];
excluded=[];

for k=1:numel(participants_CF)
    if ismember(k,excluded), continue, end
    THW=participants_CF(k).THW(:);
    thisWork=participants_CF(k).pupil_diameter_scaled(:);

    Work_n=local_normalize01(thisWork);
    keep=local_buildKeepMask(THW,Work_n,trim);
    Work_keep=Work_n(keep);

    if strcmpi(participants_CF(k).drivingExperience,'L')
        all_work_L=[all_work_L;Work_keep];
        all_pid_L=[all_pid_L;repmat(k,numel(Work_keep),1)];
    elseif strcmpi(participants_CF(k).drivingExperience,'H')
        all_work_H=[all_work_H;Work_keep];
        all_pid_H=[all_pid_H;repmat(k,numel(Work_keep),1)];
    end
end

meanPD_L=splitapply(@mean,all_work_L,findgroups(all_pid_L));
meanPD_H=splitapply(@mean,all_work_H,findgroups(all_pid_H));

figure
boxplot([meanPD_L; meanPD_H], ...
        [repelem({'Low'},numel(meanPD_L))' ; ...
         repelem({'High'},numel(meanPD_H))' ], ...
        'Symbol','')
ylabel('Pupil diameter')
title('Pupil diameter – Low vs High Driving Experience')
set(gca,'XLim',[0.5 2.5])

sd_L=std(meanPD_L); sd_H=std(meanPD_H);
[H,p,CI,stats]=ttest2(meanPD_L,meanPD_H);
fprintf('\nLow  mean = %.3f (SD = %.4f)  |  High mean = %.3f (SD = %.4f)\n', ...
        mean(meanPD_L),sd_L,mean(meanPD_H),sd_H)
fprintf('t(%d) = %.2f,  p = %.4f\n\n',stats.df,stats.tstat,p)

if ~exist('./Figure','dir'), mkdir('./Figure'); end
exportgraphics(gcf,'./Figure/Pupil_Experience_Low_vs_High_trimmed.png','Resolution',300)

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
