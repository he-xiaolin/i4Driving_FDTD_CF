%% EEG workload – Age Group A vs B (rescaled + trimmed)
clc; clear; close all
load participants_CF.mat

trim.method='sd'; trim.centralPct=95; trim.edgesTHW=0:0.25:10; trim.minPts=5;

all_work_A=[]; all_pid_A=[];
all_work_B=[]; all_pid_B=[];
excluded=[];

for k=1:numel(participants_CF)
    if ismember(k,excluded), continue, end
    THW=participants_CF(k).THW(:);
    thisWork=participants_CF(k).eeg_workload_scaled(:);

    Work_n=local_normalize01(thisWork);
    keep=local_buildKeepMask(THW,Work_n,trim);
    Work_keep=Work_n(keep);

    if strcmpi(participants_CF(k).ageGroup,'A')      % (<31)
        all_work_A=[all_work_A;Work_keep];
        all_pid_A=[all_pid_A;repmat(k,numel(Work_keep),1)];
    elseif strcmpi(participants_CF(k).ageGroup,'B')  % (>=31)
        all_work_B=[all_work_B;Work_keep];
        all_pid_B=[all_pid_B;repmat(k,numel(Work_keep),1)];
    end
end

meanHR_A=splitapply(@mean,all_work_A,findgroups(all_pid_A));
meanHR_B=splitapply(@mean,all_work_B,findgroups(all_pid_B));

figure
boxplot([meanHR_A; meanHR_B], ...
        [repelem({'Age A (<31)'},numel(meanHR_A))' ; ...
         repelem({'Age B (>=31)'},numel(meanHR_B))' ], ...
        'Symbol','')

ylabel('EEG workload')
title('EEG workload – Age Group A vs B')
set(gca,'XLim',[0.5 2.5])

sd_A=std(meanHR_A); sd_B=std(meanHR_B);
[H,p,CI,stats]=ttest2(meanHR_A,meanHR_B);
fprintf('\nAge A mean = %.3f (SD = %.4f)  |  Age B mean = %.3f (SD = %.4f)\n', ...
        mean(meanHR_A),sd_A,mean(meanHR_B),sd_B)
fprintf('t(%d) = %.2f,  p = %.4f\n\n',stats.df,stats.tstat,p)

if ~exist('./Figure','dir'), mkdir('./Figure'); end
exportgraphics(gcf,'./Figure/HeartRate_Age_A_vs_B_trimmed.png','Resolution',300)

function keep=local_buildKeepMask(x,y,t)
    keep=false(size(x));
    if strcmpi(t.method,'sd'), k=sqrt(2)*erfinv(t.centralPct/100); end
    for b=1:numel(t.edgesTHW)-1
        in=x>=t.edgesTHW(b) & x<t.edgesTHW(b+1);
        if nnz(in)<t.minPts, keep(in)=true; continue, end
        switch lower(t.method)
            case 'percentile'
                tail=(100-t.centralPct)/2; lo=prctile(y(in),tail); hi=prctile(y(in),100-tail);
                keep(in)=y(in)>=lo & y(in)<=hi;
            case 'sd'
                mu=mean(y(in)); sd=std(y(in),0);
                keep(in)=abs(y(in)-mu)<=k*sd;
        end
    end
end
function y_n=local_normalize01(y)
    d=max(y)-min(y);
    if d>0, y_n=(y-min(y))./d; else, y_n=0.5*ones(size(y)); end
end
