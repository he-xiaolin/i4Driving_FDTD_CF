%% EEG workload – Male vs Female (rescaled + trimmed)
clc; clear; close all
load participants_CF.mat

trim.method='sd'; trim.centralPct=95; trim.edgesTHW=0:0.25:10; trim.minPts=5;

all_work_M=[]; all_pid_M=[];
all_work_F=[]; all_pid_F=[];
excluded=[];

for k=1:numel(participants_CF)
    if ismember(k,excluded), continue, end
    THW=participants_CF(k).THW(:);
    thisWork=participants_CF(k).eeg_workload_scaled(:);

    Work_n=local_normalize01(thisWork);
    keep=local_buildKeepMask(THW,Work_n,trim);
    Work_keep=Work_n(keep);

    if strcmpi(participants_CF(k).genderGroup,'M')
        all_work_M=[all_work_M;Work_keep];
        all_pid_M=[all_pid_M;repmat(k,numel(Work_keep),1)];
    elseif strcmpi(participants_CF(k).genderGroup,'F')
        all_work_F=[all_work_F;Work_keep];
        all_pid_F=[all_pid_F;repmat(k,numel(Work_keep),1)];
    end
end

meanHR_M=splitapply(@mean,all_work_M,findgroups(all_pid_M));
meanHR_F=splitapply(@mean,all_work_F,findgroups(all_pid_F));

figure
boxplot([meanHR_M; meanHR_F], ...
        [repelem({'Male'},numel(meanHR_M))' ; ...
         repelem({'Female'},numel(meanHR_F))' ], ...
        'Symbol','')

ylabel('EEG workload')
title('EEG workload – Male vs Female')
set(gca,'XLim',[0.5 2.5])

sd_M=std(meanHR_M); sd_F=std(meanHR_F);
[H,p,CI,stats]=ttest2(meanHR_M,meanHR_F);
fprintf('\nMale   mean = %.3f (SD = %.4f)  |  Female mean = %.3f (SD = %.4f)\n', ...
        mean(meanHR_M),sd_M,mean(meanHR_F),sd_F)
fprintf('t(%d) = %.2f,  p = %.4f\n\n',stats.df,stats.tstat,p)

if ~exist('./Figure','dir'), mkdir('./Figure'); end
exportgraphics(gcf,'./Figure/HeartRate_Gender_Male_vs_Female_trimmed.png','Resolution',300)

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
