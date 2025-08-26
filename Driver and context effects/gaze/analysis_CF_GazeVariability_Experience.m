%% Gaze variability  – Experience  (rescaled (within participant) + trimmed)

clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ───────── trimming settings (same as Effort scripts) ─────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central band
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ───────── collectors ─────────
all_Y_L = [];  all_pid_L = [];
all_Y_H = [];  all_pid_H = [];


        excluded = [];

        for k = 1:numel(participants_CF)
            if ismember(k,excluded),  continue,  end

            THW = participants_CF(k).THW(:);
            Y   = participants_CF(k).gaze_variability_scaled(:);

            Y_n = local_normalize01(Y);
            keep = local_buildKeepMask(THW,Y_n,trim);
            Y_keep = Y_n(keep);

            if strcmpi(participants_CF(k).drivingExperience,'L')      % Low
                all_Y_L = [all_Y_L ; Y_keep];
                all_pid_L = [all_pid_L ; repmat(k,numel(Y_keep),1)];
            elseif strcmpi(participants_CF(k).drivingExperience,'H')  % High
                all_Y_H = [all_Y_H ; Y_keep];
                all_pid_H = [all_pid_H ; repmat(k,numel(Y_keep),1)];
            end
        end

        %% Participant-level means
        meanY_L = splitapply(@mean, all_Y_L, findgroups(all_pid_L));
        meanY_H = splitapply(@mean, all_Y_H, findgroups(all_pid_H));

        %% Box-plot + statistics
        figure
        labels = [repmat({'Low'},  numel(meanY_L), 1); ...
                  repmat({'High'}, numel(meanY_H), 1)];
        boxplot([meanY_L; meanY_H], labels, 'Symbol','')
        ylabel('Gaze variability')
        title('Gaze variability  – Low vs High Driving Experience')
        set(gca,'XLim',[0.5 2.5])

        [H,p,CI,stats] = ttest2(meanY_L,meanY_H);
        std_L = std(meanY_L);  std_H = std(meanY_H);

        fprintf('\nLow  mean = %.3f (SD = %.4f)   |  High mean = %.3f (SD = %.4f)\n', ...
                mean(meanY_L), std_L, mean(meanY_H), std_H)
        fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

        if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
        exportgraphics(gcf,'./Figure/GazeVariability_Experience_Low_vs_High_trimmed.png','Resolution',300)
        

% ───────────────────── helper functions ─────────────────────
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
    if d>0
        y_n=(y-min(y))./d;      % map to [0,1]
    else
        y_n=0.5*ones(size(y));  % constant vector→mid-level
    end
end
