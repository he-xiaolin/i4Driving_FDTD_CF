%% Gaze variability  – Gender  (rescaled (within participant) + trimmed)

clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ───────── trimming settings (same as Effort scripts) ─────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central band
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ───────── collectors ─────────
all_Y_M = [];  all_pid_M = [];
all_Y_F = [];  all_pid_F = [];


        excluded = [];

        for k = 1:numel(participants_CF)
            if ismember(k,excluded),  continue,  end

            THW = participants_CF(k).THW(:);
            Y   = participants_CF(k).gaze_variability_scaled(:);

            Y_n = local_normalize01(Y);
            keep = local_buildKeepMask(THW,Y_n,trim);
            Y_keep = Y_n(keep);

            if strcmpi(participants_CF(k).genderGroup,'M')      % Male
                all_Y_M = [all_Y_M ; Y_keep];
                all_pid_M = [all_pid_M ; repmat(k,numel(Y_keep),1)];
            elseif strcmpi(participants_CF(k).genderGroup,'F')  % Female
                all_Y_F = [all_Y_F ; Y_keep];
                all_pid_F = [all_pid_F ; repmat(k,numel(Y_keep),1)];
            end
        end

        %% Participant-level means
        meanY_M = splitapply(@mean, all_Y_M, findgroups(all_pid_M));
        meanY_F = splitapply(@mean, all_Y_F, findgroups(all_pid_F));

        %% Box-plot + statistics
        figure
        labels = [repmat({'Male'},   numel(meanY_M), 1); ...
                  repmat({'Female'}, numel(meanY_F), 1)];
        boxplot([meanY_M; meanY_F], labels, 'Symbol','')
        ylabel('Gaze variability')
        title('Gaze variability  – Male vs Female')
        set(gca,'XLim',[0.5 2.5])

        [H,p,CI,stats] = ttest2(meanY_M,meanY_F);
        std_M = std(meanY_M);  std_F = std(meanY_F);

        fprintf('\nMale mean = %.3f (SD = %.4f)   |  Female mean = %.3f (SD = %.4f)\n', ...
                mean(meanY_M), std_M, mean(meanY_F), std_F)
        fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

        if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
        exportgraphics(gcf,'./Figure/GazeVariability_Gender_Male_vs_Female_trimmed.png','Resolution',300)
        

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
