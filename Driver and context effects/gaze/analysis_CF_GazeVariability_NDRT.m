%% Gaze variability  – NDRT  (rescaled (within participant) + trimmed)

clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ───────── trimming settings (same as Effort scripts) ─────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central band
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ───────── collectors ─────────
all_Y_NDRT   = [];   all_pid_NDRT   = [];
all_Y_noNDRT = [];   all_pid_noNDRT = [];


        excluded = [];

        isNDRTfile = @(fname) (numel(fname) >= 13 && fname(13) == '1');

        for k = 1:numel(participants_CF)
            if ismember(k,excluded),  continue,  end

            THW = participants_CF(k).THW(:);
            Y   = participants_CF(k).gaze_variability_scaled(:);

            Y_n = local_normalize01(Y);
            keep = local_buildKeepMask(THW,Y_n,trim);
            Y_keep = Y_n(keep);

            if isNDRTfile(participants_CF(k).filename)
                all_Y_NDRT   = [all_Y_NDRT ;   Y_keep];
                all_pid_NDRT = [all_pid_NDRT ; repmat(k,numel(Y_keep),1)];
            else
                all_Y_noNDRT   = [all_Y_noNDRT ;   Y_keep];
                all_pid_noNDRT = [all_pid_noNDRT ; repmat(k,numel(Y_keep),1)];
            end
        end

        %% Participant-level means
        meanY_NDRT   = splitapply(@mean, all_Y_NDRT  , findgroups(all_pid_NDRT  ));
        meanY_noNDRT = splitapply(@mean, all_Y_noNDRT, findgroups(all_pid_noNDRT));

        %% Box-plot + statistics
        figure
        labels = [repmat({'NDRT'},   numel(meanY_NDRT), 1); ...
                  repmat({'No NDRT'}, numel(meanY_noNDRT), 1)];
        boxplot([meanY_NDRT; meanY_noNDRT], labels, 'Symbol','')
        ylabel('Gaze variability')
        title('Gaze variability  – NDRT vs NoNDRT')
        set(gca,'XLim',[0.5 2.5])

        [H,p,CI,stats] = ttest2(meanY_NDRT,meanY_noNDRT);
        std_NDRT   = std(meanY_NDRT);
        std_noNDRT = std(meanY_noNDRT);

        fprintf('\nNDRT mean = %.3f (SD = %.4f)   |  NoNDRT mean = %.3f (SD = %.4f)\n', ...
                mean(meanY_NDRT), std_NDRT, mean(meanY_noNDRT), std_noNDRT)
        fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

        if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
        exportgraphics(gcf,'./Figure/GazeVariability_NDRT_vs_NoNDRT_trimmed.png','Resolution',300)
        

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
