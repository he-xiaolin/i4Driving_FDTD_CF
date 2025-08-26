%% Eyelid  – DayNight  (rescaled (within participant) + trimmed)

clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ───────── trimming settings (same as Effort scripts) ─────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central band
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ───────── collectors ─────────
all_Y_Day   = [];   all_pid_Day   = [];
all_Y_Night = [];   all_pid_Night = [];


        excluded = [];

        for k = 1:numel(participants_CF)
            if ismember(k,excluded),  continue,  end

            THW = participants_CF(k).THW(:);
            Y   = participants_CF(k).eyelid_scaled(:);

            Y_n = local_normalize01(Y);
            keep = local_buildKeepMask(THW,Y_n,trim);
            Y_keep = Y_n(keep);

            if contains(participants_CF(k).filename,'i4Driving_D')      % Day
                all_Y_Day   = [all_Y_Day ;   Y_keep];
                all_pid_Day = [all_pid_Day ; repmat(k,numel(Y_keep),1)];
            elseif contains(participants_CF(k).filename,'i4Driving_N')  % Night
                all_Y_Night   = [all_Y_Night ;   Y_keep];
                all_pid_Night = [all_pid_Night ; repmat(k,numel(Y_keep),1)];
            end
        end

        %% Participant-level means
        meanY_Day   = splitapply(@mean, all_Y_Day  , findgroups(all_pid_Day  ));
        meanY_Night = splitapply(@mean, all_Y_Night, findgroups(all_pid_Night));

        %% Box-plot + statistics
        figure
        labels = [repmat({'Day'},   numel(meanY_Day), 1); ...
                  repmat({'Night'}, numel(meanY_Night), 1)];
        boxplot([meanY_Day; meanY_Night], labels, 'Symbol','')
        ylabel('Eyelid')
        title('Eyelid  – Day vs Night')
        set(gca,'XLim',[0.5 2.5])

        [H,p,CI,stats] = ttest2(meanY_Day,meanY_Night);
        std_Day   = std(meanY_Day);
        std_Night = std(meanY_Night);

        fprintf('\nDay  mean = %.3f (SD = %.4f)   |  Night mean = %.3f (SD = %.4f)\n', ...
                mean(meanY_Day), std_Day, mean(meanY_Night), std_Night)
        fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

        if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
        exportgraphics(gcf,'./Figure/Eyelid_Day_vs_Night_trimmed.png','Resolution',300)
        

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
