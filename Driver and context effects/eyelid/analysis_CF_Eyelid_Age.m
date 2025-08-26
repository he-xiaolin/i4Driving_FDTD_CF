%% Eyelid  – Age  (rescaled (within participant) + trimmed)

clc; clear; close all
load participants_CF.mat          % → struct array (1×N)

% ───────── trimming settings (same as Effort scripts) ─────────
trim.method      = 'sd';          % 'sd'  or  'percentile'
trim.centralPct  = 95;            % keep central band
trim.edgesTHW    = 0 : 0.25 : 10; % THW bins for local trimming
trim.minPts      = 5;             % leave tiny bins untrimmed

% ───────── collectors ─────────
all_Y_A = [];  all_pid_A = [];
all_Y_B = [];  all_pid_B = [];


        excluded = [];

        for k = 1:numel(participants_CF)
            if ismember(k,excluded),  continue,  end

            THW = participants_CF(k).THW(:);
            Y   = participants_CF(k).eyelid_scaled(:);

            Y_n = local_normalize01(Y);
            keep = local_buildKeepMask(THW,Y_n,trim);
            Y_keep = Y_n(keep);

            if strcmpi(participants_CF(k).ageGroup,'A')      % Age A (<31)
                all_Y_A = [all_Y_A ; Y_keep];
                all_pid_A = [all_pid_A ; repmat(k,numel(Y_keep),1)];
            elseif strcmpi(participants_CF(k).ageGroup,'B')  % Age B (>=31)
                all_Y_B = [all_Y_B ; Y_keep];
                all_pid_B = [all_pid_B ; repmat(k,numel(Y_keep),1)];
            end
        end

        %% Participant-level means
        meanY_A = splitapply(@mean, all_Y_A, findgroups(all_pid_A));
        meanY_B = splitapply(@mean, all_Y_B, findgroups(all_pid_B));

        %% Box-plot + statistics
        figure
        labels = [repmat({'Age A (<31)'}, numel(meanY_A), 1); repmat({'Age B (>=31)'}, numel(meanY_B), 1)];
        boxplot([meanY_A; meanY_B], labels, 'Symbol','')
        ylabel('Eyelid')
        title('Eyelid  – Age Group A vs B')
        set(gca,'XLim',[0.5 2.5])

        [H,p,CI,stats] = ttest2(meanY_A,meanY_B);
        std_A = std(meanY_A);  std_B = std(meanY_B);

        fprintf('\nAge A mean = %.3f (SD = %.4f)   |  Age B mean = %.3f (SD = %.4f)\n', ...
                mean(meanY_A), std_A, mean(meanY_B), std_B)
        fprintf('t(%d) = %.2f,  p = %.4f\n\n', stats.df, stats.tstat, p)

        if ~exist('./Figure','dir'),  mkdir('./Figure'),  end
        exportgraphics(gcf,'./Figure/Eyelid_Age_A_vs_B_trimmed.png','Resolution',300)
        

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
