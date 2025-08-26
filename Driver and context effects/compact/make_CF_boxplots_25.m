function make_CF_boxplots_25_compact_v2
%% 25 compact boxplots: 5 proxies × 5 groupings
clc; close all

outdir = './Figure/FeatureEffect/Compact';
if ~exist(outdir,'dir'), mkdir(outdir); end
load participants_CF.mat   % struct array participants_CF(1×N)

% ----- proxies: label + short + candidate field names（大小写兼容） -----
proxies = {
    struct('label','Effort','short','Effort', ...
           'cands',{{'effort_scaled','Effort_Scaled','effort','Effort'}})
    struct('label','EEG-workload','short','EEG', ...
           'cands',{{'eeg_workload_scaled','wl_eeg_scaled','WL_EEG_Scaled','EEG_Workload_Scaled'}})
    struct('label','Pupil','short','Pupil', ...
           'cands',{{'pupil_scaled','pupil_diameter_scaled','Pupil_Scaled','Pupil_Diameter_Scaled'}})
    struct('label','Eyelid','short','Eyelid', ...
           'cands',{{'eyelid_scaled','Eyelid_Scaled'}})
    struct('label','Gaze variability','short','GazeVar', ...
           'cands',{{'gaze_variability_scaled','GazeVariability_Scaled'}})
    };

% ----- groupings -----
groups = ["Age","Gender","Experience","DayNight","NDRT"];

% ----- trimming settings（与你原脚本一致） -----
trim.method      = 'sd';
trim.centralPct  = 95;
trim.edgesTHW    = 0 : 0.25 : 10;
trim.minPts      = 5;

% ----- figure 样式 -----
figW = 6.0; figH = 4.2;      % cm，更宽以容纳水平 x 标签
pos  = [1, 1.10];            % 两盒中心更近
wid  = 0.10;                 % 盒更窄

for gi = 1:numel(groups)
    G = groups(gi);
    for pi = 1:numel(proxies)
        Pdef = proxies{pi};
        proxy_field = resolve_field(participants_CF, Pdef.cands);

        % 数据：每参与者均值（先按 THW 分箱裁剪、再 [0,1] 归一）
        [Y1,Y2,g1,g2] = collect_means(participants_CF, proxy_field, string(G), trim);

        % 自适应 y 轴范围（含 15% 边距，裁到 [0,1] 内）
        yall = [Y1; Y2];
        rng  = max(yall) - min(yall);
        pad  = 0.15 * (rng + eps);
        ylo  = max(0, min(yall) - pad);
        yhi  = min(1, max(yall) + pad);

        % 统计打印
        [~,p,~,stats] = ttest2(Y1,Y2);
        fprintf('%-14s | %-10s: %s = %.3f±%.3f  vs  %s = %.3f±%.3f  |  t(%d)=%.2f, p=%.4f\n',...
            Pdef.label, G, g1, mean(Y1), std(Y1), g2, mean(Y2), std(Y2), stats.df, stats.tstat, p);

        % 绘图（无标题，x 标签水平）
% ------- 紧凑 boxplot（更窄 + 更近；左右边距极小；x 标签水平） -------
fig = figure('Visible','off','Units','centimeters','Position',[0 0 6.0 4.2]);

labels = [repmat({g1}, numel(Y1),1); repmat({g2}, numel(Y2),1)];

% 盒子中心的间距（越小越“靠近”），以及盒宽（越小越“变细”）
posGap = 0.05;                 % ← 把 0.35 改成 0.16（建议范围 0.12–0.20）
pos    = [1, 1 + posGap];      % 两个盒子的中心横坐标
wid    = 0.025;                % ← 变细（可再试 0.04）

boxplot([Y1;Y2], labels, ...
        'Positions', pos, ...
        'Widths',    wid, ...
        'Symbol','', 'Notch','off');

ax = gca;
ax.XTick      = pos;
ax.XTickLabel = {g1,g2};
ax.TickLabelInterpreter = 'none';
try, xtickangle(ax,0); end

% 左右留白：跟随盒宽给一个很小的 padding，这样不会离 y 轴/右边太远
xpad   = max(0.06, 0.9*wid);     % ← 进一步减小空白，可按需调 0.05–0.08
ax.XLim = [pos(1)-xpad, pos(2)+xpad];

% y 轴：根据数据自适应（含 15% 边距，并裁在 [0,1] 内）
yall = [Y1; Y2];
yrng = max(yall) - min(yall);
ypad = 0.15 * (yrng + eps);
ax.YLim  = [max(0, min(yall)-ypad), min(1, max(yall)+ypad)];
ax.YTick = nice_ticks(ax.YLim, 5);

% 其余样式
ax.TickDir   = 'out';
ax.LineWidth = 0.9;
ax.FontSize  = 9;
ylabel(Pdef.label,'FontSize',11,'FontWeight','bold');
box off
set(ax,'LooseInset', max(get(ax,'TightInset'), [0.02 0.26 0.02 0.02]));

% 导出
fname = sprintf('%s_%s_compact.png', char(G), Pdef.short);
exportgraphics(fig, fullfile(outdir,fname), 'Resolution',300,'ContentType','vector');
close(fig)

    end
end

fprintf('\nDone. Files saved in: %s\n', outdir);
end

%% ======================= helpers =======================
function fieldname = resolve_field(S, candidates)
    fn = fieldnames(S);
    for i = 1:numel(candidates)
        idx = find(strcmpi(fn, candidates{i}), 1);
        if ~isempty(idx), fieldname = fn{idx}; return, end
    end
    error('None of the candidate fields found: %s', strjoin(candidates, ', '));
end

function [Y1_means, Y2_means, g1_name, g2_name] = collect_means(P, proxy_field, group_name, trim)
    all_Y1 = []; all_id1 = [];
    all_Y2 = []; all_id2 = [];

    for k = 1:numel(P)
        THW = P(k).THW(:);
        Y   = P(k).(proxy_field)(:);

        Y_n   = normalize01(Y);                % 每人 [0,1]
        keep  = buildKeepMask(THW, Y_n, trim); % THW 分箱 95% 中央带
        Y_use = Y_n(keep);

        switch lower(group_name)
            case 'age'
                condA   = strcmpi(P(k).ageGroup,'A');
                g1_name = 'Age A';
                g2_name = 'Age B';
            case 'gender'
                condA   = strcmpi(P(k).genderGroup,'M');
                g1_name = 'Male';   g2_name = 'Female';
            case 'experience'
                condA   = strcmpi(P(k).drivingExperience,'L');
                g1_name = 'Low';    g2_name = 'High';
            case 'daynight'
                fname   = get_if_field(P(k),'filename','');
                condA   = contains(fname,'i4Driving_D'); % Day
                g1_name = 'Day';    g2_name = 'Night';
            case 'ndrt'
                fname   = get_if_field(P(k),'filename','');
                condA   = (numel(fname) >= 13 && fname(13) == '1');  % NDRT
                g1_name = 'NDRT';   g2_name = 'No NDRT';
            otherwise
                error('Unknown group_name: %s', group_name);
        end

        if condA
            all_Y1 = [all_Y1;  Y_use];
            all_id1 = [all_id1; repmat(k,numel(Y_use),1)];
        else
            all_Y2 = [all_Y2;  Y_use];
            all_id2 = [all_id2; repmat(k,numel(Y_use),1)];
        end
    end

    Y1_means = splitapply(@mean, all_Y1, findgroups(all_id1));
    Y2_means = splitapply(@mean, all_Y2, findgroups(all_id2));
end

function keep = buildKeepMask(x,y,t)
    keep = false(size(x));
    if strcmpi(t.method,'sd'), k = sqrt(2)*erfinv(t.centralPct/100); end
    for b = 1:numel(t.edgesTHW)-1
        in = (x>=t.edgesTHW(b)) & (x<t.edgesTHW(b+1));
        if nnz(in) < t.minPts, keep(in)=true; continue, end
        switch lower(t.method)
            case 'percentile'
                tail=(100-t.centralPct)/2;
                lo=prctile(y(in),tail); hi=prctile(y(in),100-tail);
                keep(in)= (y(in)>=lo & y(in)<=hi);
            case 'sd'
                mu=mean(y(in)); sd=std(y(in),0);
                keep(in)= abs(y(in)-mu) <= k*sd;
        end
    end
end

function y_n = normalize01(y)
    d = max(y)-min(y);
    if d>0, y_n=(y-min(y))./d; else, y_n=0.5*ones(size(y)); end
end

function v = get_if_field(S, fname, default)
    if isfield(S, fname), v = S.(fname); else, v = default; end
end

function ticks = nice_ticks(ylim_, maxN)
% 生成不超过 maxN 个等距刻度
    lo = ylim_(1); hi = ylim_(2); 
    if hi<=lo, ticks = lo; return, end
    span = hi - lo;
    step = 10^floor(log10(span/maxN));
    for m = [1,2,5,10]
        cand = step*m;
        if ceil(span/cand) <= maxN
            step = cand; break
        end
    end
    first = ceil(lo/step)*step;
    last  = floor(hi/step)*step;
    ticks = first:step:last;
    if isempty(ticks) || numel(ticks)==1
        ticks = [lo, hi];
    end
end
