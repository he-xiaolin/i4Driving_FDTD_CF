clc; clear; close all; 
load('participants_CF.mat');   % 需要字段: participants_CF(k).TTC, .effort_scaled

%% ---------------- 参数区 ----------------
numP  = numel(participants_CF);
edges = linspace(0,10,121);                 % TTC 分箱（按需改）
binCtr = edges(1:end-1) + diff(edges)/2;

normalizeEffort = true;                     % 是否对每个被试的 effort(TD) 做 min-max
bandMode   = 'sd';                          % 'sd' | 'sem' | 'iqr' | 'p90'
smoothWin  = 1;                             % 跨被试聚合后的移动平均窗口 (bin)；1=不平滑
minProfilesPerBin = 5;                      % 至少多少个被试在该 bin 有数据才绘制
exclude = [22 52 57 64];                    % 合并曲线时排除的人(可设为空 [])

% 统一的图像输出文件夹
figDir = fullfile(pwd,'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

%% -------- 容器（Figure 3 用；以及跨被试聚合用） --------
r_TTC_TD = nan(numP,1);
p_TTC_TD = nan(numP,1);
n_pairs  = zeros(numP,1);

B = numel(edges)-1;
muPerP = nan(numP, B);                      % 每个被试在每个 bin 的均值 (min-max 后)

%% -------- Figure 1（隐藏）：8×8 单人 ribbon，min-max 后 --------
f1 = figure('Name','Per-profile TTC vs TD (8x8)','Units','normalized', ...
            'OuterPosition',[0 0 1 1],'Visible','off');

for k = 1:numP
    X = participants_CF(k).TTC(:);
    TD = participants_CF(k).effort_scaled(:);

    % 相关系数（Pearson；与是否 min-max 无关）
    v = isfinite(X) & isfinite(TD);
    n_pairs(k) = sum(v);
    if n_pairs(k) >= 2
        [R,P] = corrcoef(X(v), TD(v));
        r_TTC_TD(k) = R(1,2);
        p_TTC_TD(k) = P(1,2);
    end

    % 每个被试对 TD 做 min-max（仅用于绘图与跨被试聚合）
    if normalizeEffort
        eMin = min(TD,[],'omitnan'); eMax = max(TD,[],'omitnan');
        if ~isnan(eMin) && ~isnan(eMax) && eMax > eMin
            TD = (TD - eMin) / (eMax - eMin);
        else
            TD(:) = NaN;
        end
    end

    % 分箱统计（个人）
    [Xs, idx] = sort(X);
    TDs = TD(idx);
    mu = nan(1,B); sd = nan(1,B);
    for i = 1:B
        m = (Xs >= edges(i)) & (Xs < edges(i+1));
        if any(m)
            mu(i) = mean(TDs(m), 'omitnan');
            sd(i) = std( TDs(m), 'omitnan');
        end
    end
    muPerP(k,:) = mu;

    % 子图（个人 ribbon）
    ax = subplot(8,8,k,'Parent',f1); hold(ax,'on');
    keep = ~isnan(mu);
    x = binCtr(keep); muPlot = mu(keep); sdPlot = sd(keep);
    if ~isempty(x)
        fill(ax,[x, fliplr(x)], [muPlot - sdPlot, fliplr(muPlot + sdPlot)], ...
             'b','FaceAlpha',0.18,'EdgeColor','none');
        plot(ax,x,muPlot,'b-','LineWidth',1.1);
    end
    title(ax,sprintf('P%02d',k));
    xlim(ax,[0 10]); ylim(ax,[0 1]); grid(ax,'on');
end

exportgraphics(f1, fullfile(figDir,'Fig1_PerProfile_TTC_vs_TD_8x8.png'), 'Resolution',300);
close(f1);  % 不在屏幕上显示

%% -------- Figure 2：合并曲线（先“个人分箱”，再“跨人聚合”） --------
useRows = true(numP,1);
useRows(exclude) = false;                   % 若无需排除，设 exclude = []
muP = muPerP(useRows,:);

meanAcross = nan(1,B);
lo = nan(1,B); hi = nan(1,B);

for i = 1:B
    col = muP(:,i);
    col = col(~isnan(col));
    nProf = numel(col);
    if nProf >= minProfilesPerBin
        m = mean(col);
        meanAcross(i) = m;
        switch lower(bandMode)
            case 'sd'
                s = std(col);
                lo(i) = m - s; hi(i) = m + s;
            case 'sem'
                s = std(col) / sqrt(nProf);
                lo(i) = m - s; hi(i) = m + s;
            case 'iqr'
                q = prctile(col,[25 75]);
                lo(i) = q(1);  hi(i) = q(2);
            case 'p90'
                q = prctile(col,[5 95]);
                lo(i) = q(1);  hi(i) = q(2);
            otherwise
                error('bandMode must be sd/sem/iqr/p90');
        end
    end
end

% 可选平滑（这里默认 smoothWin=1，不执行）
if smoothWin > 1
    meanAcross = movmean(meanAcross, smoothWin, 'omitnan');
    lo         = movmean(lo,         smoothWin, 'omitnan');
    hi         = movmean(hi,         smoothWin, 'omitnan');
end

keep = ~isnan(meanAcross) & ~isnan(lo) & ~isnan(hi);
x = binCtr(keep); mu = meanAcross(keep); L = lo(keep); H = hi(keep);

figure(2); clf;
plot(x, mu, 'b-', 'LineWidth', 2); hold on;
fill([x, fliplr(x)], [L, fliplr(H)], 'b', ...
     'FaceAlpha', 0.18, 'EdgeColor', 'none');
xlabel('TTC (s)'); ylabel('Effort');
title('TTC vs Effort');
xlim([0 10]); ylim([0 1]); 
% grid on;

exportgraphics(gcf, fullfile(figDir, sprintf('Fig2_TTC_vs_TD_Combined_%s.png', bandMode)), ...
               'Resolution',300);

%% -------- Figure 3：每人相关系数箱线图 --------
validR = ~isnan(r_TTC_TD);
figure(3); clf;
boxplot(r_TTC_TD(validR),'Symbol','', 'Labels', {''}); hold on;
yline(0,'r--','LineWidth',1);
ylabel('Correlation (r)'); title('Correlation of TTC vs Effort per Profile');
% grid on;

exportgraphics(gcf, fullfile(figDir,'Fig3_Boxplot_r_TTC_TD.png'), 'Resolution',300);

%% -------- 相关系数明细另存 --------
T = table((1:numP).', r_TTC_TD, p_TTC_TD, n_pairs, ...
          'VariableNames',{'Profile','r','p','N_pairs'});
writetable(T,'TTC_TD_correlation_per_profile.csv');
save('TTC_TD_correlation_per_profile.mat','T');
