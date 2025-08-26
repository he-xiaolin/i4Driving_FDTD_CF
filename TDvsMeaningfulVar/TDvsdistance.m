clc; clear; close all; 
load('participants_CF.mat');   % 需要字段: participants_CF(k).distance_leader, .effort_scaled

%% ---------------- 参数区 ----------------
numP  = numel(participants_CF);
edges = linspace(0,10,121);                 % 距离分箱（按需改，例如 0~60）
binCtr = edges(1:end-1) + diff(edges)/2;

normalizeEffort = true;                     % 是否对每个被试的 effort(TD) 做 min-max
bandMode   = 'sd';                          % 'sd' | 'sem' | 'iqr' | 'p90'
smoothWin  = 1;                             % 1=不平滑
minProfilesPerBin = 5;                      % 至少多少个被试在该 bin 有数据才绘制
exclude = [22 52 57 64];                    % 合并曲线时排除的人(可设为空 [])

% 统一的图像输出文件夹
figDir = fullfile(pwd,'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

%% -------- 容器 --------
r_DIST_TD = nan(numP,1);
p_DIST_TD = nan(numP,1);
n_pairs   = zeros(numP,1);

B = numel(edges)-1;
muPerP = nan(numP, B);

%% -------- Figure 1（隐藏）：8×8 单人 ribbon，min-max 后 --------
f1 = figure('Name','Per-profile Distance vs TD (8x8)','Units','normalized', ...
            'OuterPosition',[0 0 1 1],'Visible','off');

for k = 1:numP
    X  = participants_CF(k).distance_leader(:);
    TD = participants_CF(k).effort_scaled(:);

    % 相关系数
    v = isfinite(X) & isfinite(TD);
    n_pairs(k) = sum(v);
    if n_pairs(k) >= 2
        [R,P] = corrcoef(X(v), TD(v));
        r_DIST_TD(k) = R(1,2);
        p_DIST_TD(k) = P(1,2);
    end

    % 每个被试对 TD 做 min-max
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
    xlim(ax,[0 10]); ylim(ax,[0 1]); grid(ax,'on');   % 距离范围不同可改
end

exportgraphics(f1, fullfile(figDir,'Fig1_PerProfile_Distance_vs_TD_8x8.png'), 'Resolution',300);
close(f1);

%% -------- Figure 2: 合并曲线（按样本池聚合） --------
% 可选：把每人TD先min-max再合并
allX  = [];
allTD = [];
for k = 1:numP
    if ismember(k, exclude), continue; end
    X  = participants_CF(k).distance_leader(:);   % 换成 TTC 则用 .TTC
    TD = participants_CF(k).effort_scaled(:);

    % min-max（仅用于绘图；相关分析仍用原始TD）
    if normalizeEffort
        eMin = min(TD,[],'omitnan'); eMax = max(TD,[],'omitnan');
        if ~isnan(eMin) && ~isnan(eMax) && eMax > eMin
            TD = (TD - eMin) / (eMax - eMin);
        else
            TD(:) = NaN;
        end
    end

    allX  = [allX;  X];
    allTD = [allTD; TD];
end

% 设置距离分箱（确保 0~100m）
edges  = linspace(0,100,121);
binCtr = edges(1:end-1) + diff(edges)/2;

% 直接在“样本池”上分箱
[allXs, idx] = sort(allX);
allTDs = allTD(idx);

mu = nan(size(binCtr));
sd = nan(size(binCtr));

for i = 1:numel(edges)-1
    m = (allXs >= edges(i)) & (allXs < edges(i+1));
    if any(m)
        mu(i) = mean(allTDs(m), 'omitnan');
        switch lower(bandMode)
            case 'sd'
                s  = std(allTDs(m), 'omitnan');
                lo = mu(i) - s;  hi = mu(i) + s;
            case 'sem'
                s  = std(allTDs(m), 'omitnan') / sqrt(sum(m));
                lo = mu(i) - s;  hi = mu(i) + s;
            case 'iqr'
                q  = prctile(allTDs(m), [25 75]);
                lo = q(1); hi = q(2);
            case 'p90'
                q  = prctile(allTDs(m), [5 95]);
                lo = q(1); hi = q(2);
        end
        L(i) = lo; H(i) = hi;   %#ok<AGROW>
    end
end

keep = ~isnan(mu);
x = binCtr(keep); mu = mu(keep); 
L = L(keep); H = H(keep);

figure(2); clf;
plot(x, mu, 'b-', 'LineWidth', 2); hold on;
fill([x, fliplr(x)], [L, fliplr(H)], 'b', 'FaceAlpha',0.18, 'EdgeColor','none');
xlabel('Distance to Leader (m)');
ylabel('Effort');
title('Distance vs Effort');
xlim([0 100]); ylim([0 1]); 
% grid on;

exportgraphics(gcf, fullfile(figDir, sprintf('Fig2_Distance_vs_TD_Combined_%s_pool.png', bandMode)), ...
               'Resolution',300);


%% -------- Figure 3：每人相关系数箱线图 --------
validR = ~isnan(r_DIST_TD);
figure(3); clf;
boxplot(r_DIST_TD(validR),'Symbol','','Labels', {''}); hold on;
yline(0,'r--','LineWidth',1);
ylabel('Correlation (r)'); title('Correlation of Distance vs Effort per Profile');
% grid on;

exportgraphics(gcf, fullfile(figDir,'Fig3_Boxplot_r_Distance_TD.png'), 'Resolution',300);

%% -------- 相关系数明细另存 --------
T = table((1:numP).', r_DIST_TD, p_DIST_TD, n_pairs, ...
          'VariableNames',{'Profile','r','p','N_pairs'});
writetable(T,'Distance_TD_correlation_per_profile.csv');
save('Distance_TD_correlation_per_profile.mat','T');
