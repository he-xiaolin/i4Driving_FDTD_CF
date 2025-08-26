clc; clear; close all;
load('participants_CF.mat');   % 需要字段: participants_CF(k).THW, .effort_scaled

%% ---------------- 参数区 ----------------
numP  = numel(participants_CF);
edges = linspace(0,10,121);                 % THW 分箱
binCtr = edges(1:end-1) + diff(edges)/2;

normalizeEffort = true;                     % 是否对每个被试的 effort 做 min-max
bandMode   = 'sd';                         % 'sd' | 'sem' | 'iqr' | 'p90'   <<< 改这里即可切换到 SD
smoothWin  = 1;                             % 跨被试聚合后的移动平均窗口 (bin)
minProfilesPerBin = 5;                      % 至少多少个被试在该 bin 有数据才绘制
exclude = [22 52 57 64];                    % 合并曲线时排除的人(可设为空 [])

% 统一的图像输出文件夹
figDir = fullfile(pwd,'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

%% -------- 容器（Figure 3 用；以及跨被试聚合用） --------
r_THW_effort = nan(numP,1);
p_THW_effort = nan(numP,1);
n_pairs      = zeros(numP,1);

B = numel(edges)-1;
muPerP = nan(numP, B);                      % 每个被试在每个 bin 的均值 (min-max 后)

%% -------- Figure 1（隐藏）：8×8 单人 ribbon，min-max 后 --------
f1 = figure('Name','Per-profile THW vs Effort (8x8)','Units','normalized', ...
            'OuterPosition',[0 0 1 1],'Visible','off');

for k = 1:numP
    THW = participants_CF(k).THW(:);
    Eff = participants_CF(k).effort_scaled(:);

    % 相关系数（Pearson；与是否 min-max 无关）
    v = isfinite(THW) & isfinite(Eff);
    n_pairs(k) = sum(v);
    if n_pairs(k) >= 2
        [R,P] = corrcoef(THW(v), Eff(v));
        r_THW_effort(k) = R(1,2);
        p_THW_effort(k) = P(1,2);
    end

    % 每个被试对 effort 做 min-max（仅用于绘图与跨被试聚合）
    if normalizeEffort
        eMin = min(Eff,[],'omitnan'); eMax = max(Eff,[],'omitnan');
        if ~isnan(eMin) && ~isnan(eMax) && eMax > eMin
            Eff = (Eff - eMin) / (eMax - eMin);
        else
            Eff(:) = NaN;
        end
    end

    % 分箱统计（个人）
    [THWs, idx] = sort(THW);
    Effs = Eff(idx);
    mu = nan(1,B); sd = nan(1,B);
    for i = 1:B
        m = (THWs >= edges(i)) & (THWs < edges(i+1));
        if any(m)
            mu(i) = mean(Effs(m), 'omitnan');
            sd(i) = std(  Effs(m), 'omitnan');
        end
    end
    muPerP(k,:) = mu;                       % 存起来供后续“被试层面”聚合

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

exportgraphics(f1, fullfile(figDir,'Fig1_PerProfile_THW_vs_Effort_8x8.png'), 'Resolution',300);
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

% 可选平滑
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
xlabel('THW (s)');
ylabel('Effort');

title('THW vs Effort');
xlim([0 10]); ylim([0 1]); 
% grid on;

exportgraphics(gcf, fullfile(figDir, sprintf('Fig2_THW_vs_Effort_Combined_%s.png', bandMode)), ...
               'Resolution',300);

%% -------- Figure 3：每人相关系数箱线图 --------
validR = ~isnan(r_THW_effort);
figure(3); clf;
boxplot(r_THW_effort(validR),'Symbol','', 'Labels', {''}); hold on;
yline(0,'r--','LineWidth',1);
ylabel('Correlation (r)'); title('Correlation of THW vs Effort per Profile');
% grid on;

exportgraphics(gcf, fullfile(figDir,'Fig3_Boxplot_r_THW_Effort.png'), 'Resolution',300);

%% -------- 相关系数明细另存（非图片，可选放根目录） --------
T = table((1:numP).', r_THW_effort, p_THW_effort, n_pairs, ...
          'VariableNames',{'Profile','r','p','N_pairs'});
writetable(T,'THW_Effort_correlation_per_profile.csv');
save('THW_Effort_correlation_per_profile.mat','T');
