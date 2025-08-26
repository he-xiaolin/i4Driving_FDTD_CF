clc; clear; close all;
load('participants_CF.mat');   % fields: THW, effort_scaled

%% --------- 参数 ----------
numP  = numel(participants_CF);
edges = linspace(0,10,121);
binCenters = edges(1:end-1) + diff(edges)/2;

% 每人相关系数（Figure 3 用）
r_THW_effort = nan(numP,1);
p_THW_effort = nan(numP,1);
n_pairs      = zeros(numP,1);

%% --------- Figure 1（隐藏）：64 个 profile 的 ribbon 全放在一张图 ----------
f1 = figure('Name','Per-profile THW vs Effort (8x8)','Units','normalized', ...
            'OuterPosition',[0 0 1 1],'Visible','off');

for k = 1:numP
    THW = participants_CF(k).THW(:);
    Eff = participants_CF(k).effort_scaled(:);

    % 相关系数
    v = ~isnan(THW) & ~isnan(Eff);
    n_pairs(k) = sum(v);
    if n_pairs(k) >= 2
        [r, p] = corr(THW(v), Eff(v), 'Type','Pearson');
        r_THW_effort(k) = r;  p_THW_effort(k) = p;
    end

    % 分箱统计
    [THWs, idx] = sort(THW);
    Effs = Eff(idx);
    mu = nan(size(binCenters)); sd = mu;
    for i = 1:numel(edges)-1
        m = (THWs >= edges(i)) & (THWs < edges(i+1));
        if any(m)
            mu(i) = mean(Effs(m), 'omitnan');
            sd(i) = std(Effs(m),  'omitnan');
        end
    end
    keep = ~isnan(mu); x = binCenters(keep); mu = mu(keep); sd = sd(keep);

    % 子图
    ax = subplot(8,8,k,'Parent',f1); hold(ax,'on');
    if ~isempty(x)
        fill(ax,[x, fliplr(x)], [mu - sd, fliplr(mu + sd)], ...
            'b','FaceAlpha',0.18,'EdgeColor','none');
        plot(ax,x,mu,'b-','LineWidth',1.2);
    end
    title(ax,sprintf('P%02d',k));
    xlim(ax,[0 10]); ylim(ax,[0 1]); grid(ax,'on');
end

exportgraphics(f1,'Figure1_PerProfile_THW_vs_Effort_8x8.png','Resolution',300);
close(f1);  % 不在屏幕上显示

%% --------- Figure 2：合并曲线（仅此图显示） ----------
exclude = [22 52 57 64]; % 不排除就设为 []
all_THW = []; all_Eff = [];
for k = 1:numP
    if ismember(k, exclude), continue; end
    all_THW = [all_THW; participants_CF(k).THW(:)];
    all_Eff = [all_Eff; participants_CF(k).effort_scaled(:)];
end

[THWs, idx] = sort(all_THW);
Effs = all_Eff(idx);

mu = nan(size(binCenters)); sd = mu;
for i = 1:numel(edges)-1
    m = (THWs >= edges(i)) & (THWs < edges(i+1));
    if any(m)
        mu(i) = mean(Effs(m), 'omitnan');
        sd(i) = std(Effs(m),  'omitnan');   % 注意：不除以 2
    end
end
keep = ~isnan(mu); x = binCenters(keep); mu = mu(keep); sd = sd(keep);

figure(2); clf;
plot(x,mu,'b-','LineWidth',2); hold on;
fill([x, fliplr(x)], [mu - sd, fliplr(mu + sd)], ...
     'b','FaceAlpha',0.2,'EdgeColor','none');
xlabel('THW'); ylabel('Effort Scaled');
title('THW vs Effort');
xlim([0 10]); ylim([0.1 0.8]); grid on;
exportgraphics(gcf,'Figure2_THW_vs_Effort_Combined.png','Resolution',300);

%% --------- Figure 3：箱线图（仅此图显示） ----------
validR = ~isnan(r_THW_effort);
figure(3); clf;
boxplot(r_THW_effort(validR),'Symbol','o'); hold on;
yline(0,'r--','LineWidth',1);
ylabel('Correlation (r)'); title('Correlation of THW vs Effort per Profile');
grid on;
exportgraphics(gcf,'Figure3_Boxplot_r_THW_Effort.png','Resolution',300);

% 相关系数明细保存
T = table((1:numP).', r_THW_effort, p_THW_effort, n_pairs, ...
          'VariableNames',{'Profile','r','p','N_pairs'});
writetable(T,'THW_Effort_correlation_per_profile.csv');
save('THW_Effort_correlation_per_profile.mat','T');
