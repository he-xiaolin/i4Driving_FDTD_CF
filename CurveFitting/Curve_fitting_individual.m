%% THW vs Effort — fast fit (exp vs linear), no inverse, no parallel
% - Per-driver: exponential (nonlinear, 1 param) via fminbnd + linear (1 param, OLS)
% - Pooled/global fit: same两模型
% - 模型比较：AIC/BIC + RMSE（R²仅作参考）
% - 图形：保持原8×8叠加曲线、R²>0筛选、残差网格、直方图/QQ、参数箱线图
% -------------------------------------------------------------------------
clc; clear; close all
load participants_CF_without_filter.mat   % 需要 participants_CF 结构，含 THW, effort_scaled
numParticipants = numel(participants_CF);

% ------------------------- trimming（仅用于可视化散点） -------------------------
trim.method      = 'sd';        % 'sd' 或 'percentile'
trim.centralPct  = 50;          % 保留中间 50%（≈ ±0.674σ）
trim.edgesTHW    = 0 : 0.25 : 10;
trim.minPts      = 5;           % 每个THW小bin内最少点数，不足则不trim
trim.minDrivers  = 5;           % 绘带所需的最少driver数
trim.ciType      = 'ci95';      % 'sem' 或 'ci95'

% ------------------------- 构建 DriverPlot & DriverCalc -------------------
DriverPlot = struct('THW',[],'Effort',[]);
DriverCalc = DriverPlot;  % 同字段，无trim，用于统计/拟合

for k = 1:numParticipants
    x = participants_CF(k).THW;
    y = participants_CF(k).effort_scaled;

    keep = local_buildKeepMask(x,y,trim);     % 仅用于可视化

    % 可视化用（trim 后 + 0-1 归一）
    DriverPlot(k).THW    = x(keep);
    DriverPlot(k).Effort = local_normalize01(y(keep));

    % 计算用（不 trim，只做 0-1 归一）
    DriverCalc(k).THW    = x;
    DriverCalc(k).Effort = local_normalize01(y);
end

% ------------------------- 绘图全局配置 -----------------------------------
cfg.markerSize = 6;
cfg.deepBlue   = 0.85*[0 0.4470 0.7410];
cfg.alphaBoost = 1.4;
cfg.xLim       = [0 10];
cfg.yLim       = [0  1];

% ------------------------- 汇总所有点（pooled） ---------------------------
allTHW    = vertcat(DriverCalc.THW);
allEffort = vertcat(DriverCalc.Effort);

%% =====================  每名被试：exp + linear 拟合  =====================
if ~exist('DriverCalc','var') || ~exist('DriverPlot','var')
    error('DriverCalc/DriverPlot 不存在。请检查上方构建步骤。');
end

nDrv = numel(DriverCalc);

% 结果容器
hExp   = nan(nDrv,1);  r2Exp  = nan(nDrv,1);  rmseExp = nan(nDrv,1);
aLin   = nan(nDrv,1);  r2Lin  = nan(nDrv,1);  rmseLin = nan(nDrv,1);

logLexp = nan(nDrv,1); AICexp = nan(nDrv,1); BICexp = nan(nDrv,1);
logLlin = nan(nDrv,1); AIClin = nan(nDrv,1); BIClin = nan(nDrv,1);

% 模型与辅助
expPred = @(h,x) exp(-x./h);
sse     = @(yhat,y) sum((y - yhat).^2);
stats   = @(yhat,y) deal( ...
              1 - sse(yhat,y)/max(eps, sum((y-mean(y)).^2)), ...
              sqrt(mean((y - yhat).^2)) );

lb = 0.05;  ub = 20;                       % 指数时间常数搜索范围 [s]
opts = optimset('TolX',1e-4,'Display','off');

for k = 1:nDrv
    x = DriverCalc(k).THW(:);
    y = DriverCalc(k).Effort(:);

    % --- exponential：单参数最小化 SSE
    obj = @(h) sse(expPred(h,x), y);
    hExp(k) = fminbnd(obj, lb, ub, opts);
    yHatExp = expPred(hExp(k), x);
    [r2Exp(k), rmseExp(k)] = stats(yHatExp, y);
    [logLexp(k), AICexp(k), BICexp(k)] = normalLogLike(y, yHatExp, 1); % k=1

    % --- linear：E = 1 + a*x （闭式解）
    denom = sum(x.^2);
    if denom<=eps, aLin(k)=0; else, aLin(k) = sum((y-1).*x) / denom; end
    yHatLin = 1 + aLin(k)*x;
    [r2Lin(k), rmseLin(k)] = stats(yHatLin, y);
    [logLlin(k), AIClin(k), BIClin(k)] = normalLogLike(y, yHatLin, 1); % k=1
end

% 汇总结果表
FitTable = table((1:nDrv).', ...
    hExp, r2Exp, rmseExp, logLexp, AICexp, BICexp, ...
    aLin, r2Lin, rmseLin, logLlin, AIClin, BIClin, ...
    'VariableNames',{'Driver', ...
        'h0_exp','R2_exp','RMSE_exp','logL_exp','AIC_exp','BIC_exp', ...
        'a_lin','R2_lin','RMSE_lin','logL_lin','AIC_lin','BIC_lin'});

%% =====================  模型比较（每名被试）  ============================
nDrv = height(FitTable);

% 差值：lin - exp；AIC/BIC 越小越好 ⇒ Δ>0 表示 exp 更优
FitTable.dRMSE = FitTable.RMSE_lin - FitTable.RMSE_exp;
FitTable.dAIC  = FitTable.AIC_lin  - FitTable.AIC_exp;
FitTable.dBIC  = FitTable.BIC_lin  - FitTable.BIC_exp;

% AIC 判优（|ΔAIC|<2 视为证据不足 → tie）
FitTable.winner_AIC = repmat("linear", nDrv, 1);
FitTable.winner_AIC(FitTable.dAIC > 0) = "exp";
FitTable.winner_AIC(abs(FitTable.dAIC) < 2) = "tie";

% AIC 证据强度分箱（|ΔAIC|）
edges  = [0 2 4 10 inf];
labels = {'≤2 (dubious)','2–4 (weak)','4–10 (considerable)','>10 (strong)'};

idx = discretize(abs(FitTable.dAIC), edges);                 % N×1
FitTable.AIC_strength = categorical(idx(:), 1:numel(labels), labels);  % N×1


% 打印小结
nExp_AIC = sum(FitTable.dAIC> 2);
nLin_AIC = sum(FitTable.dAIC<-2);
nTie_AIC = sum(abs(FitTable.dAIC)<2);
fprintf('AIC (|Δ|>2): exp %d, linear %d, ties %d (N=%d)\n', ...
    nExp_AIC, nLin_AIC, nTie_AIC, nDrv);
fprintf('RMSE: exp better %d, linear better %d, ties %d (|ΔRMSE|<=1e-6)\n', ...
    sum(FitTable.dRMSE>1e-6), sum(FitTable.dRMSE<-1e-6), sum(abs(FitTable.dRMSE)<=1e-6));

% 保存结果表
save   effort_model_fits.mat  FitTable
writetable(FitTable,'effort_model_fits.csv')
fprintf('Saved fit results (exp vs linear) to effort_model_fits.(mat|csv)\n');

%% ==================  8×8 叠加曲线（两模型） =============================
figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(8,8,'Padding','compact','TileSpacing','compact');
for k = 1:nDrv
    x = DriverPlot(k).THW;   y = DriverPlot(k).Effort;
    alfa = fastAlpha(x,y)*1.4; alfa(alfa>1)=1;
    ax = nexttile(tl);
    h1 = scatter(x,y,6,[0 0.45 0.74],'filled', ...
            'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
            'AlphaData',alfa,'AlphaDataMapping','none'); hold on
    xx = linspace(0,10,200);
    h2 = plot(xx, exp(-xx./hExp(k)) ,'b' ,'LineWidth',1);
    h4 = plot(xx, 1 + aLin(k)*xx    ,'k--','LineWidth',1);
    xlim([0 10]); ylim([0 1]); grid on
    title(sprintf('P %d',k));
end
lgd = legend([h1 h2 h4], {'samples','exp','linear'}, ...
             'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'north';
sgtitle(tl, 'THW vs Effort curve fitting', 'FontWeight','bold');
xlabel(tl,  'Time Headway (s)');
ylabel(tl,  'Effort (normalized)');
tl.Padding = 'loose'; tl.TileSpacing = 'compact';

%% ============== R²>0 的被试（只画保留者，仍为两模型） ====================
% bestR2  = max([FitTable.R2_exp, FitTable.R2_lin],[],2);
% keepMask = bestR2 > 0;
keepMask = FitTable.R2_exp>0;

figure('units','normalized','outerposition',[0 0 1 1]);
tlKeep = tiledlayout(8,8,'Padding','compact','TileSpacing','compact');
hScatter=[]; hExpLine=[]; hLinLine=[];
for k = 1:64
    ax = nexttile(tlKeep);
    if keepMask(k)
        x = DriverPlot(k).THW;   y = DriverPlot(k).Effort;
        alfa = fastAlpha(x,y)*1.4;  alfa(alfa>1)=1;

        h1 = scatter(ax,x,y,6,[0 0.45 0.74],'filled', ...
                'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
                'AlphaData',alfa,'AlphaDataMapping','none'); hold(ax,'on')
        xx = linspace(0,10,200);
        h2 = plot(ax,xx, exp(-xx./hExp(k)),'b' ,'LineWidth',1);
        h4 = plot(ax,xx, 1 + aLin(k)*xx  ,'k--','LineWidth',1);
        xlim(ax,[0 10]); ylim(ax,[0 1]); grid(ax,'on')
        title(ax,sprintf('P %d',k));

        if isempty(hScatter), hScatter=h1; hExpLine=h2; hLinLine=h4; end
    else
        axis(ax,'off')
        text(0.5,0.5,{'excluded','(R^2 ≤ 0)'}, ...
             'Units','normalized','HorizontalAlignment','center', ...
             'FontSize',7,'Color',[.5 .5 .5],'Parent',ax);
        title(ax,sprintf('P %d',k));
    end
end
lgd = legend([hScatter hExpLine hLinLine], {'samples','exp','linear'}, ...
             'Orientation','horizontal','NumColumns',3);
lgd.Layout.Tile = 'north';
sgtitle(tlKeep,'THW vs Effort – profiles with positive R²','FontWeight','bold');
xlabel(tlKeep,'Time Headway (s)'); ylabel(tlKeep,'Effort (normalized)');
tlKeep.Padding='loose'; tlKeep.TileSpacing='compact';

%% ============================ 残差与正态性 ===============================
% 1) 每被试残差
resExp = cell(nDrv,1); resLin = cell(nDrv,1);
for k = 1:nDrv
    x = DriverCalc(k).THW(:);
    y = DriverCalc(k).Effort(:);
    resExp{k} = y - exp(-x./hExp(k));
    resLin{k} = y - (1 + aLin(k)*x);
end
allExp = vertcat(resExp{:}); allLin = vertcat(resLin{:});

% 8×8 残差散点格
plotResidualGrid(resExp, '\epsilon  =  y -  y_{exp}',   cfg);
plotResidualGrid(resLin, '\epsilon  =  y -  y_{linear}',cfg);

% 全局直方图 + QQ
figure('units','normalized','outerposition',[0 0 1 0.4]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
nexttile, histogram(allExp,40,'FaceColor',[0.2 0.6 1]); hold on
xline(mean(allExp),'k','LineWidth',1.5); title('Residuals – exp model')
nexttile, histogram(allLin,40,'FaceColor',[0.2 0.2 0.2]); hold on
xline(mean(allLin),'w','LineWidth',1.5); title('Residuals – linear model')

figure('units','normalized','outerposition',[0 0 1 0.4]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
qqplotPanel(allExp,'exp'); qqplotPanel(allLin,'linear');

% R²>0 的被试：残差可视化
keepIdx = (FitTable.R2_exp > 0) | (FitTable.R2_lin > 0);
if any(keepIdx)
    resExpFull = cell(64,1); resLinFull = cell(64,1);
    for k = 1:64
        if keepIdx(k)
            x = DriverCalc(k).THW(:);
            y = DriverCalc(k).Effort(:);
            resExpFull{k} = y - exp(-x ./ hExp(k));
            resLinFull{k} = y - (1 + aLin(k) * x);
        end
    end
    plotResidualGridFull(resExpFull, keepIdx, '\epsilon = y - y_{exp}',   cfg);
    plotResidualGridFull(resLinFull, keepIdx, '\epsilon = y - y_{linear}',cfg);

    allExpK = vertcat(resExpFull{keepIdx}); allLinK = vertcat(resLinFull{keepIdx});
    figure('units','normalized','outerposition',[0 0 1 0.4]);
    tlH = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile(tlH), histogram(allExpK,40,'FaceColor',[0.2 0.6 1]); hold on
    xline(mean(allExpK),'k','LineWidth',1.5); title('Residuals – exp model'); grid on
    nexttile(tlH), histogram(allLinK,40,'FaceColor',[0.2 0.2 0.2]); hold on
    xline(mean(allLinK),'w','LineWidth',1.5); title('Residuals – linear model'); grid on
    sgtitle(tlH,'Histograms of residuals  (kept drivers)','FontWeight','bold')

    figure('units','normalized','outerposition',[0 0 1 0.4]);
    tlQ = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    qqplotPanelnew(allExpK,'exp', tlQ);
    qqplotPanelnew(allLinK,'linear', tlQ);
    sgtitle(tlQ,'QQ-plots of residuals  (kept drivers)','FontWeight','bold')
end

%% ====================== 参数分布（箱线图，仅R²>0） =======================
mask = keepMask;
h0exp = FitTable.h0_exp(mask);
aLinK = FitTable.a_lin(mask);

figure('units','normalized','outerposition',[0 0 0.6 0.45]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile(tl,1)
boxplot(h0exp,'Notch','off','Labels',{'h_{exp}'},'Whisker',1.5)
set(gca,'YScale','log'); ylabel('h_{exp}  [s]'); title('exp model')
nexttile(tl,2)
boxplot(aLinK ,'Notch','off','Labels',{'slope  k'})
yline(0,'k--'); ylabel('k'); title('linear model')
sgtitle(tl,'Fitted-parameter distributions','FontWeight','bold');

%% =============================  保存全部图  ==============================
if ~exist('./curveFittingFigures','dir'), mkdir('./curveFittingFigures'); end
allFigs = findall(0,'Type','figure');
for i = 1:length(allFigs)
    f = ancestor(allFigs(i),'figure');
    if ~isempty(f)
        filename = sprintf('./curveFittingFigures/Figure_%d.jpg', i);
        exportgraphics(f, filename, 'Resolution',300);
    end
end

%% ============================ GLOBAL（pooled） ===========================
sse   = @(yhat,y) sum((y - yhat).^2);
expFx = @(h,x) exp(-x./h);
stats = @(yhat,y) deal(1 - sse(yhat,y)/max(eps,sum((y-mean(y)).^2)), ...
                       sqrt(mean((y - yhat).^2)));
lb=0.05; ub=20; opts=optimset('TolX',1e-4,'Display','off');

% (1) exp
hExpAll = fminbnd(@(h) sse(expFx(h,allTHW),allEffort), lb, ub, opts);
[R2ExpAll, RMSEExpAll] = stats(expFx(hExpAll,allTHW), allEffort);

% (2) linear
aLinAll = sum((allEffort-1).*allTHW) / max(eps,sum(allTHW.^2));
[R2LinAll, RMSELinAll] = stats(1 + aLinAll*allTHW, allEffort);

fprintf('\nGLOBAL fits (all participants)\n');
fprintf('  exp:    h0 = %.3g  |  R² = %.3f  |  RMSE = %.3f\n', hExpAll, R2ExpAll, RMSEExpAll);
fprintf('  linear:   a = %.3g |  R² = %.3f  |  RMSE = %.3f\n', aLinAll, R2LinAll, RMSELinAll);

% 全局 AIC/BIC 比较
[logLexpAll, AICexpAll, BICexpAll] = normalLogLike(allEffort, expFx(hExpAll,allTHW), 1);
[logLlinAll, AIClinAll, BIClinAll] = normalLogLike(allEffort, 1 + aLinAll*allTHW, 1);
dAIC_all = AIClinAll - AICexpAll;   % >0 => exp 更好
dBIC_all = BIClinAll - BICexpAll;
whoBetter = "tie";
if     dAIC_all >  2, whoBetter = "exp";
elseif dAIC_all < -2, whoBetter = "linear"; end
fprintf('GLOBAL ΔAIC (lin-exp) = %.2f  -> %s better by AIC (|Δ|>2 rule)\n', dAIC_all, whoBetter);
fprintf('GLOBAL ΔBIC (lin-exp) = %.2f\n', dBIC_all);

% 全局散点 + 两模型曲线
figure('units','normalized','outerposition',[0 0 1 1]);
sh = randperm(numel(allTHW), min(6e4,numel(allTHW)));   % 可视化下采样
alphaAll = fastAlpha(allTHW, allEffort);
scatter(allTHW(sh), allEffort(sh), 6, cfg.deepBlue, 'filled', ...
        'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
        'AlphaData',alphaAll(sh),'AlphaDataMapping','none'); hold on
xx = linspace(0,10,300);
plot(xx, expFx(hExpAll,xx) ,'b' ,'LineWidth',2);
plot(xx, 1 + aLinAll*xx    ,'k--','LineWidth',2);
xlim([0 10]); ylim([0 1]); grid on
legend({'samples','exp','linear'},'Location','southwest');
title('THW vs Effort – global fits (all drivers)');
xlabel('Time Headway (s)');  ylabel('Effort (normalized)');

%% ================================ 函数区 ================================
function plotUniformAlpha(D, cfg)
    figure('units','normalized','outerposition',[0 0 1 1]);
    for k = 1:numel(D)
        subplot(8,8,k)
        x = D(k).THW;  y = D(k).Effort;
        alpha = fastAlpha(x,y)*cfg.alphaBoost; alpha(alpha>1)=1;
        scatter(x,y,cfg.markerSize,'filled', ...
            'MarkerFaceColor',cfg.deepBlue,'MarkerEdgeColor',cfg.deepBlue, ...
            'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
            'AlphaData',alpha,'AlphaDataMapping','none');
        xlim(cfg.xLim); ylim(cfg.yLim); grid on; title(sprintf('P %d',k));
    end
    set(axes('Visible','off'), ...
        'Title',text('String','THW vs Effort (opacity ∝ density)','Visible','on'), ...
        'XLabel',text('String','Time Headway (s)','Visible','on'), ...
        'YLabel',text('String','Effort (normalized)','Visible','on'));
end

function plotColourScatter(D, cfg, field, labelTxt)
    figure('units','normalized','outerposition',[0 0 1 1]);
    allVals = [];
    for k = 1:numel(D)
        subplot(8,8,k)
        x = D(k).THW;   y = D(k).Effort;
        cRaw = evalin('caller',['participants_CF(' num2str(k) ').' field]);
        c    = cRaw(1:numel(x));
        alpha = fastAlpha(x,y);
        scatter(x,y,cfg.markerSize,c,'filled', ...
            'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
            'AlphaData',alpha,'AlphaDataMapping','none');
        xlim(cfg.xLim); ylim(cfg.yLim); grid on; title(sprintf('P %d',k));
        allVals = [allVals ; c];
    end
    han = axes('Visible','off');
    set([han.Title han.XLabel han.YLabel],'Visible','on');
    title(han,['THW vs Effort coloured by ' labelTxt]);
    xlabel(han,'Time Headway (s)'); ylabel(han,'Effort (normalized)');
    colormap(parula); clim([min(allVals) max(allVals)]);
    cb = colorbar(han,'eastoutside'); cb.Label.String = labelTxt;
    cb.Position = [0.935 0.05 0.015 0.90];
end

function plotGroupMeanRibbon(Dplot, Dcalc, trim, cfg)
    % 可选：群体均值带（未在主流程调用，保留以备需要）
    allX = cat(1,Dplot.THW);   allY = cat(1,Dplot.Effort);
    alpha = fastAlpha(allX, allY);

    edges   = trim.edgesTHW;  centers = edges(1:end-1)+diff(edges)/2;
    nBins   = numel(centers);  nDrv = numel(Dcalc);
    meanMat = nan(nDrv,nBins);

    for k = 1:nDrv
        x = Dcalc(k).THW;  y = Dcalc(k).Effort;
        for i = 1:nBins
            m = x>=edges(i) & x<edges(i+1);
            if any(m), meanMat(k,i)=mean(y(m)); end
        end
    end

    nPerBin = sum(~isnan(meanMat),1);
    gMean   = nanmean(meanMat,1);
    gSEM    = nanstd(meanMat,0,1) ./ sqrt(nPerBin);
    ciMult  = 1.96*strcmpi(trim.ciType,'ci95') + strcmpi(trim.ciType,'sem');
    lo      = gMean - ciMult*gSEM;    hi = gMean + ciMult*gSEM;
    valid   = nPerBin >= trim.minDrivers;

    figure('units','normalized','outerposition',[0 0 1 1]);
    sh = randperm(numel(allX), min(4e4,numel(allX)));
    scatter(allX(sh), allY(sh), 8, cfg.deepBlue,'filled', ...
        'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
        'AlphaData',alpha(sh),'AlphaDataMapping','none'); hold on
    fill([centers(valid) fliplr(centers(valid))], ...
         [lo(valid)      fliplr(hi(valid))], ...
         'b','FaceAlpha',0.15,'EdgeColor','none');
    plot(centers(valid), gMean(valid),'b-','LineWidth',2);
    xlabel('THW (s)'); ylabel('Effort (normalized)');
    title(['Group mean Effort ± ' upper(trim.ciType)]);
    xlim(cfg.xLim); ylim(cfg.yLim); grid on
    legend({'samples','precision band','group mean'},'Location','southwest');
end

function keep = local_buildKeepMask(x,y,t)
    keep = false(size(x));
    if strcmpi(t.method,'sd'), k=sqrt(2)*erfinv(t.centralPct/100); end
    for b=1:numel(t.edgesTHW)-1
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
    d = max(y) - min(y);
    if d > 0
        y_n = (y - min(y)) ./ d;   % 0–1
    else
        y_n = 0.5*ones(size(y));   % 常量向量 → 0.5
    end
end

function a = fastAlpha(x,y)
    ok = ~(isnan(x) | isnan(y)); x = x(ok); y = y(ok);
    if numel(x) < 2, a = 0.5*ones(size(x)); return, end
    nb = 50;
    if range(x)==0, ex=[x(1)-0.01, x(1)+0.01]; else, ex=linspace(min(x),max(x),nb+1); end
    if range(y)==0, ey=[y(1)-0.01, y(1)+0.01]; else, ey=linspace(min(y),max(y),nb+1); end
    [N,~,~,ix,iy] = histcounts2(x,y,ex,ey);
    lin = sub2ind(size(N), ix, iy);
    a   = rescale(N(lin), 0.05, 1.0);
end

function plotResidualGrid(resCell, yLabelTxt, cfg)
    figure('units','normalized','outerposition',[0 0 1 1]);
    for k = 1:numel(resCell)
        subplot(8,8,k)
        n = numel(resCell{k});
        scatter(1:n, resCell{k}, 6, cfg.deepBlue, 'filled',...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
        yline(0,'k--'); ylim([-0.6 0.6]); xlim([0 max(1,n)+1]); grid on
        title(sprintf('P %d',k));
    end
    axHost = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5, 0.98, ['Residuals grid  (' yLabelTxt ')'], ...
         'Units','normalized','HorizontalAlignment','center','FontWeight','bold','Parent',axHost);
    text(0.5, -0.04, 'Sample index','Units','normalized','HorizontalAlignment','center','Parent',axHost);
    text(-0.04, 0.5, yLabelTxt,'Units','normalized','HorizontalAlignment','center','Rotation',90,'Parent',axHost);
end

function qqplotPanel(r, modelName)
    nexttile, qqplot(r); grid on, title(['QQ-plot  (' modelName ')'])
end

function plotResidualGridFull(resCell, keepMask, yLabelTxt, cfg)
    figure('units','normalized','outerposition',[0 0 1 1]);
    nRows=8; nCols=8;
    for k = 1:64
        subplot(nRows, nCols, k)
        if keepMask(k)
            r = resCell{k};
            scatter(1:numel(r), r, 6, cfg.deepBlue, 'filled', ...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6); hold on
            yline(0,'k--'); ylim([-0.6 0.6]); xlim([0 max(1,numel(r))+1]); grid on
        else
            axis off
            text(0.5,0.5,{'excluded','(R^2 \le 0)'},'HorizontalAlignment','center','FontSize',7,'Color',[.4 .4 .4]);
        end
        title(sprintf('P %d', k));
    end
    axHost = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5, 0.98, ['Residuals grid  (' yLabelTxt ')'], ...
         'Units','normalized','HorizontalAlignment','center','FontWeight','bold','Parent',axHost);
    text(0.5, -0.04, 'Sample index','Units','normalized','HorizontalAlignment','center','Parent',axHost);
    text(-0.04, 0.5, yLabelTxt,'Units','normalized','HorizontalAlignment','center','Rotation',90,'Parent',axHost);
end

function qqplotPanelnew(resVec, labelStr, parentTL)
    nexttile(parentTL), qqplot(resVec); grid on, title(['QQ-plot (' labelStr ')'])
end

function [logL, AIC, BIC] = normalLogLike(y, yHat, k)
    n      = numel(y);
    SSE    = sum((y - yHat).^2);
    sigma2 = SSE / max(1,n);                  % MLE of variance
    logL   = -0.5*n * ( log(2*pi*sigma2) + 1 );
    AIC    = 2*k - 2*logL;                    % k = #parameters (=1 here)
    BIC    = k*log(max(1,n)) - 2*logL;
end
