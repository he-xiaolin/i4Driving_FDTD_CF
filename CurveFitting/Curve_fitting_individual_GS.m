%% plot_THW_vs_effort_density_alpha_FAST.m
% Clean visuals vs. full-data calculations
% ---------------------------------------
% • DriverPlot : trimmed (central 50 %) + 0-1 scaling – for scatter density
% • DriverCalc : no trimming, 0-1 scaling only        – for ribbon / stats
% Uses a fast 2-D histogram for alpha (no mvksdensity).

clc; clear; close all
load participants_CF_without_filter.mat
numParticipants = numel(participants_CF);

% -------------------------  trimming policy (for visuals only)
trim.method      = 'sd';      % 'sd'  or 'percentile'
trim.centralPct  = 50;        % keep central 50 %  (≈ ±0.674 σ)
trim.edgesTHW    = 0 : 0.25 : 10;
trim.minPts      = 5;         % leave tiny bins untrimmed
trim.minDrivers  = 5;         % bins in ribbon need ≥ this many drivers
trim.ciType      = 'ci95';    % 'sem' or 'ci95'

% -------------------------  build DriverPlot & DriverCalc --------------
DriverPlot = struct('THW',[],'Effort',[]);
DriverCalc = DriverPlot;                      % same fields

for k = 1:numParticipants
    x = participants_CF(k).THW;
    y = participants_CF(k).effort_scaled;

    keep = local_buildKeepMask(x,y,trim);     % ±0.67 σ   (visual only)

    % ----- for scatter visuals (trimmed)
    DriverPlot(k).THW    = x(keep);
    DriverPlot(k).Effort = local_normalize01(y(keep));

    % ----- for calculations (no trimming)
    DriverCalc(k).THW    = x;
    DriverCalc(k).Effort = local_normalize01(y);
end

% -------------------------  plotting settings --------------------------
cfg.markerSize = 6;
cfg.deepBlue   = 0.85*[0 0.4470 0.7410];
cfg.alphaBoost = 1.4;
cfg.xLim       = [0 10];
cfg.yLim       = [0  1];

% % -------------------------  1) 8×8 uniform scatter ---------------------
% plotUniformAlpha(DriverPlot, cfg);
% 
% % -------------------------  2) 8×8 colour scatters ---------------------
% plotColourScatter(DriverPlot, cfg,'relative_speed',      'Relative speed (m/s)');
% plotColourScatter(DriverPlot, cfg,'speed_ego',           'Ego speed (m/s)');
% plotColourScatter(DriverPlot, cfg,'distance_leader',     'Distance to leader (m)');
% plotColourScatter(DriverPlot, cfg,'acceleration_leader', 'Leader acceleration (m/s²)');
% plotColourScatter(DriverPlot, cfg,'relative_acceleration','Relative acceleration (m/s²)');
% 
% % -------------------------  3) combined cloud + ribbon -----------------
% plotGroupMeanRibbon(DriverPlot, DriverCalc, trim, cfg);
%% all points collection
% ---------- pooled data --------------------------------------------------
allTHW    = vertcat(DriverCalc.THW);     %  N×1  vector
allEffort = vertcat(DriverCalc.Effort);  %  N×1  vector (already 0–1 scaled)


%% Curve fitting Grid search
%% fit_effort_models.m  (grid-search version)
% -------------------------------------------------------------------------
% Fits three effort–THW laws to every participant.
%   Model 1 : E = exp(-THW / h0)         (grid search)
%   Model 2 : E = h0 / (h0 + THW)        (grid search)
%   Model 3 : E = 1 + a·THW              (analytic OLS)
% Uses DriverCalc (un-trimmed, 0-1 scaled) for fitting.
% -------------------------------------------------------------------------

if ~exist('DriverCalc','var') || ~exist('DriverPlot','var')
    error('Run Scatters_with_coloured_features_updated.m first.');
end

nDrv = numel(DriverCalc);

% ---------- allocate result arrays --------------------------------------
hExp = nan(nDrv,1);  r2Exp = nan(nDrv,1);  rmseExp = nan(nDrv,1);
hRat = nan(nDrv,1);  r2Rat = nan(nDrv,1);  rmseRat = nan(nDrv,1);
aLin = nan(nDrv,1);  r2Lin = nan(nDrv,1);  rmseLin = nan(nDrv,1);

logLexp = nan(nDrv,1); AICexp = nan(nDrv,1); BICexp = nan(nDrv,1);
logLrat = nan(nDrv,1); AITrat = nan(nDrv,1); BITrat = nan(nDrv,1);
logLlin = nan(nDrv,1); AIClin = nan(nDrv,1); BIClin = nan(nDrv,1);

% ---------- coarse grid (log) & helper lambdas --------------------------
gridCoarse = logspace(log10(0.05), log10(20), 40);   % 0.05–20 s

expPred = @(h,x) exp(-x./h);
ratPred = @(h,x) h./(h + x);
sse     = @(yhat,y) sum((y - yhat).^2);
stats   = @(yhat,y) deal(1 - sse(yhat,y)/sum((y-mean(y)).^2), ...
                         sqrt(mean((y - yhat).^2)));          % R2, RMSE

% ---------- loop over drivers -------------------------------------------
for k = 1:nDrv
    x = DriverCalc(k).THW(:);
    y = DriverCalc(k).Effort(:);

    % ----- exponential : coarse search ----------------------------------
    sseExp = arrayfun(@(h) sse(expPred(h,x),y), gridCoarse);
    [~,idx] = min(sseExp);  h0Start = gridCoarse(idx);

    % refine ±25 % linearly around best coarse point
    gridFine = linspace(0.75*h0Start, 1.25*h0Start, 60);
    sseFine  = arrayfun(@(h) sse(expPred(h,x),y), gridFine);
    [~,idx2] = min(sseFine);  hExp(k) = gridFine(idx2);
    [r2Exp(k), rmseExp(k)] = stats(expPred(hExp(k),x), y);

    yHatExp              = expPred(hExp(k),x);
[r2Exp(k), rmseExp(k)] = stats(yHatExp, y);
[logLexp(k), AICexp(k), BICexp(k)] = normalLogLike(y, yHatExp, 1);
    % ----- rational : same two-stage search -----------------------------
    sseRat = arrayfun(@(h) sse(ratPred(h,x),y), gridCoarse);
    [~,idx] = min(sseRat);   h0Start = gridCoarse(idx);
    gridFine = linspace(0.75*h0Start, 1.25*h0Start, 60);
    sseFine  = arrayfun(@(h) sse(ratPred(h,x),y), gridFine);
    [~,idx2] = min(sseFine); hRat(k) = gridFine(idx2);
    [r2Rat(k), rmseRat(k)] = stats(ratPred(hRat(k),x), y);
    yHatRat              = ratPred(hRat(k),x);
[r2Rat(k), rmseRat(k)] = stats(yHatRat, y);
[logLrat(k), AITrat(k), BITrat(k)] = normalLogLike(y, yHatRat, 1);
    % ----- linear : analytic OLS ---------------------------------------
    aLin(k) = sum((y-1).*x) / sum(x.^2);
    yHatLin = 1 + aLin(k)*x;
    [r2Lin(k), rmseLin(k)] = stats(yHatLin, y);
    yHatLin              = 1 + aLin(k)*x;
[r2Lin(k), rmseLin(k)] = stats(yHatLin, y);
    [logLlin(k), AIClin(k), BIClin(k)] = normalLogLike(y, yHatLin, 1);
end

% ---------- compile & save table ----------------------------------------
FitTable = table((1:nDrv).', ...
    hExp, r2Exp, rmseExp, logLexp, AICexp, BICexp, ...
    hRat, r2Rat, rmseRat, logLrat, AITrat, BITrat, ...
    aLin, r2Lin, rmseLin, logLlin, AIClin, BIClin, ...
    'VariableNames',{'Driver', ...
        'h0_exp','R2_exp','RMSE_exp','logL_exp','AIC_exp','BIC_exp', ...
        'h0_rat','R2_rat','RMSE_rat','logL_rat','AIC_rat','BIC_rat', ...
        'a_lin','R2_lin','RMSE_lin','logL_lin','AIC_lin','BIC_lin'});

save   effort_model_fits.mat  FitTable
writetable(FitTable,'effort_model_fits.csv')
fprintf('Saved fit results (grid search) to effort_model_fits.(mat|csv)\n');

% ---------- overlay curves on 8×8 scatter grid --------------------------
figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(8,8,'Padding','compact','TileSpacing','compact');
for k = 1:nDrv
    % subplot(8,8,k)
    x = DriverPlot(k).THW;   y = DriverPlot(k).Effort;
    alfa = fastAlpha(x,y)*1.4; alfa(alfa>1)=1;
    ax = nexttile(tl);
    h1 = scatter(x,y,6,[0 0.45 0.74],'filled', ...
            'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
            'AlphaData',alfa,'AlphaDataMapping','none'); hold on
    xx = linspace(0,10,200);
    h2 = plot(xx, expPred(hExp(k),xx),'b' ,'LineWidth',1);
    h3 = plot(xx, ratPred(hRat(k),xx),'r' ,'LineWidth',2);
    h4 = plot(xx, 1+aLin(k)*xx      ,'k--','LineWidth',1);
    xlim([0 10]); ylim([0 1]); grid on
    title(sprintf('P %d',k));
end
% set(axes('Visible','off'), ...
%     'XLabel',text('String','Time Headway (s)','Visible','on'), ...
%     'YLabel',text('String','Effort (normalized)','Visible','on'));
% ONE call – pass only the graphic handles and labels
lgd = legend([h1 h2 h3 h4], {'samples','exp','inverse','linear'}, ...
             'Orientation','horizontal','NumColumns',4);

% tell MATLAB to place the legend in the reserved row of the tiled layout
lgd.Layout.Tile = 'north';
% sgtitle('THW vs Effort curve fitting')
    % 'Title',text('String','THW vs Effort – fitted curves (grid search)','Visible','on'), ...
% global title and labels
sgtitle(tl, 'THW vs Effort curve fitting', 'FontWeight','bold');
xlabel(tl,  'Time Headway (s)');
ylabel(tl,  'Effort (normalized)');

% nudge their apparent position by adjusting layout padding
tl.Padding = 'loose';          % more breathing room
tl.TileSpacing = 'compact';
%%
%% ------------------------------------------------------------------------
% Extra figure – keep only drivers whose BEST R² > 0
% ------------------------------------------------------------------------
bestR2 = max([FitTable.R2_exp, FitTable.R2_rat, FitTable.R2_lin],[],2);
keepMask = bestR2 > 0;          % 64×1 logical

figure('units','normalized','outerposition',[0 0 1 1]);
tlKeep = tiledlayout(8,8,'Padding','compact','TileSpacing','compact');

% store handles once for the legend
hScatter = []; hExpLine = []; hRatLine = []; hLinLine = [];

for k = 1:64
    ax = nexttile(tlKeep);

    if keepMask(k)                         % -------- kept driver --------
        x = DriverPlot(k).THW;   y = DriverPlot(k).Effort;
        alfa = fastAlpha(x,y)*1.4;  alfa(alfa>1)=1;

        h1 = scatter(ax,x,y,6,[0 0.45 0.74],'filled', ...
                'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
                'AlphaData',alfa,'AlphaDataMapping','none'); hold(ax,'on')

        xx = linspace(0,10,200);
        h2 = plot(ax,xx, exp(-xx./hExp(k))       ,'b' ,'LineWidth',1);
        h3 = plot(ax,xx, hRat(k)./(hRat(k)+xx)   ,'r' ,'LineWidth',2);
        h4 = plot(ax,xx, 1 + aLin(k)*xx          ,'k--','LineWidth',1);

        xlim(ax,[0 10]); ylim(ax,[0 1]); grid(ax,'on')
        title(ax,sprintf('P %d',k));

        % store legend handles once
        if isempty(hScatter)
            hScatter = h1; hExpLine = h2; hRatLine = h3; hLinLine = h4;
        end
    else                                    % -------- excluded driver ----
        axis(ax,'off')
        text(0.5,0.5,{'excluded','(R^2 ≤ 0)'}, ...
             'Units','normalized','HorizontalAlignment','center', ...
             'FontSize',7,'Color',[.5 .5 .5],'Parent',ax);
        title(ax,sprintf('P %d',k));
    end
end

% ---------------- global legend / title / labels ------------------------
% --- global legend (no layout handle!) ----------------------------------
lgd = legend([hScatter hExpLine hRatLine hLinLine], ...
             {'samples','exp','inverse','linear'}, ...
             'Orientation','horizontal','NumColumns',4);

lgd.Layout.Tile = 'north';      % dock it above the 8 × 8 grid


sgtitle(tlKeep,'THW vs Effort – profiles with positive R²','FontWeight','bold');
xlabel(tlKeep,'Time Headway (s)');
ylabel(tlKeep,'Effort (normalized)');

% tighten spacing a touch
tlKeep.Padding    = 'loose';
tlKeep.TileSpacing= 'compact';

%% P value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------------------------------------------
% EXTRA:  bootstrap p-values for each fitted parameter
%         (add directly after the for-loop that fills hExp … rmseLin )
% ------------------------------------------------------------------------
% ================================================================
%  ADD THIS SECTION *AFTER* you have obtained hExp, hRat, aLin etc.
% ================================================================

B = 1000;                       % # permutations  (adjust if you wish)

% allocate
pFitExp = nan(nDrv,1);
pFitRat = nan(nDrv,1);
pFitLin = nan(nDrv,1);

% ---- open the biggest pool the profile allows ------------------
if isempty(gcp('nocreate'))
    parpool('local');           % uses NumWorkers in the profile
end

parfor k = 1:nDrv
    % -------------------------------------------------------------
    x = DriverCalc(k).THW(:);
    y = DriverCalc(k).Effort(:);
    n = numel(x);

    % observed SSEs ----------------------------------------------
    SSEexp_obs = sum( (y - expPred(hExp(k),x)).^2 );
    SSErat_obs = sum( (y - ratPred(hRat(k),x)).^2 );
    SSElin_obs = sum( (y - (1 + aLin(k)*x)).^2 );

    % permutation holders
    sseExp_perm = zeros(B,1);
    sseRat_perm = zeros(B,1);
    sseLin_perm = zeros(B,1);

    % one permuted fit uses the *same* fine grid you already made
    for b = 1:B
        yPerm = y(randperm(n));          % shuffle Effort

        % --- exponential ----------------------------------------
        SSE = arrayfun(@(h) sse(expPred(h,x),yPerm), gridFine);
        sseExp_perm(b) = min(SSE);

        % --- ratio ----------------------------------------------
        SSE = arrayfun(@(h) sse(ratPred(h,x),yPerm), gridFine);
        sseRat_perm(b) = min(SSE);

        % --- linear (closed form) -------------------------------
        aTmp = sum((yPerm-1).*x) / sum(x.^2);
        sseLin_perm(b) = sum( (yPerm - (1 + aTmp*x)).^2 );
    end

    % two-sided permutation p-values  ( +1 to avoid zero )
    pFitExp(k) = (1 + nnz(sseExp_perm <= SSEexp_obs)) / (B + 1);
    pFitRat(k) = (1 + nnz(sseRat_perm <= SSErat_obs)) / (B + 1);
    pFitLin(k) = (1 + nnz(sseLin_perm <= SSElin_obs)) / (B + 1);
end

% ---------- append to results table ------------------------------
FitTable.pFit_exp  = pFitExp;
FitTable.pFit_rat  = pFitRat;
FitTable.pFit_lin  = pFitLin;

% save updated table
save   effort_model_fits.mat  FitTable
writetable(FitTable,'effort_model_fits.csv')


%% ---------- add to results table & save again --------------------------
% ---------- append permutation p-values to results table --------------
FitTable.pFit_exp = pFitExp;
FitTable.pFit_rat = pFitRat;
FitTable.pFit_lin = pFitLin;


save   effort_model_fits.mat  FitTable
writetable(FitTable,'effort_model_fits.csv')
fprintf('Saved fit + p-values to effort_model_fits.(mat|csv)\n');

%% Residual plot

%%------------------------------------------------------------
%  1) residual vectors for every participant / every model
% ------------------------------------------------------------
nDrv = numel(DriverCalc);
resExp = cell(nDrv,1);     % cell arrays keep lengths separate
resRat = cell(nDrv,1);
resLin = cell(nDrv,1);

for k = 1:nDrv
    x = DriverCalc(k).THW(:);
    y = DriverCalc(k).Effort(:);

    resExp{k} = y - exp(-x./hExp(k));
    resRat{k} = y - hRat(k)./(hRat(k)+x);
    resLin{k} = y - (1 + aLin(k)*x);
end

% concatenate for global plots
allExp = vertcat(resExp{:});
allRat = vertcat(resRat{:});
allLin = vertcat(resLin{:});

%% ------------------------------------------------------------
% 2) helper for an 8×8 residual scatter grid
% ------------------------------------------------------------
plotResidualGrid(resExp, '\epsilon  =  y -  y_{exp}', cfg);
plotResidualGrid(resRat, '\epsilon  =  y -  y_{ratio}', cfg);
plotResidualGrid(resLin, '\epsilon  =  y -  y_{linear}', cfg);

%% ------------------------------------------------------------
% 3) histogram + QQ-plot for each model
% ------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 0.4]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact')

nexttile
histogram(allExp,40,'FaceColor',[0.2 0.6 1]); hold on
xline(mean(allExp),'k','LineWidth',1.5)
title('Residuals – exp model')

nexttile
histogram(allRat,40,'FaceColor',[1 0 0]); hold on
xline(mean(allRat),'k','LineWidth',1.5)
title('Residuals – ratio model')

nexttile
histogram(allLin,40,'FaceColor',[0 0 0]); hold on
xline(mean(allLin),'w','LineWidth',1.5)
title('Residuals – linear model')

figure('units','normalized','outerposition',[0 0 1 0.4]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact')
qqplotPanel(allExp,'exp');
qqplotPanel(allRat,'inverse');
qqplotPanel(allLin,'linear');

%% select those with positive R-Squaired
% --- 0) Select drivers whose *any* model has R² > 0 ----------------------
keepIdx = (FitTable.R2_exp > 0) | ...
          (FitTable.R2_rat > 0) | ...
          (FitTable.R2_lin > 0);

if ~any(keepIdx)
    warning('No drivers satisfy R² > 0 criterion. Nothing to plot.'); return
end

DriverCalcSel = DriverCalc(keepIdx);
DriverPlotSel = DriverPlot(keepIdx);      % for alpha (optional)
hExpSel = hExp(keepIdx);
hRatSel = hRat(keepIdx);
aLinSel = aLin(keepIdx);

nDrv = numel(DriverCalcSel);

% -------------------------------------------------------------------------
% --- 1)  gather residuals for the selected drivers -----------------------

keepMask = (FitTable.R2_exp > 0) | (FitTable.R2_rat > 0) | (FitTable.R2_lin > 0);
resExp = cell(64,1);
resRat = cell(64,1);
resLin = cell(64,1);

expPred = @(h,x) exp(-x./h);
ratPred = @(h,x) h ./ (h + x);

for k = 1:64
    if keepMask(k)                      % driver passed the R² > 0 filter
        x = DriverCalc(k).THW(:);
        y = DriverCalc(k).Effort(:);

        resExpFull{k} = y - exp(-x ./ hExp(k));
        resRatFull{k} = y - hRat(k) ./ (hRat(k) + x);
        resLinFull{k} = y - (1 + aLin(k) * x);
        % keepMask(k)==false → corresponding cell stays {}
    end
end

% -------------------------------------------------------------------------
% --- 2)  plot residual grids (function defined at end) -------------------
plotResidualGridFull(resExpFull, keepMask, '\epsilon = y - y_{exp}',   cfg);
plotResidualGridFull(resRatFull, keepMask, '\epsilon = y - y_{ratio}', cfg);
plotResidualGridFull(resLinFull, keepMask, '\epsilon = y - y_{linear}',cfg);

% -------- concatenate residuals from kept drivers ----------------------
allExp = vertcat(resExpFull{keepMask});
allRat = vertcat(resRatFull{keepMask});
allLin = vertcat(resLinFull{keepMask});

% ---------- 1) Histograms ----------------------------------------------
figure('units','normalized','outerposition',[0 0 1 0.4]);
tlH = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile(tlH)
histogram(allExp,40,'FaceColor',[0.2 0.6 1]); hold on
xline(mean(allExp),'k','LineWidth',1.5)
title('Residuals – exp model'); grid on

nexttile(tlH)
histogram(allRat,40,'FaceColor',[1 0 0]); hold on
xline(mean(allRat),'k','LineWidth',1.5)
title('Residuals – ratio model'); grid on

nexttile(tlH)
histogram(allLin,40,'FaceColor',[0.2 0.2 0.2]); hold on
xline(mean(allLin),'w','LineWidth',1.5)
title('Residuals – linear model'); grid on

sgtitle(tlH,'Histograms of residuals  (kept drivers)','FontWeight','bold')

% ---------- 2) QQ-plots -------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 0.4]);
tlQ = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

qqplotPanelnew(allExp,'exp', tlQ);
qqplotPanelnew(allRat,'inverse', tlQ);
qqplotPanelnew(allLin,'linear', tlQ);

sgtitle(tlQ,'QQ-plots of residuals  (kept drivers)','FontWeight','bold')

%% ---------------------------------------------------------------------
%  Box-plots of fitted parameters  (drivers selected by keepMask)
% ---------------------------------------------------------------------
mask = keepMask;            % 64-by-1 logical you already created

% vectors of the parameters you want to visualise
h0exp = FitTable.h0_exp(mask);
h0rat = FitTable.h0_rat(mask);
aLin  = FitTable.a_lin(mask);

figure('units','normalized','outerposition',[0 0 0.75 0.45]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ---- (1)  h0  of the exponential model  -------------------------------
nexttile(tl,1)
boxplot(h0exp,'Notch','off','Labels',{'h_exp '})
set(gca,'YScale','log')            % logarithmic is more informative
ylabel('h_{exp}  [s]')
% title(sprintf('exp model   (n = %d)',numel(h0exp)))
title(sprintf('exp model'))

% ---- (2)  h_inv  of the ratio model  -------------------------------------
nexttile(tl,2)
boxplot(h0rat,'Notch','off','Labels',{'h_inv '})
set(gca,'YScale','log')
ylabel('h_{inv}  [s]')
title(sprintf('inverse model'))

% ---- (3)  slope  k  of the linear model  ------------------------------
nexttile(tl,3)
boxplot(aLin ,'Notch','off','Labels',{'slope  k'})
yline(0,'k--');                    % zero reference
ylabel('k')
title(sprintf('linear model'))
% grid on

sgtitle(tl,'Fitted-parameter distributions ','FontWeight','bold');

%% save all figures
if ~exist('./curveFittingFigures','dir')
    mkdir('./curveFittingFigures');
end

allFigs = findall(0, 'Type','figure');
for i = 1:length(allFigs)
    figure(allFigs(i));
    filename = sprintf('./curveFittingFigures/Figure_%d.jpg', i);
    exportgraphics(allFigs(i), filename, 'Resolution',300);
end


%%

% ========================================================================
% GLOBAL fits  (all dots from all participants)
% ========================================================================

% ----- (1) exponential:  E = exp(-x / h0) -------------------------------
sseExpAll = arrayfun(@(h) sse(expPred(h,allTHW), allEffort), gridCoarse);
[~,idx]   = min(sseExpAll);  h0Start   = gridCoarse(idx);

gridFine  = linspace(0.75*h0Start, 1.25*h0Start, 60);
sseFine   = arrayfun(@(h) sse(expPred(h,allTHW), allEffort), gridFine);
[~,idx2]  = min(sseFine);    hExpAll   = gridFine(idx2);
[R2ExpAll, RMSEExpAll] = stats(expPred(hExpAll,allTHW), allEffort);

% ----- (2) ratio:  E = h0 / (h0 + x) ------------------------------------
sseRatAll = arrayfun(@(h) sse(ratPred(h,allTHW), allEffort), gridCoarse);
[~,idx]   = min(sseRatAll);  h0Start   = gridCoarse(idx);

gridFine  = linspace(0.75*h0Start, 1.25*h0Start, 60);
sseFine   = arrayfun(@(h) sse(ratPred(h,allTHW), allEffort), gridFine);
[~,idx2]  = min(sseFine);    hRatAll   = gridFine(idx2);
[R2RatAll, RMSERatAll] = stats(ratPred(hRatAll,allTHW), allEffort);

% ----- (3) linear:  E = 1 + a·x   (closed-form) -------------------------
aLinAll   = sum((allEffort-1).*allTHW) / sum(allTHW.^2);
[R2LinAll, RMSELinAll] = stats(1 + aLinAll*allTHW, allEffort);

fprintf('\nGLOBAL fits (all participants)\n');
fprintf('  exp:  h0 = %.3g  |  R² = %.3f  |  RMSE = %.3f\n', hExpAll, R2ExpAll, RMSEExpAll);
fprintf(' ratio: h0 = %.3g  |  R² = %.3f  |  RMSE = %.3f\n', hRatAll, R2RatAll, RMSERatAll);
fprintf(' linear:  a = %.3g |  R² = %.3f  |  RMSE = %.3f\n', aLinAll, R2LinAll, RMSELinAll);

% ---------- pooled scatter + global curves ------------------------------
figure('units','normalized','outerposition',[0 0 1 1]);

sh = randperm(numel(allTHW), min(6e4,numel(allTHW)));   % down-sample dots
alphaAll = fastAlpha(allTHW, allEffort);

scatter(allTHW(sh), allEffort(sh), 6, cfg.deepBlue, 'filled', ...
        'MarkerFaceAlpha','flat','MarkerEdgeAlpha','flat', ...
        'AlphaData',alphaAll(sh),'AlphaDataMapping','none'); hold on

xx = linspace(0,10,300);
plot(xx, expPred(hExpAll,xx) ,'b' ,'LineWidth',2);
plot(xx, ratPred(hRatAll,xx) ,'r' ,'LineWidth',2);
plot(xx, 1 + aLinAll*xx      ,'k--','LineWidth',2);

xlim([0 10]); ylim([0 1]); grid on
legend({'samples','exp','inverse','linear'},'Location','southwest');
title('THW vs Effort – global fits (all drivers)');
xlabel('Time Headway (s)');  ylabel('Effort (normalized)');




%% ======================================================================
%                              FUNCTIONS
% ======================================================================

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
        c    = cRaw(1:numel(x));                % align with trimmed len
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
    gap   = 0.05;                   % horizontal gap from right edge
    width = 0.015;                  % bar width
    bottom= 0.05; height = 0.90;    % vertical span

cb.Position = [1-gap-width, bottom, width, height];
end

function plotGroupMeanRibbon(Dplot, Dcalc, trim, cfg)
    % ----- density cloud (Dplot) ----------------------------------
    allX = cat(1,Dplot.THW);   allY = cat(1,Dplot.Effort);
    alpha = fastAlpha(allX, allY);

    % ----- stats on full data (Dcalc) -----------------------------
    edges   = trim.edgesTHW;                 centers = edges(1:end-1)+diff(edges)/2;
    nBins   = numel(centers);   nDrv = numel(Dcalc);
    meanMat = nan(nDrv,nBins);

    for k = 1:nDrv
        x = Dcalc(k).THW;  y = Dcalc(k).Effort;
        for i = 1:nBins
            m = x>=edges(i) & x<edges(i+1);  if any(m), meanMat(k,i)=mean(y(m)); end
        end
    end

    nPerBin = sum(~isnan(meanMat),1);
    gMean   = nanmean(meanMat,1);
    gSEM    = nanstd(meanMat,0,1) ./ sqrt(nPerBin);
    ciMult  = 1.96*strcmpi(trim.ciType,'ci95') + strcmpi(trim.ciType,'sem');
    lo      = gMean - ciMult*gSEM;    hi = gMean + ciMult*gSEM;
    valid   = nPerBin >= trim.minDrivers;

    % ----- draw ----------------------------------------------------
    figure('units','normalized','outerposition',[0 0 1 1]);
    sh = randperm(numel(allX), min(4e4,numel(allX)));   % down-sample dots
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

% ---------- utilities ---------------------------------------------------
% --- build safe edges -----------------------------------------


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
        y_n = (y - min(y)) ./ d;   % rescale to 0–1
    else
        y_n = 0.5 * ones(size(y)); % constant vector → mid-level
    end
end

% ------------------------------------------------------------------
%  fastAlpha  –  quick 2-D histogram → opacity in [0.05  1]
% ------------------------------------------------------------------
function a = fastAlpha(x,y)
    ok = ~(isnan(x) | isnan(y));
    x = x(ok);  y = y(ok);

    if numel(x) < 2               % one (or zero) point → constant alpha
        a = 0.5 * ones(size(x));
        return
    end

    nb = 50;                      % 50×50 grid
    % --- safe edges (avoid identical min/max crash) ----------------
    if range(x)==0
        ex = [x(1)-0.01, x(1)+0.01];
    else
        ex = linspace(min(x), max(x), nb+1);
    end
    if range(y)==0
        ey = [y(1)-0.01, y(1)+0.01];
    else
        ey = linspace(min(y), max(y), nb+1);
    end
    % --- histogram density -----------------------------------------
    [N,~,~,ix,iy] = histcounts2(x,y,ex,ey);
    lin = sub2ind(size(N), ix, iy);
    a   = rescale(N(lin), 0.05, 1.0);   % map counts → alpha
end


% ========================================================================
%                      helper functions used above
% ========================================================================
% function [R2,RMSE]=stats(~,~); end %#ok<DEFNU>  % (placeholder for lint)



function plotResidualGrid(resCell, yLabelTxt, cfg)
    % resCell : 64×1 cell, each contains residual vector for that driver
    figure('units','normalized','outerposition',[0 0 1 1]);
    for k = 1:numel(resCell)
        subplot(8,8,k)
        n = numel(resCell{k});
        scatter(1:n, resCell{k}, 6, cfg.deepBlue, 'filled',...
                'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6);
        yline(0,'k--');
        ylim([-0.6 0.6]);               % adjust as needed
        xlim([0 n+1]); grid on
        title(sprintf('P %d',k));
    end
    axHost = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5, 0.98, ['Residuals grid  (' yLabelTxt ')'], ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontWeight','bold','Parent',axHost);
    text(0.5, -0.04, 'Sample index','Units','normalized', ...
         'HorizontalAlignment','center','Parent',axHost);
    text(-0.04, 0.5, yLabelTxt,'Units','normalized', ...
         'HorizontalAlignment','center','Rotation',90,'Parent',axHost);
end

function qqplotPanel(r, modelName)
    nexttile
    qqplot(r); grid on
    title(['QQ-plot  (' modelName ')'])
end


% ========================================================================
function plotResidualGridFull(resCell, keepMask, yLabelTxt, cfg)
% resCell  : 64×1 cell array of residual vectors or {} for excluded drv
% keepMask : logical 64×1, true if driver kept
% yLabelTxt: e.g. '\epsilon = y - y_{exp}'
%
% Uses the original 8 × 8 grid so driver IDs match your earlier plots.

    figure('units','normalized','outerposition',[0 0 1 1]);

    nRows = 8; nCols = 8;

    for k = 1:64
        subplot(nRows, nCols, k)

        if keepMask(k)                     % ---------- kept driver -------
            r = resCell{k};
            scatter(1:numel(r), r, 6, cfg.deepBlue, 'filled', ...
                    'MarkerFaceAlpha',0.6,'MarkerEdgeAlpha',0.6); hold on
            yline(0,'k--');
            ylim([-0.6 0.6]); xlim([0 numel(r)+1]); grid on

        else                               % ---------- excluded driver ---
            axis off
            text(0.5,0.5,{'excluded','(R^2 \le 0)'}, ...
                 'HorizontalAlignment','center', ...
                 'FontSize',7,'Color',[.4 .4 .4]);
        end

        title(sprintf('P %d', k));
    end

    % -------- global labels via invisible host axes -----------------
    axHost = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5, 0.98, ['Residuals grid  (' yLabelTxt ')'], ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontWeight','bold','Parent',axHost);
    text(0.5, -0.04, 'Sample index', ...
         'Units','normalized','HorizontalAlignment','center', ...
         'Parent',axHost);
    text(-0.04, 0.5, yLabelTxt, ...
         'Units','normalized','HorizontalAlignment','center', ...
         'Rotation',90,'Parent',axHost);
end


function qqplotPanelnew(resVec, labelStr, parentTL)
    nexttile(parentTL)
    qqplot(resVec); grid on
    title(['QQ-plot (' labelStr ')'])
end

function [logL, AIC, BIC] = normalLogLike(y, yHat, k)
    n      = numel(y);
    SSE    = sum((y - yHat).^2);
    sigma2 = SSE / n;
    logL   = -0.5*n * ( log(2*pi*sigma2) + 1 );   % Gaussian log-likelihood
    AIC    = 2*k - 2*logL;
    BIC    = k*log(n) - 2*logL;
end
