% 取出第二列 r（若无变量名，用 T{:,2}）
x = T.r;              % 等价：x = T{:, 'r'};  或  x = T{:,2};

% 去掉 NaN
x = x(:);
x = x(~isnan(x));

% 1) 均值
mu = mean(x);
mid = median(x);

% 2) IQR：四分位数区间和宽度
q = prctile(x,[25 75]);   % q(1)=Q1, q(2)=Q3
IQR_width = q(2) - q(1);  % IQR 宽度

% 3) 小于 0 的百分比
pct_lt0 = 100 * mean(x < 0);

% 打印结果
fprintf('mean = %.6f\n', mu);
fprintf('median = %.6f\n', mid);
fprintf('IQR range = [%.6f, %.6f] (width = %.6f)\n', q(1), q(2), IQR_width);
fprintf('%%<0 = %.2f%%  (n=%d / %d)\n', pct_lt0, sum(x<0), numel(x));
