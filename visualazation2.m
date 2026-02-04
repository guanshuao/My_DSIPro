%% 将Stack.mat文件中的数据可视化
clc;clear;close all;

% 检查Stack变量是否存在，否则尝试从 .mat 文件加载
if ~exist('Stack', 'var')
    if exist('Stack.mat', 'file')
        disp('Loading Stack.mat...');
        load('Stack.mat');
    else
        error('Stack variable not in workspace and Stack.mat not found.');
    end
end

% - Stack结构体包含多个字段，具体如下：
% - TrueCoh: 真实的相干矩阵
% - BoxCoh_XX: 使用窗口大小为XX的盒式相干估计，例如BoxCoh_3表示窗口大小为3x3估计出的相干矩阵
% - SHPNum: 每个位置的统计同质性像素（Statistically Homogeneous Pixels）数矩阵，例如某元素值为50表示该位置周围有50个同质性像素。选取SHP时使用的窗口大小为11x11，共121个像素
% - AdpCph: 自适应相干估计矩阵。即估计相干性时只使用统计同质性像素的数值
% - CV：由RCS计算得到的变异系数矩阵，衡量地表的均匀程度（改进后仅保留一个CV字段）
% - Best_Est: 同一个位置中，所有BoxCoh_XX估计中与TrueCoh最接近的估计矩阵
% - Best_Win: 对应Best_Est的窗口大小矩阵，表示每个位置使用的最佳窗口大小，如某一位置Best_Win值为5表示该位置使用5x5窗口估计出的相干矩阵最接近真实值TrueCoh
% - ModifiedCoh: 修正后的相干矩阵
% - Calwin: 用于选取SHP的窗口大小矩阵，如 [11 11]
Calwin = Stack.Calwin;
PixelNum = Calwin(1) * Calwin(2);  % 计算窗口内像素总数

%% 图一：横轴为CV，纵轴比较各种估计方法的均值，包括BoxCoh_11、AdpCph、ModifiedCoh、TrueCoh
% CV是连续变量，使用滑动窗口方法计算均值

% 定义CV的采样点和窗口
cv_all = Stack.CV(:);
cv_all = cv_all(~isnan(cv_all));
cv_min = min(cv_all);
cv_max = max(cv_all);
cv_sample_points = linspace(cv_min, cv_max, 100);  % 100个采样点
cv_window_size = (cv_max - cv_min) / 50;  % 窗口宽度

% 初始化存储数组
mean_TrueCoh_cv = zeros(size(cv_sample_points));
mean_BoxCoh11_cv = zeros(size(cv_sample_points));
mean_AdpCoh_cv = zeros(size(cv_sample_points));
mean_ModifiedCoh_cv = zeros(size(cv_sample_points));
mean_BestEst_cv = zeros(size(cv_sample_points));

% 对每个CV采样点计算均值
for i = 1:length(cv_sample_points)
    cv_val = cv_sample_points(i);
    % 找到CV在当前窗口内的所有位置
    idx = (Stack.CV >= cv_val - cv_window_size/2) & ...
          (Stack.CV <= cv_val + cv_window_size/2);
    
    if sum(idx(:)) > 0  % 如果存在该窗口内的点
        mean_TrueCoh_cv(i) = mean(Stack.TrueCoh(idx), 'omitnan');
        mean_BoxCoh11_cv(i) = mean(Stack.BoxCoh_11(idx), 'omitnan');
        mean_AdpCoh_cv(i) = mean(Stack.AdpCoh(idx), 'omitnan');
        mean_ModifiedCoh_cv(i) = mean(Stack.ModifiedCoh(idx), 'omitnan');
        mean_BestEst_cv(i) = mean(Stack.Best_Est(idx), 'omitnan');
    else
        mean_TrueCoh_cv(i) = NaN;
        mean_BoxCoh11_cv(i) = NaN;
        mean_AdpCoh_cv(i) = NaN;
        mean_ModifiedCoh_cv(i) = NaN;
        mean_BestEst_cv(i) = NaN;
    end
end

% 绘制对比图
figure;
plot(cv_sample_points, mean_TrueCoh_cv, 'k-', 'LineWidth', 2, 'DisplayName', 'TrueCoh'); hold on;
plot(cv_sample_points, mean_BoxCoh11_cv, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11');
plot(cv_sample_points, mean_AdpCoh_cv, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh');
plot(cv_sample_points, mean_ModifiedCoh_cv, 'g-', 'LineWidth', 1.5, 'DisplayName', 'ModifiedCoh');
plot(cv_sample_points, mean_BestEst_cv, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Best\_Est');
hold off;

xlabel('CV (变异系数)');
ylabel('相干性均值');
ylim([0 1]);
title('不同CV下各估计方法的相干性均值对比');
legend('Location', 'best');
grid on;

%% 图二：横轴为CV，纵轴比较各种估计方法之比。1: Box/Adp, 2: Box/Mod, 3: Adp/Mod
% 初始化存储比值的数组
ratio_BoxAdp_cv = zeros(size(cv_sample_points));      % BoxCoh_11 / AdpCoh
ratio_BoxMod_cv = zeros(size(cv_sample_points));      % BoxCoh_11 / ModifiedCoh
ratio_AdpMod_cv = zeros(size(cv_sample_points));      % AdpCoh / ModifiedCoh

% 对每个CV采样点计算比值
for i = 1:length(cv_sample_points)
    cv_val = cv_sample_points(i);
    % 找到CV在当前窗口内的所有位置
    idx = (Stack.CV >= cv_val - cv_window_size/2) & ...
          (Stack.CV <= cv_val + cv_window_size/2);
    
    if sum(idx(:)) > 0  % 如果存在该窗口内的点
        % 计算均值
        mean_box = mean(Stack.BoxCoh_11(idx), 'omitnan');
        mean_adp = mean(Stack.AdpCoh(idx), 'omitnan');
        mean_mod = mean(Stack.ModifiedCoh(idx), 'omitnan');
        
        % 计算比值（避免除以0）
        if mean_adp > 0
            ratio_BoxAdp_cv(i) = mean_box / mean_adp;
        else
            ratio_BoxAdp_cv(i) = NaN;
        end
        
        if mean_mod > 0
            ratio_BoxMod_cv(i) = mean_box / mean_mod;
        else
            ratio_BoxMod_cv(i) = NaN;
        end
        
        if mean_mod > 0
            ratio_AdpMod_cv(i) = mean_adp / mean_mod;
        else
            ratio_AdpMod_cv(i) = NaN;
        end
    else
        ratio_BoxAdp_cv(i) = NaN;
        ratio_BoxMod_cv(i) = NaN;
        ratio_AdpMod_cv(i) = NaN;
    end
end

% 绘制比值对比图
figure;
plot(cv_sample_points, ones(size(cv_sample_points)), 'k--', 'LineWidth', 1, 'DisplayName', '参考线 (y=1)'); hold on;
plot(cv_sample_points, ratio_BoxAdp_cv, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 / AdpCoh');
plot(cv_sample_points, ratio_BoxMod_cv, 'r-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 / ModifiedCoh');
plot(cv_sample_points, ratio_AdpMod_cv, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh / ModifiedCoh');
hold off;

xlabel('CV (变异系数)');
ylabel('相干性估计比值');
ylim([0 4]);
title('不同CV下各估计方法之间的比值对比');
legend('Location', 'best');
grid on;

%% 图三：横轴为CV，纵轴为相干性估计值，绘制箱线图，包括BoxCoh_11、AdpCoh、ModifiedCoh、TrueCoh、BestEst
% CV是连续变量，需要先划分区间

% 定义CV区间（按分位数划分，每个区间包含12.5%的像元）
num_cv_intervals = 8;
percentiles = linspace(0, 100, num_cv_intervals + 1);
cv_edges = prctile(cv_all, percentiles);
cv_intervals = [cv_edges(1:end-1)', cv_edges(2:end)'];

% 为每个估计方法准备数据
data_TrueCoh_cv = cell(num_cv_intervals, 1);
data_BoxCoh11_cv = cell(num_cv_intervals, 1);
data_AdpCoh_cv = cell(num_cv_intervals, 1);
data_BestEst_cv = cell(num_cv_intervals, 1);

% 对每个区间收集数据
for i = 1:num_cv_intervals
    interval_start = cv_intervals(i, 1);
    interval_end = cv_intervals(i, 2);
    
    % 找到CV在当前区间内的所有位置
    idx = (Stack.CV >= interval_start) & (Stack.CV <= interval_end);
    
    if sum(idx(:)) > 0  % 如果存在该区间内的点
        data_TrueCoh_cv{i} = Stack.TrueCoh(idx);
        data_BoxCoh11_cv{i} = Stack.BoxCoh_11(idx);
        data_AdpCoh_cv{i} = Stack.AdpCoh(idx);
        data_BestEst_cv{i} = Stack.Best_Est(idx);
    else
        data_TrueCoh_cv{i} = [];
        data_BoxCoh11_cv{i} = [];
        data_AdpCoh_cv{i} = [];
        data_BestEst_cv{i} = [];
    end
end

% 准备箱线图的标签
cv_interval_labels = arrayfun(@(i) sprintf('%.2f-%.2f', cv_intervals(i,1), cv_intervals(i,2)), ...
                          1:num_cv_intervals, 'UniformOutput', false);

% 绘制综合箱线图
figure;
hold on;

% 定义颜色
color_TrueCoh = [0 0 0];  % 黑色
color_BoxCoh11 = [0 0.447 0.741];  % 蓝色
color_AdpCoh = [0.850 0.325 0.098];  % 红色
color_BestEst = [0.929 0.694 0.125];  % 黄色

% 定义每组箱线图的间距
group_width = 4;  % 每个CV区间占据的x轴宽度
box_width = 1;    % 每个箱线图的宽度

% 为每个区间绘制四个箱线图
for i = 1:num_cv_intervals
    % 计算当前区间四个箱线图的中心位置
    center_pos = i * group_width;
    positions = center_pos + [-1.5*box_width, -0.5*box_width, 0.5*box_width, 1.5*box_width];
    
    % 绘制TrueCoh箱线图
    if ~isempty(data_TrueCoh_cv{i})
        h1 = boxplot(data_TrueCoh_cv{i}, 'Positions', positions(1), 'Widths', box_width*0.8, ...
                     'Colors', color_TrueCoh, 'Symbol', '');
        set(h1, 'LineWidth', 1.5);
    end
    
    % 绘制BoxCoh_11箱线图
    if ~isempty(data_BoxCoh11_cv{i})
        h2 = boxplot(data_BoxCoh11_cv{i}, 'Positions', positions(2), 'Widths', box_width*0.8, ...
                     'Colors', color_BoxCoh11, 'Symbol', '');
        set(h2, 'LineWidth', 1.5);
    end
    
    % 绘制AdpCoh箱线图
    if ~isempty(data_AdpCoh_cv{i})
        h3 = boxplot(data_AdpCoh_cv{i}, 'Positions', positions(3), 'Widths', box_width*0.8, ...
                     'Colors', color_AdpCoh, 'Symbol', '');
        set(h3, 'LineWidth', 1.5);
    end
    
    % 绘制BestEst箱线图
    if ~isempty(data_BestEst_cv{i})
        h4 = boxplot(data_BestEst_cv{i}, 'Positions', positions(4), 'Widths', box_width*0.8, ...
                     'Colors', color_BestEst, 'Symbol', '');
        set(h4, 'LineWidth', 1.5);
    end
end

% 设置x轴刻度和标签
x_tick_positions = (1:num_cv_intervals) * group_width;
set(gca, 'XTick', x_tick_positions);
set(gca, 'XTickLabel', cv_interval_labels);
xtickangle(45);

% 设置坐标轴标签和标题
xlabel('CV区间');
ylabel('相干性估计值');
title('不同CV区间下各估计方法的分布对比（综合箱线图）');
grid on;

% 添加图例
% 创建虚拟的线条用于图例
legend_handles = [
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_TrueCoh, 'MarkerEdgeColor', color_TrueCoh, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_BoxCoh11, 'MarkerEdgeColor', color_BoxCoh11, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_AdpCoh, 'MarkerEdgeColor', color_AdpCoh, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_BestEst, 'MarkerEdgeColor', color_BestEst, 'MarkerSize', 10);
];
legend(legend_handles, {'TrueCoh', 'BoxCoh\_11', 'AdpCoh', 'Best\_Est'}, 'Location', 'best');
hold off;

%% 图四：SHPNum与CV的联合分布及其对估计误差的影响（二维热力图）
% 设计思路：
%   - 横轴：SHPNum（统计同质性像素数）
%   - 纵轴：CV（变异系数）
%   - 颜色：|BoxCoh_11 - TrueCoh|（估计误差绝对值）

figure;

% 定义SHPNum和CV的分箱数
num_shp_bins = 30;  % SHPNum的分箱数
num_cv_bins = 30;   % CV的分箱数

% 获取SHPNum的范围
shp_all = Stack.SHPNum(:);
shp_all = shp_all(~isnan(shp_all));
shp_min = double(min(shp_all));
shp_max = double(max(shp_all));

% 创建分箱边界
shp_edges = linspace(shp_min, shp_max, num_shp_bins + 1);
cv_edges_heatmap = linspace(cv_min, cv_max, num_cv_bins + 1);

% 计算每个bin的中心点
shp_centers = (shp_edges(1:end-1) + shp_edges(2:end)) / 2;
cv_centers = (cv_edges_heatmap(1:end-1) + cv_edges_heatmap(2:end)) / 2;

% 初始化热力图矩阵（存储每个bin的平均误差）
error_matrix = NaN(num_cv_bins, num_shp_bins);
count_matrix = zeros(num_cv_bins, num_shp_bins);

% 计算每个bin的平均估计误差
for i = 1:num_shp_bins
    for j = 1:num_cv_bins
        % 找到同时满足SHPNum和CV条件的位置
        idx = (Stack.SHPNum >= shp_edges(i)) & (Stack.SHPNum < shp_edges(i+1)) & ...
              (Stack.CV >= cv_edges_heatmap(j)) & (Stack.CV < cv_edges_heatmap(j+1));
        
        count = sum(idx(:));
        count_matrix(j, i) = count;
        
        if count > 0
            % 计算该bin内的平均误差绝对值 |BoxCoh_11 - TrueCoh|
            error_matrix(j, i) = mean(abs(Stack.BoxCoh_11(idx) - Stack.TrueCoh(idx)), 'omitnan');
        end
    end
end

% 绘制热力图
imagesc(shp_centers, cv_centers, error_matrix);
set(gca, 'YDir', 'normal');  % 使Y轴从下到上递增
colormap(jet);
c = colorbar;
c.Label.String = '|BoxCoh_{11} - TrueCoh|';
c.Label.FontSize = 11;

% 设置坐标轴
xlabel('SHPNum (统计同质性像素数)');
ylabel('CV (变异系数)');
title('SHPNum与CV联合分布下的BoxCoh估计误差绝对值热力图');
grid off;

% 添加等值线以增强可读性
hold on;
[X, Y] = meshgrid(shp_centers, cv_centers);
% 只在有足够数据的区域绘制等值线
valid_mask = count_matrix >= 5;  % 至少5个样本点
error_matrix_contour = error_matrix;
error_matrix_contour(~valid_mask) = NaN;
contour(X, Y, error_matrix_contour, 8, 'k', 'LineWidth', 0.5);
hold off;

%% 图五：SHPNum与CV的联合分布下的最佳窗口大小热力图
% 设计思路：
%   - 横轴：SHPNum（统计同质性像素数）
%   - 纵轴：CV（变异系数）
%   - 颜色：Best_Win（最佳窗口大小）

figure;

% 使用与图四相同的分箱设置
% num_shp_bins = 30;  % 已在图四中定义
% num_cv_bins = 30;   % 已在图四中定义

% 初始化热力图矩阵（存储每个bin的最佳窗口平均值）
best_win_matrix = NaN(num_cv_bins, num_shp_bins);
count_matrix_win = zeros(num_cv_bins, num_shp_bins);

% 计算每个bin内最佳窗口大小的平均值
for i = 1:num_shp_bins
    for j = 1:num_cv_bins
        % 找到同时满足SHPNum和CV条件的位置
        idx = (Stack.SHPNum >= shp_edges(i)) & (Stack.SHPNum < shp_edges(i+1)) & ...
              (Stack.CV >= cv_edges_heatmap(j)) & (Stack.CV < cv_edges_heatmap(j+1));
        
        count = sum(idx(:));
        count_matrix_win(j, i) = count;
        
        if count > 0
            % 获取该bin内的Best_Win值（排除0和NaN）
            best_wins = Stack.Best_Win(idx);
            best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
            
            if ~isempty(best_wins)
                % 使用平均值表示该bin的最佳窗口大小
                best_win_matrix(j, i) = mean(best_wins);
            end
        end
    end
end

% 绘制热力图
imagesc(shp_centers, cv_centers, best_win_matrix);
set(gca, 'YDir', 'normal');  % 使Y轴从下到上递增

% 使用连续颜色映射（因为平均值是连续值）
colormap(jet);

% 设置colorbar
c = colorbar;
c.Label.String = '最佳窗口大小平均值 (Best\_Win)';
c.Label.FontSize = 11;

% 设置坐标轴
xlabel('SHPNum (统计同质性像素数)');
ylabel('CV (变异系数)');
title('SHPNum与CV联合分布下的最佳窗口大小平均值热力图');
grid off;

% 添加等值线以增强可读性
hold on;
% 只在有足够数据的区域绘制等值线
valid_mask_win = count_matrix_win >= 5;  % 至少5个样本点
best_win_contour = best_win_matrix;
best_win_contour(~valid_mask_win) = NaN;
contour(X, Y, best_win_contour, 8, 'k', 'LineWidth', 0.5);
hold off;

