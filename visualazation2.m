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

%% 图一：横轴为CV等样本量分级（level1-level100），纵轴比较各种估计方法的均值
% 将CV从低到高排序，按等样本量划分为100个level，每个level内样本数相同

num_levels = 100;  % 分级数量

% 获取有效CV值及其线性索引
cv_all = Stack.CV(:);
valid_mask = ~isnan(cv_all);
valid_indices = find(valid_mask);  % 有效像元的线性索引
cv_valid = cv_all(valid_mask);

% 按CV从低到高排序，获取排序后的索引
[~, sort_order] = sort(cv_valid);
sorted_indices = valid_indices(sort_order);  % 排序后像元的线性索引

% 将排序后的像元等分为num_levels组
num_valid = length(sorted_indices);
group_size = floor(num_valid / num_levels);
level_labels = 1:num_levels;

% 初始化存储数组
mean_TrueCoh_lv = zeros(1, num_levels);
mean_BoxCoh11_lv = zeros(1, num_levels);
mean_AdpCoh_lv = zeros(1, num_levels);
mean_ModifiedCoh_lv = zeros(1, num_levels);
mean_BestEst_lv = zeros(1, num_levels);

% 对每个level计算均值
for i = 1:num_levels
    if i < num_levels
        idx = sorted_indices((i-1)*group_size+1 : i*group_size);
    else
        idx = sorted_indices((i-1)*group_size+1 : end);  % 最后一组包含剩余像元
    end
    
    mean_TrueCoh_lv(i) = mean(Stack.TrueCoh(idx), 'omitnan');
    mean_BoxCoh11_lv(i) = mean(Stack.BoxCoh_11(idx), 'omitnan');
    mean_AdpCoh_lv(i) = mean(Stack.AdpCoh(idx), 'omitnan');
    mean_ModifiedCoh_lv(i) = mean(Stack.ModifiedCoh(idx), 'omitnan');
    mean_BestEst_lv(i) = mean(Stack.Best_Est(idx), 'omitnan');
end

% 绘制对比图
figure;
plot(level_labels, mean_TrueCoh_lv, 'k-', 'LineWidth', 2, 'DisplayName', 'TrueCoh'); hold on;
plot(level_labels, mean_BoxCoh11_lv, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11');
plot(level_labels, mean_AdpCoh_lv, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh');
plot(level_labels, mean_ModifiedCoh_lv, 'g-', 'LineWidth', 1.5, 'DisplayName', 'ModifiedCoh');
plot(level_labels, mean_BestEst_lv, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Best\_Est');
hold off;

xlabel('CV Level (低 \rightarrow 高)');
ylabel('相干性均值');
xlim([1 num_levels]);
ylim([0 1]);
title('不同CV等级下各估计方法的相干性均值对比');
legend('Location', 'best');
grid on;

%% 图二：横轴为CV等样本量分级（level1-level100），纵轴比较各种估计方法之比
% 复用图一中已计算的level分组

% 初始化存储比值的数组
ratio_BoxAdp_lv = zeros(1, num_levels);      % BoxCoh_11 / AdpCoh
ratio_BoxMod_lv = zeros(1, num_levels);      % BoxCoh_11 / ModifiedCoh
ratio_AdpMod_lv = zeros(1, num_levels);      % AdpCoh / ModifiedCoh

% 对每个level计算比值
for i = 1:num_levels
    if i < num_levels
        idx = sorted_indices((i-1)*group_size+1 : i*group_size);
    else
        idx = sorted_indices((i-1)*group_size+1 : end);
    end
    
    % 计算均值
    mean_box = mean(Stack.BoxCoh_11(idx), 'omitnan');
    mean_adp = mean(Stack.AdpCoh(idx), 'omitnan');
    mean_mod = mean(Stack.ModifiedCoh(idx), 'omitnan');
    
    % 计算比值（避免除以0）
    if mean_adp > 0
        ratio_BoxAdp_lv(i) = mean_box / mean_adp;
    else
        ratio_BoxAdp_lv(i) = NaN;
    end
    
    if mean_mod > 0
        ratio_BoxMod_lv(i) = mean_box / mean_mod;
    else
        ratio_BoxMod_lv(i) = NaN;
    end
    
    if mean_mod > 0
        ratio_AdpMod_lv(i) = mean_adp / mean_mod;
    else
        ratio_AdpMod_lv(i) = NaN;
    end
end

% 绘制比值对比图
figure;
plot(level_labels, ones(1, num_levels), 'k--', 'LineWidth', 1, 'DisplayName', '参考线 (y=1)'); hold on;
plot(level_labels, ratio_BoxAdp_lv, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 / AdpCoh');
plot(level_labels, ratio_BoxMod_lv, 'r-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 / ModifiedCoh');
plot(level_labels, ratio_AdpMod_lv, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh / ModifiedCoh');
hold off;

xlabel('CV Level (低 \rightarrow 高)');
ylabel('相干性估计比值');
xlim([1 num_levels]);
ylim([0 4]);
title('不同CV等级下各估计方法之间的比值对比');
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
%   - 横轴：SHPNum等样本量分级（121个level）
%   - 纵轴：CV等样本量分级（100个level）
%   - 颜色：|BoxCoh_11 - TrueCoh|（估计误差绝对值）
%   - 共 100 x 121 = 12100 个格子

figure;

% 定义分箱数
num_shp_bins = 121;  % SHPNum的分级数
num_cv_bins = 100;   % CV的分级数

% --- CV: 按分位数划分为等样本量的100个区间 ---
cv_percentiles = linspace(0, 100, num_cv_bins + 1);
cv_edges_heatmap = prctile(cv_valid, cv_percentiles);
% 修正可能存在的重复边界（当大量像元CV值相同时）
cv_edges_heatmap = unique(cv_edges_heatmap);
num_cv_bins_actual = length(cv_edges_heatmap) - 1;

% --- SHPNum: 按分位数划分为等样本量的121个区间 ---
shp_all = Stack.SHPNum(:);
shp_valid_mask = ~isnan(shp_all);
shp_valid = double(shp_all(shp_valid_mask));
shp_percentiles = linspace(0, 100, num_shp_bins + 1);
shp_edges = prctile(shp_valid, shp_percentiles);
shp_edges = unique(shp_edges);
num_shp_bins_actual = length(shp_edges) - 1;

% 使用level序号作为坐标轴
shp_levels = 1:num_shp_bins_actual;
cv_levels = 1:num_cv_bins_actual;

% 初始化热力图矩阵（存储每个bin的平均误差）
error_matrix = NaN(num_cv_bins_actual, num_shp_bins_actual);
count_matrix = zeros(num_cv_bins_actual, num_shp_bins_actual);

% 计算每个bin的平均估计误差
parfor i = 1:num_shp_bins_actual
    for j = 1:num_cv_bins_actual
        % 找到同时满足SHPNum和CV条件的位置
        if i < num_shp_bins_actual
            shp_cond = (Stack.SHPNum >= shp_edges(i)) & (Stack.SHPNum < shp_edges(i+1));
        else
            shp_cond = (Stack.SHPNum >= shp_edges(i)) & (Stack.SHPNum <= shp_edges(i+1));
        end
        if j < num_cv_bins_actual
            cv_cond = (Stack.CV >= cv_edges_heatmap(j)) & (Stack.CV < cv_edges_heatmap(j+1));
        else
            cv_cond = (Stack.CV >= cv_edges_heatmap(j)) & (Stack.CV <= cv_edges_heatmap(j+1));
        end
        idx = shp_cond & cv_cond;
        
        count = sum(idx(:));
        count_matrix(j, i) = count;
        
        if count > 0
            % 计算该bin内的平均误差绝对值 |BoxCoh_11 - TrueCoh|
            error_matrix(j, i) = mean(abs(Stack.BoxCoh_11(idx) - Stack.TrueCoh(idx)), 'omitnan');
        end
    end
end

% 绘制热力图
imagesc(shp_levels, cv_levels, error_matrix);
set(gca, 'YDir', 'normal');  % 使Y轴从下到上递增
colormap(jet);
c = colorbar;
c.Label.String = '|BoxCoh_{11} - TrueCoh|';
c.Label.FontSize = 11;

% 设置坐标轴
xlabel('SHPNum Level (低 \rightarrow 高)');
ylabel('CV Level (低 \rightarrow 高)');
title('SHPNum与CV等样本量分级下的BoxCoh估计误差绝对值热力图');
grid off;

% 添加等值线以增强可读性
hold on;
[X, Y] = meshgrid(shp_levels, cv_levels);
% 只在有足够数据的区域绘制等值线
valid_mask = count_matrix >= 10;  % 至少3个样本点
error_matrix_contour = error_matrix;
error_matrix_contour(~valid_mask) = NaN;
contour(X, Y, error_matrix_contour, 8, 'k', 'LineWidth', 0.5);
hold off;

%% 图五：SHPNum与CV的联合分布下的最佳窗口大小热力图
% 设计思路：
%   - 横轴：SHPNum等样本量分级（121个level）
%   - 纵轴：CV等样本量分级（100个level）
%   - 颜色：Best_Win（最佳窗口大小）
%   - 共 100 x 121 = 12100 个格子

figure;

% 复用图四中已计算的等样本量分箱边界

% 初始化热力图矩阵（存储每个bin的最佳窗口平均值）
best_win_matrix = NaN(num_cv_bins_actual, num_shp_bins_actual);
count_matrix_win = zeros(num_cv_bins_actual, num_shp_bins_actual);

% 计算每个bin内最佳窗口大小的平均值
parfor i = 1:num_shp_bins_actual
    for j = 1:num_cv_bins_actual
        % 找到同时满足SHPNum和CV条件的位置
        if i < num_shp_bins_actual
            shp_cond = (Stack.SHPNum >= shp_edges(i)) & (Stack.SHPNum < shp_edges(i+1));
        else
            shp_cond = (Stack.SHPNum >= shp_edges(i)) & (Stack.SHPNum <= shp_edges(i+1));
        end
        if j < num_cv_bins_actual
            cv_cond = (Stack.CV >= cv_edges_heatmap(j)) & (Stack.CV < cv_edges_heatmap(j+1));
        else
            cv_cond = (Stack.CV >= cv_edges_heatmap(j)) & (Stack.CV <= cv_edges_heatmap(j+1));
        end
        idx = shp_cond & cv_cond;
        
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
imagesc(shp_levels, cv_levels, best_win_matrix);
set(gca, 'YDir', 'normal');  % 使Y轴从下到上递增

% 使用连续颜色映射（因为平均值是连续值）
colormap(jet);

% 设置colorbar
c = colorbar;
c.Label.String = '最佳窗口大小平均值 (Best\_Win)';
c.Label.FontSize = 11;

% 设置坐标轴
xlabel('SHPNum Level (低 \rightarrow 高)');
ylabel('CV Level (低 \rightarrow 高)');
title('SHPNum与CV等样本量分级下的最佳窗口大小平均值热力图');
grid off;

% 添加等值线以增强可读性
hold on;
[X5, Y5] = meshgrid(shp_levels, cv_levels);
% 只在有足够数据的区域绘制等值线
valid_mask_win = count_matrix_win >= 10;  % 至少3个样本点
best_win_contour = best_win_matrix;
best_win_contour(~valid_mask_win) = NaN;
contour(X5, Y5, best_win_contour, 8, 'k', 'LineWidth', 0.5);
hold off;

