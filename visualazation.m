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
% - RCS：雷达散射截面（Radar Cross Section, RCS）矩阵，由强度图像经过滤波得到
% - CV：由RCS矩阵计算得到的变异系数矩阵，衡量地表的均匀程度（改进后仅保留一个CV字段）
% - Best_Est: 同一个位置中，所有BoxCoh_XX估计中与TrueCoh最接近的估计矩阵
% - Best_Win: 对应Best_Est的窗口大小矩阵，表示每个位置使用的最佳窗口大小，如某一位置Best_Win值为5表示该位置使用5x5窗口估计出的相干矩阵最接近真实值TrueCoh
% - ModifiedCoh: 修正后的相干矩阵
% - Calwin: 用于选取SHP的窗口大小矩阵，如 [11 11]
Calwin = Stack.Calwin;
PixelNum = Calwin(1) * Calwin(2);  % 计算窗口内像素总数

%% 图一：横轴为SHPNum，纵轴比较各种估计方法的均值，包括BoxCoh_11、AdpCph、ModifiedCoh、TrueCoh
% 如SHPNum=10，则寻找所有SHPNum=10的点，计算这些点在TrueCoh、AdpCph、Best_Est上的均值并绘制对比图，以此类推，将SHPNum全部的均值都计算出来并用平滑的曲线绘制出来

% 初始化存储数组
shp_values = 0:PixelNum;  % SHPNum的取值范围
mean_TrueCoh = zeros(size(shp_values));
mean_BoxCoh11 = zeros(size(shp_values));
mean_AdpCoh = zeros(size(shp_values));
mean_ModifiedCoh = zeros(size(shp_values));
mean_BestEst = zeros(size(shp_values));

% 对每个SHPNum值计算均值
for i = 1:length(shp_values)
    shp_val = shp_values(i);
    % 找到SHPNum等于当前值的所有位置
    idx = (Stack.SHPNum == shp_val);
    
    if sum(idx(:)) > 0  % 如果存在该SHPNum值的点
        mean_TrueCoh(i) = mean(Stack.TrueCoh(idx), 'omitnan');
        mean_BoxCoh11(i) = mean(Stack.BoxCoh_11(idx), 'omitnan');
        mean_AdpCoh(i) = mean(Stack.AdpCoh(idx), 'omitnan');
        mean_ModifiedCoh(i) = mean(Stack.ModifiedCoh(idx), 'omitnan');
        mean_BestEst(i) = mean(Stack.Best_Est(idx), 'omitnan');
    else
        mean_TrueCoh(i) = NaN;
        mean_BoxCoh11(i) = NaN;
        mean_AdpCoh(i) = NaN;
        mean_ModifiedCoh(i) = NaN;
        mean_BestEst(i) = NaN;
    end
end

% 绘制对比图
plot(shp_values, mean_TrueCoh, 'k-', 'LineWidth', 2, 'DisplayName', 'TrueCoh'); hold on;
plot(shp_values, mean_BoxCoh11, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11');
plot(shp_values, mean_AdpCoh, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh');
plot(shp_values, mean_ModifiedCoh, 'g-', 'LineWidth', 1.5, 'DisplayName', 'ModifiedCoh');
plot(shp_values, mean_BestEst, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Best\_Est');
hold off;

xlabel('SHPNum (统计同质性像素数)');
ylabel('相干性均值');
xlim([0 PixelNum]);
title('不同SHPNum下各估计方法的相干性均值对比');
legend('Location', 'best');
grid on;

%% 图二：横轴为SHPNum，纵轴比较各种估计方法之比。1: Box/Adp, 2: Box/Mod, 3: Adp/Mod
% 初始化存储比值和差值的数组
ratio_BoxAdp = zeros(size(shp_values));      % BoxCoh_11 / AdpCoh
ratio_BoxMod = zeros(size(shp_values));      % BoxCoh_11 / ModifiedCoh
ratio_AdpMod = zeros(size(shp_values));      % AdpCoh / ModifiedCoh
diff_BoxAdp = zeros(size(shp_values));       % BoxCoh_11 - AdpCoh
diff_BoxMod = zeros(size(shp_values));       % BoxCoh_11 - ModifiedCoh
diff_AdpMod = zeros(size(shp_values));       % AdpCoh - ModifiedCoh

% 对每个SHPNum值计算比值和差值
for i = 1:length(shp_values)
    shp_val = shp_values(i);
    % 找到SHPNum等于当前值的所有位置
    idx = (Stack.SHPNum == shp_val);
    
    if sum(idx(:)) > 0  % 如果存在该SHPNum值的点
        % 计算均值
        mean_box = mean(Stack.BoxCoh_11(idx), 'omitnan');
        mean_adp = mean(Stack.AdpCoh(idx), 'omitnan');
        mean_mod = mean(Stack.ModifiedCoh(idx), 'omitnan');
        
        % 计算比值（避免除以0）
        if mean_adp > 0
            ratio_BoxAdp(i) = mean_box / mean_adp;
        else
            ratio_BoxAdp(i) = NaN;
        end
        
        if mean_mod > 0
            ratio_BoxMod(i) = mean_box / mean_mod;
        else
            ratio_BoxMod(i) = NaN;
        end
        
        if mean_mod > 0
            ratio_AdpMod(i) = mean_adp / mean_mod;
        else
            ratio_AdpMod(i) = NaN;
        end
        
        % 计算差值
        diff_BoxAdp(i) = mean_box - mean_adp;
        diff_BoxMod(i) = mean_box - mean_mod;
        diff_AdpMod(i) = mean_adp - mean_mod;
    else
        ratio_BoxAdp(i) = NaN;
        ratio_BoxMod(i) = NaN;
        ratio_AdpMod(i) = NaN;
        diff_BoxAdp(i) = NaN;
        diff_BoxMod(i) = NaN;
        diff_AdpMod(i) = NaN;
    end
end

% 绘制比值和差值对比图
figure;

% 子图1：比值对比图
subplot(1, 2, 1);
plot(shp_values, ones(size(shp_values)), 'k--', 'LineWidth', 1, 'DisplayName', '参考线 (y=1)'); hold on;
plot(shp_values, ratio_BoxAdp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 / AdpCoh');
plot(shp_values, ratio_BoxMod, 'r-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 / ModifiedCoh');
plot(shp_values, ratio_AdpMod, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh / ModifiedCoh');
hold off;

xlabel('SHPNum (统计同质性像素数)');
ylabel('相干性估计比值');
xlim([0 PixelNum]);
ylim([0 2]);
title('不同SHPNum下各估计方法之间的比值对比');
legend('Location', 'best');
grid on;

% 子图2：差值对比图
subplot(1, 2, 2);
plot(shp_values, zeros(size(shp_values)), 'k--', 'LineWidth', 1, 'DisplayName', '参考线 (y=0)'); hold on;
plot(shp_values, diff_BoxAdp, 'b-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 - AdpCoh');
plot(shp_values, diff_BoxMod, 'r-', 'LineWidth', 1.5, 'DisplayName', 'BoxCoh\_11 - ModifiedCoh');
plot(shp_values, diff_AdpMod, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh - ModifiedCoh');
hold off;

xlabel('SHPNum (统计同质性像素数)');
ylabel('相干性估计差值');
xlim([0 PixelNum]);
title('不同SHPNum下各估计方法之间的差值对比');
legend('Location', 'best');
grid on;

%% 图三：横轴为SHPNum，纵轴为相干性估计值，绘制箱线图，包括BoxCoh_11、AdpCoh、ModifiedCoh、TrueCoh、BestEst
% SHPNum从2开始，每隔15个设置为一个区间，2-16, 17-31, ..., 107-121

% 定义SHPNum区间
intervals = [2:15:107; 16:15:121]';  % 每行是一个区间 [起始, 结束]
num_intervals = size(intervals, 1);

% 为每个估计方法准备数据
data_TrueCoh = cell(num_intervals, 1);
data_BoxCoh11 = cell(num_intervals, 1);
data_AdpCoh = cell(num_intervals, 1);
data_ModifiedCoh = cell(num_intervals, 1);
data_BestEst = cell(num_intervals, 1);

% 对每个区间收集数据
for i = 1:num_intervals
    interval_start = intervals(i, 1);
    interval_end = intervals(i, 2);
    
    % 找到SHPNum在当前区间内的所有位置
    idx = (Stack.SHPNum >= interval_start) & (Stack.SHPNum <= interval_end);
    
    if sum(idx(:)) > 0  % 如果存在该区间内的点
        data_TrueCoh{i} = Stack.TrueCoh(idx);
        data_BoxCoh11{i} = Stack.BoxCoh_11(idx);
        data_AdpCoh{i} = Stack.AdpCoh(idx);
        data_ModifiedCoh{i} = Stack.ModifiedCoh(idx);
        data_BestEst{i} = Stack.Best_Est(idx);
    else
        data_TrueCoh{i} = [];
        data_BoxCoh11{i} = [];
        data_AdpCoh{i} = [];
        data_ModifiedCoh{i} = [];
        data_BestEst{i} = [];
    end
end

% 准备箱线图的标签
interval_labels = arrayfun(@(i) sprintf('%d-%d', intervals(i,1), intervals(i,2)), ...
                          1:num_intervals, 'UniformOutput', false);

% 将数据转换为boxplot所需的格式（数据向量+分组向量）
% 先计算总数据量以预分配数组
total_len_TrueCoh = sum(cellfun(@length, data_TrueCoh));
total_len_BoxCoh11 = sum(cellfun(@length, data_BoxCoh11));
total_len_AdpCoh = sum(cellfun(@length, data_AdpCoh));
total_len_ModifiedCoh = sum(cellfun(@length, data_ModifiedCoh));
total_len_BestEst = sum(cellfun(@length, data_BestEst));

% TrueCoh - 预分配
data_vec_TrueCoh = zeros(total_len_TrueCoh, 1);
group_vec_TrueCoh = zeros(total_len_TrueCoh, 1);
idx_start = 1;
for i = 1:num_intervals
    len = length(data_TrueCoh{i});
    if len > 0
        data_vec_TrueCoh(idx_start:idx_start+len-1) = data_TrueCoh{i}(:);
        group_vec_TrueCoh(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% BoxCoh_11 - 预分配
data_vec_BoxCoh11 = zeros(total_len_BoxCoh11, 1);
group_vec_BoxCoh11 = zeros(total_len_BoxCoh11, 1);
idx_start = 1;
for i = 1:num_intervals
    len = length(data_BoxCoh11{i});
    if len > 0
        data_vec_BoxCoh11(idx_start:idx_start+len-1) = data_BoxCoh11{i}(:);
        group_vec_BoxCoh11(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% AdpCoh - 预分配
data_vec_AdpCoh = zeros(total_len_AdpCoh, 1);
group_vec_AdpCoh = zeros(total_len_AdpCoh, 1);
idx_start = 1;
for i = 1:num_intervals
    len = length(data_AdpCoh{i});
    if len > 0
        data_vec_AdpCoh(idx_start:idx_start+len-1) = data_AdpCoh{i}(:);
        group_vec_AdpCoh(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% ModifiedCoh - 预分配
data_vec_ModifiedCoh = zeros(total_len_ModifiedCoh, 1);
group_vec_ModifiedCoh = zeros(total_len_ModifiedCoh, 1);
idx_start = 1;
for i = 1:num_intervals
    len = length(data_ModifiedCoh{i});
    if len > 0
        data_vec_ModifiedCoh(idx_start:idx_start+len-1) = data_ModifiedCoh{i}(:);
        group_vec_ModifiedCoh(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% BestEst - 预分配
data_vec_BestEst = zeros(total_len_BestEst, 1);
group_vec_BestEst = zeros(total_len_BestEst, 1);
idx_start = 1;
for i = 1:num_intervals
    len = length(data_BestEst{i});
    if len > 0
        data_vec_BestEst(idx_start:idx_start+len-1) = data_BestEst{i}(:);
        group_vec_BestEst(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% 绘制综合箱线图
figure;
hold on;

% 定义颜色
color_TrueCoh = [0 0 0];  % 黑色
color_BoxCoh11 = [0 0.447 0.741];  % 蓝色
color_AdpCoh = [0.850 0.325 0.098];  % 红色
color_ModifiedCoh = [0.466 0.674 0.188];  % 绿色
color_BestEst = [0.929 0.694 0.125];  % 黄色

% 定义每组箱线图的间距
group_width = 5;  % 每个SHPNum区间占据的x轴宽度
box_width = 1;    % 每个箱线图的宽度
box_spacing = 0.2;  % 箱线图之间的间距

% 为每个区间绘制五个箱线图
for i = 1:num_intervals
    % 计算当前区间五个箱线图的中心位置
    center_pos = i * group_width;
    positions = center_pos + [-2*box_width, -1*box_width, 0, 1*box_width, 2*box_width];
    
    % 绘制TrueCoh箱线图
    if ~isempty(data_TrueCoh{i})
        h1 = boxplot(data_TrueCoh{i}, 'Positions', positions(1), 'Widths', box_width*0.8, ...
                     'Colors', color_TrueCoh, 'Symbol', '');
        set(h1, 'LineWidth', 1.5);
    end
    
    % 绘制BoxCoh_11箱线图
    if ~isempty(data_BoxCoh11{i})
        h2 = boxplot(data_BoxCoh11{i}, 'Positions', positions(2), 'Widths', box_width*0.8, ...
                     'Colors', color_BoxCoh11, 'Symbol', '');
        set(h2, 'LineWidth', 1.5);
    end
    
    % 绘制AdpCoh箱线图
    if ~isempty(data_AdpCoh{i})
        h3 = boxplot(data_AdpCoh{i}, 'Positions', positions(3), 'Widths', box_width*0.8, ...
                     'Colors', color_AdpCoh, 'Symbol', '');
        set(h3, 'LineWidth', 1.5);
    end
    
    % 绘制ModifiedCoh箱线图
    if ~isempty(data_ModifiedCoh{i})
        h4 = boxplot(data_ModifiedCoh{i}, 'Positions', positions(4), 'Widths', box_width*0.8, ...
                     'Colors', color_ModifiedCoh, 'Symbol', '');
        set(h4, 'LineWidth', 1.5);
    end
    
    % 绘制BestEst箱线图
    if ~isempty(data_BestEst{i})
        h5 = boxplot(data_BestEst{i}, 'Positions', positions(5), 'Widths', box_width*0.8, ...
                     'Colors', color_BestEst, 'Symbol', '');
        set(h5, 'LineWidth', 1.5);
    end
end

% 设置x轴刻度和标签
x_tick_positions = (1:num_intervals) * group_width;
set(gca, 'XTick', x_tick_positions);
set(gca, 'XTickLabel', interval_labels);
xtickangle(45);

% 设置坐标轴标签和标题
xlabel('SHPNum区间');
ylabel('相干性估计值');
title('不同SHPNum区间下各估计方法的分布对比（综合箱线图）');
grid on;

% 添加图例
% 创建虚拟的线条用于图例
legend_handles = [
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_TrueCoh, 'MarkerEdgeColor', color_TrueCoh, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_BoxCoh11, 'MarkerEdgeColor', color_BoxCoh11, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_AdpCoh, 'MarkerEdgeColor', color_AdpCoh, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_ModifiedCoh, 'MarkerEdgeColor', color_ModifiedCoh, 'MarkerSize', 10);
    plot(NaN, NaN, 's', 'MarkerFaceColor', color_BestEst, 'MarkerEdgeColor', color_BestEst, 'MarkerSize', 10);
];
legend(legend_handles, {'TrueCoh', 'BoxCoh\_11', 'AdpCoh', 'ModifiedCoh', 'Best\_Est'}, 'Location', 'best');

hold off;

%% 图四：横轴为SHPNum，纵轴为CV的均值（改进后仅保留一个CV）
mean_CV = zeros(size(shp_values));

% 对每个SHPNum值计算均值
for i = 1:length(shp_values)
    shp_val = shp_values(i);
    % 找到SHPNum等于当前值的所有位置
    idx = (Stack.SHPNum == shp_val);
    
    if sum(idx(:)) > 0  % 如果存在该SHPNum值的点
        mean_CV(i) = mean(Stack.CV(idx), 'omitnan');
    else
        mean_CV(i) = NaN;
    end
end

% 绘制对比图
figure;
plot(shp_values, mean_CV, 'b-', 'LineWidth', 1.5, 'DisplayName', 'CV');
xlabel('SHPNum (统计同质性像素数)');
ylabel('变异系数均值');
xlim([0 PixelNum]);
title('不同SHPNum下CV的均值');
legend('Location', 'best');
grid on;

%% 图五：最佳窗口与SHPNum的关系图
% 子图1：不同最佳窗口大小对应的SHPNum分布（箱线图）
% 子图2：SHPNum区间内各窗口大小的占比（堆叠条形图）

% 获取所有唯一的窗口大小值并排序，过滤掉0和NaN
unique_wins = unique(Stack.Best_Win(:));
unique_wins = unique_wins(~isnan(unique_wins) & unique_wins > 0);
unique_wins = sort(unique_wins);
num_wins = length(unique_wins);

figure('Position', [100, 100, 1400, 900]);

% 子图1：对于每个窗口大小，显示对应的SHPNum分布
subplot(2, 2, 1);
hold on;

% 为每个窗口大小收集SHPNum数据
win_shp_data = cell(num_wins, 1);
total_len = 0;
for i = 1:num_wins
    win_val = unique_wins(i);
    idx = (Stack.Best_Win == win_val);
    win_shp_data{i} = Stack.SHPNum(idx);
    total_len = total_len + length(win_shp_data{i});
end

% 将数据转换为boxplot所需的格式（数据向量+分组向量）
data_vec = zeros(total_len, 1);
group_vec = zeros(total_len, 1);
idx_start = 1;
for i = 1:num_wins
    len = length(win_shp_data{i});
    if len > 0
        data_vec(idx_start:idx_start+len-1) = win_shp_data{i}(:);
        group_vec(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% 绘制箱线图
h_box = boxplot(data_vec, group_vec, 'Labels', arrayfun(@(x) sprintf('%d×%d', x, x), unique_wins, 'UniformOutput', false), ...
                 'Colors', 'b', 'Symbol', '');
set(h_box, 'LineWidth', 1.5);

% 叠加均值曲线
mean_shp_per_win = zeros(num_wins, 1);
for i = 1:num_wins
    mean_shp_per_win(i) = mean(win_shp_data{i}, 'omitnan');
end
plot(1:num_wins, mean_shp_per_win, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r', 'DisplayName', '均值');

xlabel('最佳窗口大小');
ylabel('SHPNum (统计同质性像素数)');
title('不同最佳窗口大小对应的SHPNum分布');
legend('Location', 'best');
grid on;
hold off;

% 子图2：在不同SHPNum区间内，各窗口大小的占比
subplot(2, 2, 2);

% 使用与图三相同的SHPNum区间
% intervals = [2:15:107; 16:15:121]';  % 已在前面定义

% 统计每个区间内各窗口大小的数量
win_count_matrix = zeros(num_intervals, num_wins);
for i = 1:num_intervals
    interval_start = intervals(i, 1);
    interval_end = intervals(i, 2);
    
    for j = 1:num_wins
        win_val = unique_wins(j);
        % 找到同时满足SHPNum在区间内且Best_Win为当前窗口大小的位置
        idx = (Stack.SHPNum >= interval_start) & (Stack.SHPNum <= interval_end) & ...
              (Stack.Best_Win == win_val);
        win_count_matrix(i, j) = sum(idx(:));
    end
end

% 计算百分比
win_percent_matrix = win_count_matrix ./ sum(win_count_matrix, 2) * 100;
% 处理可能的NaN（某些区间可能没有数据）
win_percent_matrix(isnan(win_percent_matrix)) = 0;

% 绘制堆叠条形图
h_bar = bar(1:num_intervals, win_percent_matrix, 'stacked');

% 设置每个条形的颜色
colors = lines(num_wins);  % 使用MATLAB默认的颜色方案
for i = 1:num_wins
    h_bar(i).FaceColor = colors(i, :);
end

% 设置x轴标签
set(gca, 'XTick', 1:num_intervals);
set(gca, 'XTickLabel', interval_labels);
xtickangle(45);

xlabel('SHPNum区间');
ylabel('占比 (%)');
title('不同SHPNum区间内各窗口大小的占比分布');
legend(arrayfun(@(x) sprintf('%d×%d', x, x), unique_wins, 'UniformOutput', false), ...
       'Location', 'bestoutside');
grid on;

% 子图3：横轴为SHPNum，纵轴为对应SHPNum时最佳窗口的众数
subplot(2, 2, 3);

% 对每个SHPNum值计算最佳窗口的众数
mode_win = zeros(size(shp_values));
for i = 1:length(shp_values)
    shp_val = shp_values(i);
    % 找到SHPNum等于当前值的所有位置
    idx = (Stack.SHPNum == shp_val);
    
    if sum(idx(:)) > 0  % 如果存在该SHPNum值的点
        % 获取对应的Best_Win值（排除0和NaN）
        best_wins = Stack.Best_Win(idx);
        best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
        
        if ~isempty(best_wins)
            mode_win(i) = mode(best_wins);  % 计算众数
        else
            mode_win(i) = NaN;
        end
    else
        mode_win(i) = NaN;
    end
end

% 绘制折线图
plot(shp_values, mode_win, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, ...
     'MarkerFaceColor', 'b');
xlabel('SHPNum (统计同质性像素数)');
ylabel('最佳窗口大小（众数）');
xlim([0 PixelNum]);
title('不同SHPNum下最佳窗口大小的众数');
grid on;

% 设置y轴刻度为整数
yticks(unique_wins);
yticklabels(arrayfun(@(x) sprintf('%d', x), unique_wins, 'UniformOutput', false));

% 子图4：横轴为SHPNum，纵轴为对应SHPNum时最佳窗口的平均值
subplot(2, 2, 4);

% 对每个SHPNum值计算最佳窗口的平均值
mean_win = zeros(size(shp_values));
for i = 1:length(shp_values)
    shp_val = shp_values(i);
    % 找到SHPNum等于当前值的所有位置
    idx = (Stack.SHPNum == shp_val);
    
    if sum(idx(:)) > 0  % 如果存在该SHPNum值的点
        % 获取对应的Best_Win值（排除0和NaN）
        best_wins = Stack.Best_Win(idx);
        best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
        
        if ~isempty(best_wins)
            mean_win(i) = mean(best_wins);  % 计算平均值
        else
            mean_win(i) = NaN;
        end
    else
        mean_win(i) = NaN;
    end
end

% 绘制折线图
plot(shp_values, mean_win, 'r-', 'LineWidth', 2);
xlabel('SHPNum (统计同质性像素数)');
ylabel('最佳窗口大小（平均值）');
xlim([0 PixelNum]);
title('不同SHPNum下最佳窗口大小的平均值');
grid on;

% 整体标题
sgtitle('最佳窗口大小与SHPNum的关系分析', 'FontSize', 14, 'FontWeight', 'bold');


%% 图六：最佳窗口与CV的关系图
% 子图1：不同最佳窗口大小对应的CV分布（箱线图）
% 子图2：CV区间内各窗口大小的占比（堆叠条形图）
% 子图3：横轴为CV，纵轴为最佳窗口的众数
% 子图4：横轴为CV，纵轴为最佳窗口的平均值

figure('Position', [100, 100, 1400, 900]);

% 子图1：对于每个窗口大小，显示对应的CV分布
subplot(2, 2, 1);
hold on;

% 为每个窗口大小收集CV数据
win_cv_data = cell(num_wins, 1);
total_len_cv = 0;
for i = 1:num_wins
    win_val = unique_wins(i);
    idx = (Stack.Best_Win == win_val);
    win_cv_data{i} = Stack.CV(idx);
    total_len_cv = total_len_cv + length(win_cv_data{i});
end

% 将数据转换为boxplot所需的格式（数据向量+分组向量）
data_vec_cv = zeros(total_len_cv, 1);
group_vec_cv = zeros(total_len_cv, 1);
idx_start = 1;
for i = 1:num_wins
    len = length(win_cv_data{i});
    if len > 0
        data_vec_cv(idx_start:idx_start+len-1) = win_cv_data{i}(:);
        group_vec_cv(idx_start:idx_start+len-1) = i;
        idx_start = idx_start + len;
    end
end

% 绘制箱线图
h_box_cv = boxplot(data_vec_cv, group_vec_cv, 'Labels', arrayfun(@(x) sprintf('%d×%d', x, x), unique_wins, 'UniformOutput', false), ...
                 'Colors', 'b', 'Symbol', '');
set(h_box_cv, 'LineWidth', 1.5);

% 叠加均值曲线
mean_cv_per_win = zeros(num_wins, 1);
for i = 1:num_wins
    mean_cv_per_win(i) = mean(win_cv_data{i}, 'omitnan');
end
plot(1:num_wins, mean_cv_per_win, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r', 'DisplayName', '均值');

xlabel('最佳窗口大小');
ylabel('CV');
title('不同最佳窗口大小对应的CV分布');
legend('Location', 'best');
grid on;
hold off;

% 子图2：在不同CV区间内，各窗口大小的占比
subplot(2, 2, 2);

% 定义CV区间（按分位数划分，每个区间包含12.5%的像元）
cv_all = Stack.CV(:);
cv_all = cv_all(~isnan(cv_all));
cv_min = min(cv_all);
cv_max = max(cv_all);
% 将CV按分位数划分为8个区间，每个区间包含12.5%的像元
num_cv_intervals = 8;
percentiles = linspace(0, 100, num_cv_intervals + 1);
cv_edges = prctile(cv_all, percentiles);
cv_intervals = [cv_edges(1:end-1)', cv_edges(2:end)'];

% 统计每个区间内各窗口大小的数量
win_count_matrix_cv = zeros(num_cv_intervals, num_wins);
for i = 1:num_cv_intervals
    interval_start = cv_intervals(i, 1);
    interval_end = cv_intervals(i, 2);
    
    for j = 1:num_wins
        win_val = unique_wins(j);
        % 找到同时满足CV在区间内且Best_Win为当前窗口大小的位置
        idx = (Stack.CV >= interval_start) & (Stack.CV <= interval_end) & ...
              (Stack.Best_Win == win_val);
        win_count_matrix_cv(i, j) = sum(idx(:));
    end
end

% 计算百分比
win_percent_matrix_cv = win_count_matrix_cv ./ sum(win_count_matrix_cv, 2) * 100;
% 处理可能的NaN（某些区间可能没有数据）
win_percent_matrix_cv(isnan(win_percent_matrix_cv)) = 0;

% 绘制堆叠条形图
h_bar_cv = bar(1:num_cv_intervals, win_percent_matrix_cv, 'stacked');

% 设置每个条形的颜色
for i = 1:num_wins
    h_bar_cv(i).FaceColor = colors(i, :);
end

% 设置x轴标签
cv_interval_labels = arrayfun(@(i) sprintf('%.2f-%.2f', cv_intervals(i,1), cv_intervals(i,2)), ...
                          1:num_cv_intervals, 'UniformOutput', false);
set(gca, 'XTick', 1:num_cv_intervals);
set(gca, 'XTickLabel', cv_interval_labels);
xtickangle(45);

xlabel('CV区间');
ylabel('占比 (%)');
title('不同CV区间内各窗口大小的占比分布');
legend(arrayfun(@(x) sprintf('%d×%d', x, x), unique_wins, 'UniformOutput', false), ...
       'Location', 'bestoutside');
grid on;

% 子图3：横轴为CV，纵轴为对应CV时最佳窗口的众数
subplot(2, 2, 3);

% 创建CV的采样点（使用更细的分辨率）
cv_sample_points = linspace(cv_min, cv_max, 100);
mode_win_cv = zeros(size(cv_sample_points));

% 对每个CV采样点，使用一个小窗口内的数据计算众数
cv_window_size = (cv_max - cv_min) / 50;  % 窗口宽度
for i = 1:length(cv_sample_points)
    cv_val = cv_sample_points(i);
    % 找到CV在当前窗口内的所有位置
    idx = (Stack.CV >= cv_val - cv_window_size/2) & ...
          (Stack.CV <= cv_val + cv_window_size/2);
    
    if sum(idx(:)) > 0  % 如果存在该CV窗口内的点
        % 获取对应的Best_Win值（排除0和NaN）
        best_wins = Stack.Best_Win(idx);
        best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
        
        if ~isempty(best_wins)
            mode_win_cv(i) = mode(best_wins);  % 计算众数
        else
            mode_win_cv(i) = NaN;
        end
    else
        mode_win_cv(i) = NaN;
    end
end

% 绘制折线图
plot(cv_sample_points, mode_win_cv, 'b-', 'LineWidth', 1.5);
xlabel('CV');
ylabel('最佳窗口大小（众数）');
title('不同CV下最佳窗口大小的众数');
grid on;

% 设置y轴刻度为整数
yticks(unique_wins);
yticklabels(arrayfun(@(x) sprintf('%d', x), unique_wins, 'UniformOutput', false));

% 子图4：横轴为CV，纵轴为对应CV时最佳窗口的平均值
subplot(2, 2, 4);

% 对每个CV采样点计算最佳窗口的平均值
mean_win_cv = zeros(size(cv_sample_points));
for i = 1:length(cv_sample_points)
    cv_val = cv_sample_points(i);
    % 找到CV在当前窗口内的所有位置
    idx = (Stack.CV >= cv_val - cv_window_size/2) & ...
          (Stack.CV <= cv_val + cv_window_size/2);
    
    if sum(idx(:)) > 0  % 如果存在该CV窗口内的点
        % 获取对应的Best_Win值（排除0和NaN）
        best_wins = Stack.Best_Win(idx);
        best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
        
        if ~isempty(best_wins)
            mean_win_cv(i) = mean(best_wins);  % 计算平均值
        else
            mean_win_cv(i) = NaN;
        end
    else
        mean_win_cv(i) = NaN;
    end
end

% 绘制折线图
plot(cv_sample_points, mean_win_cv, 'r-', 'LineWidth', 2);
xlabel('CV');
ylabel('最佳窗口大小（平均值）');
title('不同CV下最佳窗口大小的平均值');
grid on;

% 整体标题
sgtitle('最佳窗口大小与CV的关系分析', 'FontSize', 14, 'FontWeight', 'bold');