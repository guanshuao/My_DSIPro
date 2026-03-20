%% 将Stack.mat文件中的数据可视化，主要是与CV相关的可视化
close all;clear;clc;

data_dir = "D:\Research_Material\TGRS2026\Beijing_ML\Stack.mat";
% data_dir = "D:\Research_Material\TGRS2026\Beijing_SL\Stack.mat";
% data_dir = "D:\Research_Material\TGRS2026\Kumamoto\Stack.mat";

disp('Loading Stack.mat...');
load(data_dir);
disp('Loaded Stack.mat...');

%% 设置想要展示的图的ID
figure_id = [1]; % 例如 [1,3,4] 表示只展示图1、3、4

% - Stack结构体包含多个字段，具体如下：
% - Calwin: 用于选取SHP、计算CV、计算RCS的窗口的大小矩阵，尺寸都一样，如 [11 11]
% - RCS：雷达散射截面（Radar Cross Section, RCS）矩阵，由强度图像经过滤波得到，滤波窗口为Calwin中指定的窗口大小
% - CV：由RCS矩阵计算得到的变异系数矩阵，衡量地表的均匀程度，CV值越大表示地表越不均匀，使用的窗口为Calwin中指定的窗口大小
% - SHPNum: 每个位置的周围的统计同质性像素（Statistically Homogeneous Pixels）数矩阵，例如某元素值为50表示该位置周围有50个同质性像素。选取SHP时使用的窗口大小为Calwin中指定的窗口大小
% - AdpCph: 自适应相干系数估计矩阵。即估计相干性时只使用统计同质性像素的数据进行估计，窗口大小为Calwin中指定的窗口大小
% - TrueCoh: 真实的相干矩阵
% - BoxCoh_XX: 使用窗口大小为XX的盒式相干估计，例如BoxCoh_3表示窗口大小为3x3估计出的相干矩阵，XX取值为3, 5……。如果Calwin中指定的窗口大小为[11 11]，则会有BoxCoh_3、BoxCoh_5、BoxCoh_7、BoxCoh_9、BoxCoh_11五个字段，分别表示使用不同窗口大小估计出的相干矩阵
% - Best_Est: 同一个位置中，所有BoxCoh_XX估计中与TrueCoh最接近的估计矩阵
% - Best_Win: 对应Best_Est的窗口大小矩阵，表示每个位置使用的最佳窗口大小，如某一位置Best_Win值为5表示该位置使用5x5窗口估计出的相干矩阵最接近真实值TrueCoh
% - ModifiedCoh: 修正后的相干矩阵
% - level_num: 不同CV和SHPNum值分布数量不均一，因此需要将CV和SHPNum的取值映射到1~level_num的整数范围内，level_num由用户指定，例如level_num=50表示将CV和SHPNum的取值映射到1~50的整数范围内
% - CV_level：将CV的取值映射到1~level_num的整数范围内的矩阵，数值越大表示CV值越大，即地表越不均匀
% - SHPNum_level：将SHPNum的取值映射到1~level_num的整数范围内的矩阵，数值越大表示SHPNum值越大，即周围同质性像素数量越多
% - 其余字段请忽略


%% 获取基础信息
level_num = Stack.level_num;
levels = 1:level_num;
Calwin = Stack.Calwin;
max_win = Calwin(1);
box_coh_field = sprintf('BoxCoh_%d', max_win);
if isfield(Stack, box_coh_field)
    BoxCoh_Max = Stack.(box_coh_field);
else
    BoxCoh_Max = NaN(size(Stack.TrueCoh));
end

% 动态获取所有的 BoxCoh_XX 字段
fields = fieldnames(Stack);
box_coh_fields = fields(startsWith(fields, 'BoxCoh_'));
% 提取窗口大小并排序
box_wins = zeros(length(box_coh_fields), 1);
for i = 1:length(box_coh_fields)
    tokens = regexp(box_coh_fields{i}, 'BoxCoh_(\d+)', 'tokens');
    box_wins(i) = str2double(tokens{1}{1});
end
[box_wins, sort_idx] = sort(box_wins);
box_coh_fields = box_coh_fields(sort_idx);

%% 图一：不同CV Level下各估计方法的相干性均值对比
if ismember(1, figure_id)
    mean_TrueCoh = zeros(1, level_num);
    mean_AdpCoh = zeros(1, level_num);
    mean_BestEst = zeros(1, level_num);
    mean_ModifiedCoh = zeros(1, level_num);
    mean_BoxCohMax = zeros(1, level_num);
    
    for i = 1:level_num
        idx = (Stack.CV_level == i);
        if any(idx(:))
            mean_TrueCoh(i) = mean(Stack.TrueCoh(idx), 'omitnan');
            if isfield(Stack, 'AdpCph')
                mean_AdpCoh(i) = mean(Stack.AdpCph(idx), 'omitnan');
            elseif isfield(Stack, 'AdpCoh')
                mean_AdpCoh(i) = mean(Stack.AdpCoh(idx), 'omitnan');
            end
            mean_BestEst(i) = mean(Stack.Best_Est(idx), 'omitnan');
            if isfield(Stack, 'ModifiedCoh')
                mean_ModifiedCoh(i) = mean(Stack.ModifiedCoh(idx), 'omitnan');
            end
            mean_BoxCohMax(i) = mean(BoxCoh_Max(idx), 'omitnan');
        else
            mean_TrueCoh(i) = NaN;
            mean_AdpCoh(i) = NaN;
            mean_BestEst(i) = NaN;
            mean_ModifiedCoh(i) = NaN;
            mean_BoxCohMax(i) = NaN;
        end
    end
    
    figure('Name', '不同CV Level下的各估计方法均值对比', 'Position', [100, 100, 800, 600]);
    plot(levels, mean_TrueCoh, 'k-', 'LineWidth', 2, 'DisplayName', 'TrueCoh'); hold on;
    plot(levels, mean_BoxCohMax, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('BoxCoh\\_%d', max_win));
    plot(levels, mean_AdpCoh, 'r-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh');
    if isfield(Stack, 'ModifiedCoh')
        plot(levels, mean_ModifiedCoh, 'g-', 'LineWidth', 1.5, 'DisplayName', 'ModifiedCoh');
    end
    plot(levels, mean_BestEst, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Best\_Est');
    hold off;
    
    xlabel('CV Level');
    ylabel('相干性均值');
    xlim([1 level_num]);
    title('不同 CV Level 下各估计方法的相干性均值对比');
    legend('Location', 'best');
    grid on;
    
    clear idx mean_TrueCoh mean_AdpCoh mean_BestEst mean_ModifiedCoh mean_BoxCohMax;
end

%% 图二：不同CV Level下的各尺寸窗口估计均值对比
if ismember(2, figure_id)
    mean_BoxCoh_all = zeros(length(box_wins), level_num);
    mean_TrueCoh = zeros(1, level_num);
    
    for i = 1:level_num
        idx = (Stack.CV_level == i);
        if any(idx(:))
            mean_TrueCoh(i) = mean(Stack.TrueCoh(idx), 'omitnan');
            for j = 1:length(box_wins)
                mean_BoxCoh_all(j, i) = mean(Stack.(box_coh_fields{j})(idx), 'omitnan');
            end
        else
            mean_TrueCoh(i) = NaN;
            mean_BoxCoh_all(:, i) = NaN;
        end
    end
    
    figure('Name', '不同CV Level下的各尺寸窗口估计均值对比', 'Position', [150, 150, 800, 600]);
    plot(levels, mean_TrueCoh, 'k-', 'LineWidth', 2, 'DisplayName', 'TrueCoh (参考线)'); hold on;
    
    colors = lines(length(box_wins));
    for j = 1:length(box_wins)
        plot(levels, mean_BoxCoh_all(j, :), '-', 'Color', colors(j, :), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('BoxCoh\\_%d', box_wins(j)));
    end
    hold off;
    
    xlabel('CV Level');
    ylabel('相干性均值');
    xlim([1 level_num]);
    title('不同尺寸窗口估计方法在各 CV Level 下的相干性均值对比');
    legend('Location', 'best');
    grid on;

    clear idx mean_BoxCoh_all mean_TrueCoh colors;
end


%% 图三：不同CV Level下各估计方法之间的比值和差值比较
if ismember(3, figure_id)
    ratio_BoxAdp = zeros(1, level_num);
    ratio_BoxMod = zeros(1, level_num);
    ratio_AdpMod = zeros(1, level_num);
    diff_BoxAdp = zeros(1, level_num);
    diff_BoxMod = zeros(1, level_num);
    diff_AdpMod = zeros(1, level_num);
    
    for i = 1:level_num
        idx = (Stack.CV_level == i);
        if any(idx(:))
            box_mean = mean(BoxCoh_Max(idx), 'omitnan');
            if isfield(Stack, 'AdpCph')
                adp_mean = mean(Stack.AdpCph(idx), 'omitnan');
            else
                adp_mean = mean(Stack.AdpCoh(idx), 'omitnan');
            end
            
            ratio_BoxAdp(i) = box_mean / adp_mean;
            diff_BoxAdp(i) = box_mean - adp_mean;
            
            if isfield(Stack, 'ModifiedCoh')
                mod_mean = mean(Stack.ModifiedCoh(idx), 'omitnan');
                ratio_BoxMod(i) = box_mean / mod_mean;
                ratio_AdpMod(i) = adp_mean / mod_mean;
                diff_BoxMod(i) = box_mean - mod_mean;
                diff_AdpMod(i) = adp_mean - mod_mean;
            end
        else
            ratio_BoxAdp(i) = NaN; ratio_BoxMod(i) = NaN; ratio_AdpMod(i) = NaN;
            diff_BoxAdp(i) = NaN; diff_BoxMod(i) = NaN; diff_AdpMod(i) = NaN;
        end
    end
    
    figure('Name', '不同CV等级下各估计方法之间的比值和差值比较', 'Position', [200, 200, 1000, 500]);
    subplot(1, 2, 1);
    plot(levels, ones(1, level_num), 'k--', 'LineWidth', 1, 'DisplayName', '参考线 (y=1)'); hold on;
    plot(levels, ratio_BoxAdp, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('BoxCoh\\_%d / AdpCoh', max_win));
    if isfield(Stack, 'ModifiedCoh')
        plot(levels, ratio_BoxMod, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('BoxCoh\\_%d / ModifiedCoh', max_win));
        plot(levels, ratio_AdpMod, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh / ModifiedCoh');
    end
    hold off;
    xlabel('CV Level'); ylabel('相干性估计比值');
    xlim([1 level_num]); title('不同CV等级下各估计方法之间的比值对比'); legend('Location', 'best'); grid on;
    
    subplot(1, 2, 2);
    plot(levels, zeros(1, level_num), 'k--', 'LineWidth', 1, 'DisplayName', '参考线 (y=0)'); hold on;
    plot(levels, diff_BoxAdp, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('BoxCoh\\_%d - AdpCoh', max_win));
    if isfield(Stack, 'ModifiedCoh')
        plot(levels, diff_BoxMod, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('BoxCoh\\_%d - ModifiedCoh', max_win));
        plot(levels, diff_AdpMod, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AdpCoh - ModifiedCoh');
    end
    hold off;
    xlabel('CV Level'); ylabel('相干性估计差值');
    xlim([1 level_num]); title('不同CV等级下各估计方法之间的差值对比'); legend('Location', 'best'); grid on;

    clear idx box_mean adp_mean mod_mean ratio_BoxAdp ratio_BoxMod ratio_AdpMod diff_BoxAdp diff_BoxMod diff_AdpMod;
end

%% 图四：最佳窗口大小与CV Level的关系图
if ismember(4, figure_id)
% 子图1：不同最佳窗口大小对应的CV分布（箱线图）
% 子图2：CV等级内各窗口大小的占比（堆叠条形图）
% 子图3：横轴为CV等级，纵轴为最佳窗口的分布（箱线图）
% 子图4：横轴为CV等级，纵轴为对应CV等级时最佳窗口的平均值

% 获取所有唯一的窗口大小值并排序，过滤掉0和NaN
unique_wins = unique(Stack.Best_Win(:));
unique_wins = unique_wins(~isnan(unique_wins) & unique_wins > 0);
unique_wins = sort(unique_wins);
num_wins = length(unique_wins);

figure('Name', '最佳窗口大小与CV Level的关系', 'Position', [100, 100, 1400, 900]);

% 子图1：对于每个窗口大小，显示对应的CV Level分布
subplot(2, 2, 1);
hold on;

% 为每个窗口大小收集CV Level数据
win_cv_data = cell(num_wins, 1);
total_len_cv = 0;
for i = 1:num_wins
    win_val = unique_wins(i);
    idx = (Stack.Best_Win == win_val);
    win_cv_data{i} = Stack.CV_level(idx);
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
ylabel('CV Level');
title('不同最佳窗口大小对应的CV Level分布');
legend('Location', 'best');
grid on;
hold off;

% 子图2：在不同CV等级内，各窗口大小的占比
subplot(2, 2, 2);

% 统计每个level内各窗口大小的数量
win_count_matrix_cv = zeros(level_num, num_wins);
for i = 1:level_num
    level_idx = (Stack.CV_level == i);
    for j = 1:num_wins
        win_val = unique_wins(j);
        win_count_matrix_cv(i, j) = sum(Stack.Best_Win(level_idx) == win_val);
    end
end

% 计算百分比
row_sums = sum(win_count_matrix_cv, 2);
win_percent_matrix_cv = win_count_matrix_cv ./ row_sums * 100;
win_percent_matrix_cv(isnan(win_percent_matrix_cv)) = 0;

% 绘制堆叠条形图
colors_bar = lines(num_wins);
h_bar_cv = bar(1:level_num, win_percent_matrix_cv, 'stacked', 'EdgeColor', 'none');
for i = 1:num_wins
    h_bar_cv(i).FaceColor = colors_bar(i, :);
end

xlabel('CV Level');
ylabel('占比 (%)');
xlim([0 level_num+1]);
title('不同CV等级内各窗口大小的占比分布');
legend(arrayfun(@(x) sprintf('%d×%d', x, x), unique_wins, 'UniformOutput', false), ...
       'Location', 'bestoutside');
grid on;

% 计算每个CV Level下的均值和收集窗口数据
mean_win_cv = zeros(1, level_num);
win_data_cv = cell(level_num, 1);
total_win_cv = 0;

for i = 1:level_num
    level_idx = (Stack.CV_level == i);
    best_wins = Stack.Best_Win(level_idx);
    best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
    
    win_data_cv{i} = best_wins;
    total_win_cv = total_win_cv + length(best_wins);
    
    if ~isempty(best_wins)
        mean_win_cv(i) = mean(best_wins);
    else
        mean_win_cv(i) = NaN;
    end
end

% 子图3：横轴为CV Level，纵轴为对应的最佳窗口分布（箱线图）
subplot(2, 2, 3);
hold on;

data_vec_cv_win = zeros(total_win_cv, 1);
group_vec_cv_win = zeros(total_win_cv, 1);
idx_s = 1;
for i = 1:level_num
    len = length(win_data_cv{i});
    if len > 0
        data_vec_cv_win(idx_s:idx_s+len-1) = win_data_cv{i}(:);
        group_vec_cv_win(idx_s:idx_s+len-1) = i;
        idx_s = idx_s + len;
    end
end

h_box_cv_win = boxplot(data_vec_cv_win, group_vec_cv_win, 'Colors', 'b', 'Symbol', '');
set(h_box_cv_win, 'LineWidth', 1.0);

xtick_idx = 1:5:level_num;
set(gca, 'XTick', xtick_idx);
set(gca, 'XTickLabel', xtick_idx);

plot(1:level_num, mean_win_cv, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', '均值');
xlabel('CV Level');
ylabel('最佳窗口大小');
title('不同CV Level下最佳窗口大小的分布');
grid on;
hold off;

% 子图4：横轴为CV等级，纵轴为对应CV等级时最佳窗口的平均值
subplot(2, 2, 4);

% 绘制折线图
plot(1:level_num, mean_win_cv, 'r-', 'LineWidth', 2);
xlabel('CV Level');
ylabel('最佳窗口大小的平均值');
xlim([1 level_num]); ylim([0 max_win]);
title('不同CV Level下最佳窗口大小的平均值');
grid on;

% 整体标题
sgtitle('最佳窗口大小与CV Level的关系分析');
    
    clear win_cv_data data_vec_cv group_vec_cv mean_cv_per_win win_count_matrix_cv row_sums win_percent_matrix_cv h_bar_cv mean_win_cv win_data_cv data_vec_cv_win group_vec_cv_win h_box_cv_win;
end