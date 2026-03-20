%% 将Stack.mat文件中的数据可视化，主要是与二者同时相关的可视化
close all;clear;clc;

data_dir = "D:\Research_Material\TGRS2026\Beijing_ML\Stack.mat";
% data_dir = "D:\Research_Material\TGRS2026\Beijing_SL\Stack.mat";
% data_dir = "D:\Research_Material\TGRS2026\Kumamoto\Stack.mat";

disp('Loading Stack.mat...');
load(data_dir);
disp('Loaded Stack.mat...');

%% 设置想要展示的图的ID
figure_id = [1, 2, 3, 4]; % 例如 [1,3,4] 表示只展示图1、3、4

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

%% 图一：左右两个子图，左边子图为不同SHP Level下CV Level的均值，右图为不同CV等级下SHP Level的均值（参考visualization.m中的图四代码）
if ismember(1, figure_id)
    figure('Name', 'SHP Level与CV Level的互相关系', 'Position', [100, 100, 1000, 500]);
    
    % 左图：横轴SHP Level，纵轴CV Level的均值
    subplot(1, 2, 1);
    mean_CV_Level = zeros(1, level_num);
    for i = 1:level_num
        idx = (Stack.SHPNum_level == i);
        if any(idx(:))
            mean_CV_Level(i) = mean(Stack.CV_level(idx), 'omitnan');
        else
            mean_CV_Level(i) = NaN;
        end
    end
    plot(levels, mean_CV_Level, 'b-', 'LineWidth', 1.5, 'DisplayName', 'CV Level 均值');
    xlabel('SHP Level');
    ylabel('CV Level (均值)');
    xlim([1 level_num]); ylim([1 level_num]);
    title('不同 SHP Level 下的 CV Level 均值');
    grid on; legend('Location', 'best');
    
    % 右图：横轴CV Level，纵轴SHP Level的均值
    subplot(1, 2, 2);
    mean_SHP_Level = zeros(1, level_num);
    for i = 1:level_num
        idx = (Stack.CV_level == i);
        if any(idx(:))
            mean_SHP_Level(i) = mean(Stack.SHPNum_level(idx), 'omitnan');
        else
            mean_SHP_Level(i) = NaN;
        end
    end
    plot(levels, mean_SHP_Level, 'r-', 'LineWidth', 1.5, 'DisplayName', 'SHP Level 均值');
    xlabel('CV Level');
    ylabel('SHP Level (均值)');
    xlim([1 level_num]); ylim([1 level_num]);
    title('不同 CV Level 下的 SHP Level 均值');
    grid on; legend('Location', 'best');
    
    clear idx mean_CV_Level mean_SHP_Level;
end

%% 图二：SHP Level与CV Level分级下的BoxCoh估计误差热力图（参考visualization2.m中的图四代码）
if ismember(2, figure_id)
    error_matrix = NaN(level_num, level_num);
    count_matrix = zeros(level_num, level_num);
    
    for i = 1:level_num % SHP Level (横轴)
        for j = 1:level_num % CV Level (纵轴)
            idx = (Stack.SHPNum_level == i) & (Stack.CV_level == j);
            count_matrix(j, i) = sum(idx(:));
            if count_matrix(j, i) > 0
                % 计算绝对误差平均值: |BoxCoh_Max - TrueCoh|
                errors = BoxCoh_Max(idx) - Stack.TrueCoh(idx);
                error_matrix(j, i) = mean(errors, 'omitnan');
            end
        end
    end
    
    figure('Name', 'BoxCoh估计误差热力图', 'Position', [150, 150, 600, 500]);
    imagesc(levels, levels, error_matrix);
    set(gca, 'YDir', 'normal'); % 使Y轴递增方向从下到上
    colormap(jet);
    c = colorbar;
    c.Label.String = sprintf('BoxCoh_{%d} - TrueCoh', max_win);
    
    
    xlabel('SHP Level');
    ylabel('CV Level');
    title(sprintf('SHP Level与CV Level下的BoxCoh\\_%d估计误差热力图', max_win));
    grid off;
    
    % 添加等值线增强可读性
    hold on;
    [X, Y] = meshgrid(levels, levels);
    valid_mask = count_matrix >= 5; % 至少几个样本点
    error_matrix_contour = error_matrix;
    error_matrix_contour(~valid_mask) = NaN;
    contour(X, Y, error_matrix_contour, 8, 'k', 'LineWidth', 0.5);
    hold off;
    
    clear error_matrix count_matrix idx errors c X Y valid_mask error_matrix_contour;
end

%% 图三：SHP Level与CV Level分级下的最佳窗口大小平均值热力图（参考visualization2.m中的图五代码）
if ismember(3, figure_id)
    best_win_matrix = NaN(level_num, level_num);
    count_matrix_win = zeros(level_num, level_num);
    
    for i = 1:level_num % SHP Level (横轴)
        for j = 1:level_num % CV Level (纵轴)
            idx = (Stack.SHPNum_level == i) & (Stack.CV_level == j);
            count_matrix_win(j, i) = sum(idx(:));
            if count_matrix_win(j, i) > 0
                best_win_matrix(j, i) = mean(Stack.Best_Win(idx), 'omitnan');
            end
        end
    end
    
    figure('Name', '最佳窗口大小平均值热力图', 'Position', [200, 200, 600, 500]);
    imagesc(levels, levels, best_win_matrix);
    set(gca, 'YDir', 'normal'); 
    colormap(jet);
    c = colorbar;
    c.Label.String = '最佳窗口大小平均值 (Best\_Win)';
    
    xlabel('SHP Level');
    ylabel('CV Level');
    title('SHP Level与CV Level分级下的最佳窗口大小平均值热力图');
    grid off;
    
    hold on;
    [X5, Y5] = meshgrid(levels, levels);
    valid_mask_win = count_matrix_win >= 5; 
    best_win_contour = best_win_matrix;
    best_win_contour(~valid_mask_win) = NaN;
    contour(X5, Y5, best_win_contour, 8, 'k', 'LineWidth', 0.5);
    hold off;
    
    clear best_win_matrix count_matrix_win idx c X5 Y5 valid_mask_win best_win_contour;
end

%% 图四：SHP Level与CV Level分级下的像元数量热力图
if ismember(4, figure_id)
    count_matrix = zeros(level_num, level_num);
    for i = 1:level_num % SHP Level (横轴)
        for j = 1:level_num % CV Level (纵轴)
            idx = (Stack.SHPNum_level == i) & (Stack.CV_level == j);
            count_matrix(j, i) = sum(idx(:));
        end
    end

    figure('Name', '像元数量热力图', 'Position', [250, 250, 600, 500]);
    imagesc(levels, levels, count_matrix);
    set(gca, 'YDir', 'normal'); 
    colormap(jet);
    c = colorbar;
    c.Label.String = '像元数量';
    
    xlabel('SHP Level');
    ylabel('CV Level');
    title('SHP Level与CV Level分级下的像元分布数量热力图');
    grid off;
    
    clear count_matrix idx c;
end