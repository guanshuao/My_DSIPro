% 将Stack.mat文件中的数据可视化，主要是与二者同时相关的可视化

close all; clear; clc;

data_dir = "Stack.mat";
% data_dir = "D:\Research_Material\TGRS2026\Beijing_SL\Stack.mat";
% data_dir = "D:\Research_Material\TGRS2026\Kumamoto\Stack.mat";

disp('Loading Stack.mat...');
load(data_dir);
disp('Loaded Stack.mat...');

%% 设置想要展示的图的ID
figure_id = [1, 2, 3, 4];

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

% 提取有效像元
isValid = ~isnan(Stack.CV_level) & ~isnan(Stack.SHPNum_level) & ...
          (Stack.CV_level >= 1) & (Stack.CV_level <= level_num) & ...
          (Stack.SHPNum_level >= 1) & (Stack.SHPNum_level <= level_num);
      
CV_valid = Stack.CV_level(isValid);
SHP_valid = Stack.SHPNum_level(isValid);

% [行, 列] = [CV Level, SHP Level]
subs2D = [CV_valid(:), SHP_valid(:)];

% 1. 计算像元数量矩阵 (对应图四)
count_matrix = accumarray(subs2D, 1, [level_num, level_num], @sum, 0);

% 2. 计算误差均值矩阵 (对应图二)
errors_valid = double(BoxCoh_Max(isValid) - Stack.TrueCoh(isValid));
error_matrix = accumarray(subs2D, errors_valid(:), [level_num, level_num], @(x) mean(x, 'omitnan'), NaN);

% 3. 计算最佳窗口尺寸均值矩阵 (对应图三)
best_win_valid = double(Stack.Best_Win(isValid));
best_win_matrix = accumarray(subs2D, best_win_valid(:), [level_num, level_num], @(x) mean(x, 'omitnan'), NaN);

% 4. 计算一维序列 (对应图一)
mean_CV_Level = accumarray(double(SHP_valid(:)), double(CV_valid(:)), [level_num, 1], @(x) mean(x, 'omitnan'), NaN)';
mean_SHP_Level = accumarray(double(CV_valid(:)), double(SHP_valid(:)), [level_num, 1], @(x) mean(x, 'omitnan'), NaN)';

%% 图一：左右两个子图
if ismember(1, figure_id)
    figure('Name', 'SHP Level与CV Level的互相关系', 'Position', [100, 100, 1000, 500]);
    
    % 左图
    subplot(1, 2, 1);
    plot(levels, mean_CV_Level, 'b-', 'LineWidth', 1.5, 'DisplayName', 'CV Level 均值');
    xlabel('SHP Level');
    ylabel('CV Level (均值)');
    xlim([1 level_num]); ylim([1 level_num]);
    title('不同 SHP Level 下的 CV Level 均值');
    grid on; legend('Location', 'best');
    
    % 右图
    subplot(1, 2, 2);
    plot(levels, mean_SHP_Level, 'r-', 'LineWidth', 1.5, 'DisplayName', 'SHP Level 均值');
    xlabel('CV Level');
    ylabel('SHP Level (均值)');
    xlim([1 level_num]); ylim([1 level_num]);
    title('不同 CV Level 下的 SHP Level 均值');
    grid on; legend('Location', 'best');
end

%% 图二：BoxCoh估计误差热力图
if ismember(2, figure_id)
    figure('Name', 'BoxCoh估计误差热力图', 'Position', [150, 150, 600, 500]);
    imagesc(levels, levels, error_matrix);
    set(gca, 'YDir', 'normal'); 
    colormap(jet);
    c = colorbar;
    c.Label.String = sprintf('BoxCoh_{%d} - TrueCoh', max_win);
    
    xlabel('SHP Level');
    ylabel('CV Level');
    title(sprintf('SHP Level与CV Level下的BoxCoh\\_%d估计误差热力图', max_win));
    grid off;
    
    % 等值线
    hold on;
    [X, Y] = meshgrid(levels, levels);
    valid_mask = count_matrix >= 5; 
    error_matrix_contour = error_matrix;
    error_matrix_contour(~valid_mask) = NaN;
    contour(X, Y, error_matrix_contour, 8, 'k', 'LineWidth', 0.5);
    hold off;
end

%% 图三：最佳窗口大小平均值热力图
if ismember(3, figure_id)
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
    valid_mask_win = count_matrix >= 5; 
    best_win_contour = best_win_matrix;
    best_win_contour(~valid_mask_win) = NaN;
    contour(X5, Y5, best_win_contour, 8, 'k', 'LineWidth', 0.5);
    hold off;
end

%% 图四：像元数量热力图
if ismember(4, figure_id)
    figure('Name', '像元数量热力图', 'Position', [250, 250, 600, 500]);
    % 可选：使用 log10(count_matrix + 1) 增强对比度
    imagesc(levels, levels, count_matrix);
    set(gca, 'YDir', 'normal'); 
    colormap(jet);
    c = colorbar;
    c.Label.String = '像元数量';
    
    xlabel('SHP Level');
    ylabel('CV Level');
    title('SHP Level与CV Level分级下的像元分布数量热力图');
    grid off;
end
