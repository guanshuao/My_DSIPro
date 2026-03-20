close all; clear; clc;

data_dir = "Stack.mat";
disp('Loading Stack.mat...');
load(data_dir);
disp('Loaded Stack.mat...');

%% 获取基础信息
level_num = Stack.level_num;
levels = (1:level_num)';
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
      
CV_valid = round(Stack.CV_level(isValid));
SHP_valid = round(Stack.SHPNum_level(isValid));

% [行, 列] = [CV Level, SHP Level]
subs2D = double([CV_valid(:), SHP_valid(:)]);

% 1. 计算像元数量矩阵 (对应图四)
count_matrix = accumarray(subs2D, 1, [level_num, level_num], @sum, 0);

% 2. 计算误差均值矩阵 (对应图二)
errors_valid = double(BoxCoh_Max(isValid) - Stack.TrueCoh(isValid));
error_matrix = accumarray(subs2D, errors_valid(:), [level_num, level_num], @(x) mean(x, 'omitnan'), NaN);

% 3. 计算最佳窗口尺寸均值矩阵 (对应图三)
best_win_valid = double(Stack.Best_Win(isValid));
best_win_matrix = accumarray(subs2D, best_win_valid(:), [level_num, level_num], @(x) mean(x, 'omitnan'), NaN);

% 4. 计算一维序列 (对应图一)
mean_CV_Level = accumarray(double(SHP_valid(:)), double(CV_valid(:)), [level_num, 1], @(x) mean(x, 'omitnan'), NaN);
mean_SHP_Level = accumarray(double(CV_valid(:)), double(SHP_valid(:)), [level_num, 1], @(x) mean(x, 'omitnan'), NaN);

%% 写出一维序列数据 (Fig 1)
T1 = table(levels, mean_CV_Level, mean_SHP_Level);
writetable(T1, 'fig1_joint_data.txt', 'Delimiter', '\t');

%% 写出二维热力图数据 (Fig 2, Fig 3, Fig 4)
fid = fopen('joint_heatmap_data.txt', 'w');
fprintf(fid, 'CV_Level\tSHP_Level\tcount\terror\tbest_win\n');
for i = 1:level_num     % CV Level (Y)
    for j = 1:level_num % SHP Level (X)
        fprintf(fid, '%d\t%d\t%f\t%f\t%f\n', i, j, count_matrix(i, j), error_matrix(i, j), best_win_matrix(i, j));
    end
    fprintf(fid, '\n'); % 空行帮助pgfplots分块读取grid
end
fclose(fid);

disp('Data exported to txt successfully (Joint). MATLAB processing complete.');