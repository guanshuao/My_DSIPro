% 寻找最佳窗口大小和最佳估计
clear; clc;
% 检查变量是否存在，否则尝试从 .mat 文件加载
if ~exist('Stack', 'var')
    if exist('Stack.mat', 'file')
        disp('Loading Stack.mat...');
        load('Stack.mat');
    else
        error('Stack variable not in workspace and BoxCoh_stack.mat not found.');
    end
end

% 获取结构体中的所有字段名
fields = fieldnames(Stack);

% 提取 TrueCoh 矩阵
if isfield(Stack, 'TrueCoh')
    TrueCoh = Stack.TrueCoh;
else
    error('TrueCoh field not found in BoxCoh_stack');
end

% 获取维度
[rows, cols] = size(TrueCoh);

% 初始化结果矩阵
Best_Est = zeros(rows, cols);
Best_Win = zeros(rows, cols);

% 初始化一个矩阵用于存储目前发现的最小差异。
% 从无穷大开始，这样第一次比较总是会替换它。
min_diff = inf(rows, cols);

% 遍历结构体中的每个字段
for i = 1:length(fields)
    fname = fields{i};
    
    % 跳过 TrueCoh 字段本身
    if strcmp(fname, 'TrueCoh')
        continue;
    end
    
    % 仅处理以 'BoxCoh_' 开头的字段
    if startsWith(fname, 'BoxCoh_')
        % 从字段名中提取窗口大小 (例如 'BoxCoh_3' -> 3)
        win_str = strrep(fname, 'BoxCoh_', '');
        win_size = str2double(win_str);
        
        % 确保窗口大小有效
        if isnan(win_size)
            warning('Skipping field %s: could not parse window size.', fname);
            continue;
        end
        
        % 获取当前的估计矩阵
        CurrentEst = Stack.(fname);
        
        % 计算与 TrueCoh 的绝对差
        current_diff = abs(CurrentEst - TrueCoh);
        
        % 更新 Best_Est 和 Best_Win 矩阵
        % 我们找出当前差值小于目前发现的最小差值的索引位置
        update_mask = current_diff < min_diff;
        
        % 更新最小差值
        min_diff(update_mask) = current_diff(update_mask);
        
        % 存储最佳估计值
        Best_Est(update_mask) = CurrentEst(update_mask);
        
        % 存储相应的窗口大小
        Best_Win(update_mask) = win_size;
    end
end

disp('计算完成。Best_Est 和 Best_Win 矩阵已生成。');
save('Best_Win.mat', 'Best_Win', '-v7.3');
save('Best_Est.mat', 'Best_Est', '-v7.3');
