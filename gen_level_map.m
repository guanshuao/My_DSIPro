function level_map = gen_level_map(map, type, level_num)
% gen_level_map: 将输入矩阵划分为若干个等级，输出等级矩阵
% 输入：
%   map : 原始二维矩阵
%   type : 'float' (浮点) 或 'integer' (整数)
%   level_num : 划分等级数，默认为 50
% 输出：
%   level_map : 与输入同样大小的矩阵，元素值为 1 到 level_num 对应的等级数

    if nargin < 3
        level_num = 50;
    end
    
    level_map = zeros(size(map), 'uint16');
    data = map(:);
    
    % 剔除可能存在的NaN值
    valid_mask = ~isnan(data);
    valid_data = data(valid_mask);
    
    if isempty(valid_data)
        warning('输入矩阵为空或全为NaN');
        return;
    end

    figure('Name', 'Histogram and Levels', 'Position', [100, 100, 800, 600]);
    
    if strcmpi(type, 'float') || strcmpi(type, 'f')
        % --- 浮点数分级算法 ---
        
        % 上图：显示值的分布
        subplot(2, 1, 1);
        histogram(valid_data, 100);
        title('原始数值的分布 (Float)');
        xlabel('原始取值');
        ylabel('像元数量');

        % 计算分位数区间
        pct_edges = prctile(valid_data, linspace(0, 100, level_num + 1));
        edges = unique(pct_edges); % 去重
        
        % 扩展边界以确保包含数据范围的两端
        edges(1) = edges(1) - eps; 
        edges(end) = edges(end) + eps;

        counts = histcounts(valid_data, edges);

        % 下图：Level包含的像元数分布
        subplot(2, 1, 2);
        bar(counts);
        actual_num_levels = length(counts);
        title(sprintf('每个Level包含的像元数分布 (实际划分为 %d 个Level)', actual_num_levels));
        xlabel('Level 编号');
        ylabel('像元数量');

        % 输出信息并在原图中映射
        fprintf('\n=== 各 Level 对应的 Float 数值区间 ===\n');
        valid_levels = discretize(valid_data, edges);
        level_map(valid_mask) = valid_levels;
        
        for i = 1:length(counts)
            if i == length(counts)
                fprintf('Level %d: 区间 [%.4f, %.4f], 总像元数 = %d\n', i, edges(i), edges(i+1), counts(i));
            else
                fprintf('Level %d: 区间 [%.4f, %.4f), 总像元数 = %d\n', i, edges(i), edges(i+1), counts(i));
            end
        end

    elseif strcmpi(type, 'integer') || strcmpi(type, 'int') || strcmpi(type, 'i')
        % --- 整数分级算法 ---
        
        min_val = floor(min(valid_data));
        max_val = ceil(max(valid_data));
        
        % 上图：显示值的分布
        subplot(2, 1, 1);
        histogram(valid_data, 'BinMethod', 'integers');
        title('原始数值的分布 (Integer)');
        xlabel(sprintf('原始取值 (%d - %d)', min_val, max_val));
        ylabel('像元数量');

        % 获取每个整数的具体频数
        hist_counts_int = histcounts(valid_data, (min_val-0.5):(max_val+0.5)); 

        val_range = min_val:max_val;
        level_of_val = zeros(size(val_range));
        
        % 使用全局累加，而不是局部清零
        cum_counts = cumsum(hist_counts_int);
        
        for k = 1:length(val_range)
            % 根据当前的累积像素占比，映射到 1~level_num
            ratio = cum_counts(k) / cum_counts(end);
            % 使用 ceil 映射分布比例，为了防止 ratio=0 变为 level 0，用 max 修正
            mapped_lvl = max(1, ceil(ratio * level_num));
            level_of_val(k) = mapped_lvl;
        end

        counts = zeros(1, level_num);

        fprintf('\n=== 各 Level 对应的 Integer 整数值 ===\n');
        for i = 1:level_num
            idx_in_lvl = find(level_of_val == i);
            counts(i) = sum(hist_counts_int(idx_in_lvl));
            
            % 获取真正存在的取值
            valid_idx = idx_in_lvl(hist_counts_int(idx_in_lvl) > 0);
            
            if isempty(valid_idx)
                val_str = '空';
            else
                val_str = sprintf('%d,', val_range(valid_idx));
                val_str(end) = []; % 删除最后一个逗号
            end
            
            fprintf('Level %d: 包含的整数有 [%s], 总像元数 = %d\n', i, val_str, counts(i));
        end

        % 下图：Level包含的像元数分布
        subplot(2, 1, 2);
        bar(counts);
        title(sprintf('每个Level包含的像元数分布 (%d 个Level)', level_num));
        xlabel('Level 编号');
        ylabel('像元数量');

        % 映射到输出图上
        for k = 1:length(val_range)
            if hist_counts_int(k) > 0
                level_map(map == val_range(k)) = level_of_val(k);
            end
        end

    else
        error('不支持的 type 类型，请输入 ''float'' 或 ''integer''.');
    end
    
    % --- 绘制图2：将map和level_map使用imagesc加上colorbar处理 --- %
    figure('Name', 'Spatial Maps', 'Position', [150, 150, 900, 450]);
    
    subplot(1, 2, 1);
    imagesc(map);
    colorbar;
    title('Original Map');
    
    subplot(1, 2, 2);
    imagesc(level_map);
    colorbar;
    title('Level Map');
    
end
