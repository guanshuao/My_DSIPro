clear;
load("C:\Users\ags\Downloads\SHP_BWS_15.mat");

% SHP_BWS_15是一个结构体，其中BroNum字段包含了15*15窗口内的BroNum值，其元素取值在1-225之间
% 上下两个子图，上图显示BroNum值的分布，下图将BroNum分为50个level，每个level中包含的像元数尽量相近，展示每个Level内的像元数分布情况

broNumData = SHP_BWS_15.BroNum(:);

figure;
% 上图：显示BroNum值的分布
subplot(2, 1, 1);
histogram(broNumData, 'BinMethod', 'integers');
title('BroNum 值的分布');
xlabel('BroNum 取值 (1-225)');
ylabel('像元数量');

% 下图：将BroNum分为50个level，每个level中包含的像元数尽量相近
% 为确保离散整数划分后像元数绝对一致，采用直接将1-225整数逐个累加分配的方法
subplot(2, 1, 2);
num_levels = 50;

% 先计算 1 到 225 的每个整数的具体频数
hist_counts_int = histcounts(broNumData, 0.5:225.5); 

% 理想状态下每个Level应包含的像元数
target_per_level = numel(broNumData) / num_levels;

level_of_b = zeros(1, 225);
current_level = 1;
current_sum = 0;

for b = 1:225
    level_of_b(b) = current_level;
    current_sum = current_sum + hist_counts_int(b);
    
    % 当当前Level累加数达到目标数时（并且不要超过最大指定的level数），开启新的Level
    if current_sum >= target_per_level && current_level < num_levels
        current_level = current_level + 1;
        current_sum = 0; 
    end
end

actual_num_levels = current_level;
counts = zeros(1, actual_num_levels);

fprintf('\n=== 各 Level 对应的 BroNum 具体值 ===\n');
for i = 1:actual_num_levels
    % 将被分配给该Level对应的整数找出来
    b_in_lvl = find(level_of_b == i);
    
    % 计算该Level下，包含的所有像元总数（严格等价）
    counts(i) = sum(hist_counts_int(b_in_lvl));
    
    % 为了显示更准确，仅提取真正存在像素点的BroNum值打印
    valid_b = b_in_lvl(hist_counts_int(b_in_lvl) > 0);
    
    if isempty(valid_b)
        val_str = '空';
    else
        val_str = sprintf('%d,', valid_b);
        val_str(end) = []; % 删除最后一个逗号
    end
    
    % 打印输出
    fprintf('Level %d: 包含的BroNum有 [%s], 总像元数 = %d\n', i, val_str, counts(i));
end

bar(counts);
title(sprintf('每个Level包含的像元数分布 (按整数重新分配，实际 %d 个Level)', actual_num_levels));
xlabel('Level 编号');
ylabel('像元数量');