clear;
load("C:\Users\ags\Downloads\CV_15.mat");
% CV_15是一个矩阵，其中每个元素都是正浮点数

cvData = CV_15(:);
% 剔除可能存在的NaN值
cvData(isnan(cvData)) = [];

figure;
% 上图：显示CV值的分布
subplot(2, 1, 1);
histogram(cvData, 100); % 对于浮点数，使用100个bin展示连续分布
title('CV 值的分布');
xlabel('CV 取值');
ylabel('像元数量');

% 下图：将CV分为50个level，每个level中包含的像元数尽量相近
subplot(2, 1, 2);
num_levels = 50;

% 由于CV是连续的浮点数，采用分位数（百分位数）进行区间划分是理想且直接的方法
pct_edges = prctile(cvData, linspace(0, 100, num_levels + 1));

% 去重以防万一（存在极多相同浮点值时的小概率事件）
edges = unique(pct_edges);

% 扩展边界以确保包含数据范围的两端
edges(1) = edges(1) - eps; 
edges(end) = edges(end) + eps;

counts = histcounts(cvData, edges);

bar(counts);
title(sprintf('每个Level包含的像元数分布 (实际划分为 %d 个Level)', length(counts)));
xlabel('Level 编号');
ylabel('像元数量');

% 输出每个Level的CV区间跨度和包含的像元数
fprintf('\n=== 各 Level 对应的 CV 浮点值区间 ===\n');
for i = 1:length(counts)
    % 打印对应的浮点区间值。histcounts 默认左闭右开，最后一个bin是全闭区间
    if i == length(counts)
        fprintf('Level %d: 区间 [%.4f, %.4f], 总像元数 = %d\n', i, edges(i), edges(i+1), counts(i));
    else
        fprintf('Level %d: 区间 [%.4f, %.4f), 总像元数 = %d\n', i, edges(i), edges(i+1), counts(i));
    end
end