% 绘制相干性图的箱线图，使用boxplot函数
% /titan/guanshuao/Kumamoto/DSIPro/BoxCoh 文件夹下存储了40个格式为yymmdd_yymmdd.mat的文件，即为不同时间对的相干性数据
% 每个mat文件体积较大，约2000MB，建议逐个读取处理
% 同时还有其他一些.mat文件，与相干性无关，读取数据时需要注意
% 每个.mat文件中每个元素均为复数，取其模值即为相干性值
% 绘制箱线图，横轴为时间对（保持文件名顺序和名称）


clear; clc;

% 数据文件夹路径
dataDir = '/titan/guanshuao/Kumamoto/Data/AdpCoh_BWS_21';

% 获取该目录下所有文件
dirOutput = dir(fullfile(dataDir, '*.mat'));
fileNames = {dirOutput.name};

% 筛选符合 yymmdd_yymmdd.mat 格式的文件
pattern = '^\d{6}_\d{6}\.mat$';
matchIdx = ~cellfun(@isempty, regexp(fileNames, pattern));
targetFiles = fileNames(matchIdx);

% 排序文件名 (按时间对顺序)
targetFiles = sort(targetFiles);

% 检查是否找到文件
if isempty(targetFiles)
    error('未找到符合格式的文件');
end

% 初始化绘图
figure('Visible', 'off'); % 无头模式下不显示窗口
hold on;
title('Coherence Boxplot');
xlabel('Time Pair');
ylabel('Coherence Value');

numFiles = length(targetFiles);

% 逐个处理文件
for i = 1:numFiles
    fileName = targetFiles{i};
    filePath = fullfile(dataDir, fileName);
    
    fprintf('正在处理 (%d/%d): %s ...\n', i, numFiles, fileName);
    
    % 加载数据
    try
        dataStruct = load(filePath);
        vars = fieldnames(dataStruct);
        if isempty(vars)
            warning('文件 %s 为空，跳过', fileName);
            continue;
        end
        
        % 假设数据在第一个变量中
        rawData = dataStruct.(vars{1});
        
        % 计算模值 (相干性)
        cohData = abs(rawData(:));
        
        % 绘制箱线图
        % 由于数据量巨大，建议不绘制异常值点 ('Symbol', '') 以避免卡顿
        boxplot(cohData, 'Positions', i, 'Widths', 0.6, 'Symbol', '');
        
        % 清理内存
        clear dataStruct rawData cohData;
        
    catch ME
        warning('处理文件 %s 时出错: %s', fileName, ME.message);
    end
end

% 设置坐标轴
set(gca, 'XTick', 1:numFiles);
% 生成标签 (去掉 .mat 后缀)
labels = strrep(targetFiles, '.mat', '');
set(gca, 'XTickLabel', labels);
xtickangle(90); % 旋转标签
xlim([0, numFiles + 1]);
grid on;
hold off;

% 保存图片
outputFile = 'coherence_boxplot_bws21.png';
fprintf('正在保存图片到 %s ...\n', outputFile);
saveas(gcf, outputFile);

fprintf('处理完成。\n');