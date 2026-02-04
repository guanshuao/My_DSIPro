function [Coh,Ph] = AdpCohEst_LowMem(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode,NumWorkers)
%AdpCohEst_LowMem 低内存版本的自适应相干性估计函数
%   本函数使用SHP（同质像元）估计相干性和干涉相位，支持N视图像。
%   针对 average 模式进行了内存优化和并行加速。
%
%   用法:
%       [Coh,Ph] = AdpCohEst_LowMem(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode,NumWorkers)
%   或:
%       [Coh] = AdpCohEst_LowMem(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode,NumWorkers)
%
%   输入参数:
%   - mlistack: 高度 × 宽度 × 页数 的（实数）矩阵，例如 SAR 单视强度序列
%   - mlilist:  n×1 的文件列表，包含以 <yyyymmdd> 命名的强度图像
%   - infstack: 高度 × 宽度 × 页数 的（复数）矩阵，例如去除相位梯度后的干涉图
%   - inflist:  n×2 的文件列表，包含以 <yyyymmdd yyyymmdd> 命名的干涉图
%   - SHP:      同质像元结构体，参见 SHP_SelPoint.m
%   - BiasCorr: 是否进行偏差校正，'y' 或 'n'（默认 'n'）
%   - OutputMode: 输出模式，'stack'（保留全部）或 'average'（时间平均，默认）
%   - NumWorkers: 并行计算核心数，默认为可用CPU核心数（仅 average 模式有效）
%
%   输出参数:
%   - Coh: 相干性幅度（average模式为二维，stack模式为三维）
%   - Ph:  干涉相位（仅 stack 模式支持）
%
%   优化说明:
%   - average 模式：采用累加策略，不存储完整相干性堆栈，大幅减少内存占用
%   - average 模式：支持按行分块并行计算，充分利用多核CPU
%   - stack 模式：保持原有算法不变
%
%   参考文献:
%   [1] Delineation of built-up land change from SAR stack by analysing 
%       the coefficient of variation. M Jiang et al., ISPRS JPRS, 2020.
%   [2] Hybrid Approach for Unbiased Coherence Estimation for Multitemporal InSAR.
%       Mi Jiang et al., IEEE TGRS, 2014.
%
%   Mi JIANG, Sun Yat-sen University
%   ======================================================================
%   02/2026 低内存版本：针对 average 模式优化内存和并行计算
%   ======================================================================

%% 参数默认值设置
if nargin < 8 || isempty(NumWorkers)
    NumWorkers = feature('numcores');
end

if nargin < 7 || isempty(OutputMode)
    OutputMode = 'average';
end

if nargin < 6 || isempty(BiasCorr)
    BiasCorr = 'n';
end

if nargin < 5
    help AdpCohEst_LowMem
    return
end

%% 根据模式选择处理方式
if strcmpi(OutputMode, 'stack')
    % stack 模式：使用原有算法
    [Coh, Ph] = AdpCohEst_Stack(mlistack, mlilist, infstack, inflist, SHP, BiasCorr);
else
    % average 模式：使用低内存并行算法
    if nargout > 1
        warning('average 模式下只输出 Coh，Ph 输出将被忽略。');
    end
    Coh = AdpCohEst_Average_LowMem(mlistack, mlilist, infstack, inflist, SHP, BiasCorr, NumWorkers);
    Ph = [];
end

end

%% ========================================================================
%  Stack 模式：保持原有算法
%  ========================================================================
function [Coh, Ph] = AdpCohEst_Stack(mlistack, mlilist, infstack, inflist, SHP, BiasCorr)

tic;
[nlines, nwidths, npages] = size(infstack);
Coh = zeros(nlines, nwidths, npages, 'single');
[~, idx] = ismember(inflist, mlilist);

CalWin = SHP.CalWin;
RadiusRow = (CalWin(1)-1)/2;
RadiusCol = (CalWin(2)-1)/2;

for ii = 1:npages
    m1 = mlistack(:,:,idx(ii,1));
    m2 = mlistack(:,:,idx(ii,2));
    Intf = sqrt(m1.*m2) .* infstack(:,:,ii);
    
    m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
    m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
    Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');
    
    nu = zeros(nlines, nwidths, 'single');
    de1 = nu;
    de2 = nu;
    num = 1;
    
    for jj = 1:nwidths
        for kk = 1:nlines
            x_global = jj + RadiusCol;
            y_global = kk + RadiusRow;
            MasterValue = m1(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            SlaveValue = m2(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            InterfValue = Intf(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            MasterValue = MasterValue(SHP.PixelInd(:,num));
            SlaveValue = SlaveValue(SHP.PixelInd(:,num));
            InterfValue = InterfValue(SHP.PixelInd(:,num));
            nu(kk,jj) = mean(InterfValue, 'omitnan');
            de1(kk,jj) = mean(MasterValue, 'omitnan');
            de2(kk,jj) = mean(SlaveValue, 'omitnan');
            num = num + 1;
        end
    end
    
    denominator = sqrt(de1.*de2);
    denominator(denominator == 0) = NaN;
    Coh(:,:,ii) = nu ./ denominator;
    fprintf('adp. coherence: %3d / %d is finished...\n', ii, npages);
end

Coh(isnan(Coh) | isinf(Coh)) = 0;
CohComplex = Coh;
Coh = abs(Coh);

% 偏差校正
if strcmpi(BiasCorr, 'y')
    tmp = padarray(Coh, [RadiusRow RadiusCol], 'symmetric');
    for ii = 1:npages
        num = 1;
        for jj = 1:nwidths
            for kk = 1:nlines
                x_global = jj + RadiusCol;
                y_global = kk + RadiusRow;
                CohValue = tmp(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol, ii);
                CohValue = CohValue(SHP.PixelInd(:,num));
                CohValue(CohValue <= 0) = NaN;
                Coh(kk,jj,ii) = mean(log(CohValue), 'omitnan');
                num = num + 1;
            end
        end
        fprintf('bias mitigation: %3d / %d is finished...\n', ii, npages);
    end
    Coh = exp(Coh);
    Coh(isnan(Coh) | isinf(Coh)) = 0;
end

Ph = angle(CohComplex);

t = toc;
disp(['AdpCohEst (stack mode) completed in ', int2str(t/60), ' min(s).']);
disp('Done!');

end

%% ========================================================================
%  Average 模式：低内存并行算法
%  ========================================================================
function Coh = AdpCohEst_Average_LowMem(mlistack, mlilist, infstack, inflist, SHP, BiasCorr, NumWorkers)

tic;
[nlines, nwidths, npages] = size(infstack);
[~, idx] = ismember(inflist, mlilist);

CalWin = SHP.CalWin;
RadiusRow = (CalWin(1)-1)/2;
RadiusCol = (CalWin(2)-1)/2;

%% 确保并行池已启动
pool = gcp('nocreate');
if isempty(pool)
    disp(['正在启动并行池（', num2str(NumWorkers), ' 个工作进程）...']);
    parpool(NumWorkers);
elseif pool.NumWorkers ~= NumWorkers
    disp(['调整并行池大小为 ', num2str(NumWorkers), ' 个工作进程...']);
    delete(pool);
    parpool(NumWorkers);
else
    disp(['并行池已就绪，工作进程数: ', num2str(pool.NumWorkers)]);
end

%% 计算分块参数（按行分块）
NumBlocks = NumWorkers;
validRowsPerBlock = ceil(nlines / NumBlocks);

blockInfo = struct('startRow', {}, 'endRow', {}, ...
                   'validStartRow', {}, 'validEndRow', {}, ...
                   'globalStartRow', {}, 'globalEndRow', {});

for b = 1:NumBlocks
    validStart = (b-1) * validRowsPerBlock + 1;
    validEnd = min(b * validRowsPerBlock, nlines);
    
    if validStart > nlines
        break;
    end
    
    % 实际读取区域（包含重叠边界）
    actualStart = max(1, validStart - RadiusRow);
    actualEnd = min(nlines, validEnd + RadiusRow);
    
    % 有效区域在分块内的相对位置
    relativeValidStart = validStart - actualStart + 1;
    relativeValidEnd = validEnd - actualStart + 1;
    
    blockInfo(b).startRow = actualStart;
    blockInfo(b).endRow = actualEnd;
    blockInfo(b).validStartRow = relativeValidStart;
    blockInfo(b).validEndRow = relativeValidEnd;
    blockInfo(b).globalStartRow = validStart;
    blockInfo(b).globalEndRow = validEnd;
end

actualNumBlocks = length(blockInfo);
disp(['实际分块数: ', num2str(actualNumBlocks)]);

%% 提取各分块的 SHP 索引
% SHP.PixelInd 是 (CalWin²) × (nlines*nwidths) 的矩阵
% 需要为每个分块提取对应列的 SHP 索引
blockSHP = cell(actualNumBlocks, 1);
for b = 1:actualNumBlocks
    globalStartRow = blockInfo(b).globalStartRow;
    globalEndRow = blockInfo(b).globalEndRow;
    numValidRows = globalEndRow - globalStartRow + 1;
    
    % 提取该分块对应的 SHP 索引列
    colIndices = zeros(1, numValidRows * nwidths);
    cnt = 0;
    for col = 1:nwidths
        for row = globalStartRow:globalEndRow
            cnt = cnt + 1;
            % 全局像素索引（列优先）
            globalIdx = (col-1)*nlines + row;
            colIndices(cnt) = globalIdx;
        end
    end
    blockSHP{b} = SHP.PixelInd(:, colIndices);
end

%% 准备各分块的强度数据
disp('准备分块数据...');
blockMLI = cell(actualNumBlocks, 1);
for b = 1:actualNumBlocks
    blockMLI{b} = mlistack(blockInfo(b).startRow:blockInfo(b).endRow, :, :);
end
clear mlistack;

%% 逐干涉对处理，累加相干性
CohSum = zeros(nlines, nwidths, 'single');
CohCount = zeros(nlines, nwidths, 'single');

for ii = 1:npages
    % 提取当前干涉对的数据
    currentIntf = infstack(:,:,ii);
    m1_idx = idx(ii, 1);
    m2_idx = idx(ii, 2);
    
    % 准备各分块的干涉图数据
    blockIntf = cell(actualNumBlocks, 1);
    for b = 1:actualNumBlocks
        blockIntf{b} = currentIntf(blockInfo(b).startRow:blockInfo(b).endRow, :);
    end
    
    % 并行处理各分块
    blockResults = cell(actualNumBlocks, 1);
    
    parfor b = 1:actualNumBlocks
        blockResults{b} = processBlockAdp(blockMLI{b}, blockIntf{b}, blockSHP{b}, ...
                                           m1_idx, m2_idx, CalWin, blockInfo(b), nwidths);
    end
    
    % 合并结果到累加器
    for b = 1:actualNumBlocks
        globalStartRow = blockInfo(b).globalStartRow;
        globalEndRow = blockInfo(b).globalEndRow;
        
        cohBlock = blockResults{b};
        validMask = ~isnan(cohBlock) & ~isinf(cohBlock);
        cohBlock(~validMask) = 0;
        
        CohSum(globalStartRow:globalEndRow, :) = CohSum(globalStartRow:globalEndRow, :) + cohBlock;
        CohCount(globalStartRow:globalEndRow, :) = CohCount(globalStartRow:globalEndRow, :) + single(validMask);
    end
    
    % 清理当前干涉对数据
    clear blockIntf blockResults;
    
    % 进度显示
    if mod(ii, 10) == 0 || ii == npages
        fprintf('adp. coherence: %3d / %d completed\n', ii, npages);
    end
end

clear infstack blockMLI;

%% 计算平均相干性
CohCount(CohCount == 0) = NaN;
Coh = CohSum ./ CohCount;
Coh(isnan(Coh) | isinf(Coh)) = 0;

%% 偏差校正（如果需要）
if strcmpi(BiasCorr, 'y')
    disp('正在进行偏差校正...');
    Coh_padded = padarray(Coh, [RadiusRow RadiusCol], 'symmetric');
    Coh_corrected = zeros(nlines, nwidths, 'single');
    
    num = 1;
    for jj = 1:nwidths
        for kk = 1:nlines
            x_global = jj + RadiusCol;
            y_global = kk + RadiusRow;
            CohValue = Coh_padded(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            CohValue = CohValue(SHP.PixelInd(:,num));
            CohValue(CohValue <= 0) = NaN;
            Coh_corrected(kk,jj) = mean(log(CohValue), 'omitnan');
            num = num + 1;
        end
    end
    
    Coh = exp(Coh_corrected);
    Coh(isnan(Coh) | isinf(Coh)) = 0;
end

t = toc;
disp(['AdpCohEst_LowMem (average mode) completed in ', num2str(t/60, '%.2f'), ' min(s).']);
disp('Done!');

end

%% ========================================================================
%  处理单个分块的自适应相干性计算
%  ========================================================================
function cohBlock = processBlockAdp(blockMLI, blockIntf, blockSHP, m1_idx, m2_idx, CalWin, blockInfo, nwidths)

RadiusRow = (CalWin(1)-1)/2;
RadiusCol = (CalWin(2)-1)/2;

% 提取主辅影像强度
m1 = blockMLI(:,:,m1_idx);
m2 = blockMLI(:,:,m2_idx);

% 计算干涉项
Intf = sqrt(m1.*m2) .* blockIntf;

% 边缘填充
m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');

% 有效处理区域
validStartRow = blockInfo.validStartRow;
validEndRow = blockInfo.validEndRow;
numValidRows = validEndRow - validStartRow + 1;

% 预分配结果
cohBlock = zeros(numValidRows, nwidths, 'single');

% 处理每个像素
pixelCount = 0;
for jj = 1:nwidths
    x_global = jj + RadiusCol;
    for kk = validStartRow:validEndRow
        y_global = kk + RadiusRow;
        pixelCount = pixelCount + 1;
        
        MasterValue = m1(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
        SlaveValue = m2(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
        InterfValue = Intf(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
        
        MasterValue = MasterValue(blockSHP(:,pixelCount));
        SlaveValue = SlaveValue(blockSHP(:,pixelCount));
        InterfValue = InterfValue(blockSHP(:,pixelCount));
        
        nu = mean(InterfValue, 'omitnan');
        de1 = mean(MasterValue, 'omitnan');
        de2 = mean(SlaveValue, 'omitnan');
        
        denominator = sqrt(de1 * de2);
        if denominator == 0
            cohBlock(kk - validStartRow + 1, jj) = NaN;
        else
            cohBlock(kk - validStartRow + 1, jj) = abs(nu / denominator);
        end
    end
end

end
