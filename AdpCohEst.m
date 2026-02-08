function [Coh,Ph] = AdpCohEst(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode,NumWorkers)
%AdpCohEst 并行优化版本的自适应相干性估计函数
%   本函数使用空间分块并行策略估计相干性和干涉相位，支持SHP（同质像元）
%   针对大数据量进行了内存和速度优化
%
%   用法:
%       [Coh,Ph] = AdpCohEst(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode,NumWorkers)
%   或:
%       [Coh] = AdpCohEst(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode,NumWorkers)
%
%   输入参数:
%   - mlistack:   高度 × 宽度 × 页数 的实数矩阵，SAR单视强度图序列
%   - mlilist:    n×1 文件列表，强度图以 <yyyymmdd> 命名
%   - infstack:   高度 × 宽度 × 页数 的复数矩阵，去除相位梯度后的单视干涉图
%   - inflist:    n×2 文件列表，复数干涉图以 <yyyymmdd yyyymmdd> 命名
%   - SHP:        同质像元结构体，详见 SHP_SelPoint.m
%   - BiasCorr:   是否启用对数矩估计进行偏差校正，'y' 或 'n'，默认 'n'
%   - OutputMode: 输出模式，'stack'（完整堆栈）或 'average'（平均值），默认 'stack'
%   - NumWorkers: 并行工作进程数，默认为 CPU 核心数；设为 1 时使用串行模式
%
%   输出参数:
%   - Coh:        相干性幅度（stack模式为3D矩阵，average模式为2D矩阵）
%   - Ph:         干涉相位（仅stack模式支持）
%
%   参考文献:
%   [1] Delineation of built-up land change from SAR stack by analysing the coefficient of variation
%        M Jiang, A Hooper, X Tian, J Xu, SN Chen, ZF Ma, X Cheng,
%        ISPRS Journal of Photogrammetry and Remote Sensing 169, Nov. 2020, pp. 93-108.
%
%   [2] Hybrid Approach for Unbiased Coherence Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding and Zhiwei Li,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 52, no. 5, May 2014, pp. 2459-2473.
%
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University
%   并行优化版本
%   ======================================================================
%   02/2026 空间并行版本：使用分块策略进行并行优化，内存优化
%   ======================================================================

%% 参数默认值设置
if nargin < 8 || isempty(NumWorkers)
    NumWorkers = 1;
end

if nargin < 7 || isempty(OutputMode)
    OutputMode = 'stack';
end

if nargin < 6 || isempty(BiasCorr)
    BiasCorr = 'n';
end

if nargin < 5
    help AdpCohEst
    return
end

%% 输入验证与初始化
tic;
[nlines, nwidths, npages] = size(infstack);
[~, idx] = ismember(inflist, mlilist);

CalWin = SHP.CalWin;
RadiusRow = (CalWin(1)-1)/2;
RadiusCol = (CalWin(2)-1)/2;

%% 根据并行数选择处理模式
if NumWorkers == 1
    % 串行模式：调用原始算法
    disp('使用串行模式...');
    if strcmpi(OutputMode, 'stack')
        [Coh, Ph] = AdpCohEst_Serial_Stack(mlistack, infstack, idx, SHP, BiasCorr, ...
                                           nlines, nwidths, npages, RadiusRow, RadiusCol);
    else
        Coh = AdpCohEst_Serial_Average(mlistack, infstack, idx, SHP, BiasCorr, ...
                                        nlines, nwidths, npages, RadiusRow, RadiusCol);
        if nargout > 1
            warning('average 模式下只输出 Coh，Ph 输出将被忽略。');
            Ph = [];
        end
    end
else
    % 并行模式：使用空间分块
    disp(['使用并行模式（', num2str(NumWorkers), ' 个工作进程）...']);
    
    % 确保并行池已启动
    pool = gcp('nocreate');
    if isempty(pool)
        disp('正在启动并行池...');
        parpool(NumWorkers);
    elseif pool.NumWorkers ~= NumWorkers
        disp(['当前并行池有 ', num2str(pool.NumWorkers), ' 个工作进程，重新创建...']);
        delete(pool);
        parpool(NumWorkers);
    else
        disp(['并行池已就绪，工作进程数: ', num2str(pool.NumWorkers)]);
    end
    
    if strcmpi(OutputMode, 'stack')
        [Coh, Ph] = AdpCohEst_Parallel_Stack(mlistack, infstack, idx, SHP, BiasCorr, ...
                                              nlines, nwidths, npages, RadiusRow, RadiusCol, NumWorkers);
    else
        Coh = AdpCohEst_Parallel_Average(mlistack, infstack, idx, SHP, BiasCorr, ...
                                          nlines, nwidths, npages, RadiusRow, RadiusCol, NumWorkers);
        if nargout > 1
            warning('average 模式下只输出 Coh，Ph 输出将被忽略。');
            Ph = [];
        end
    end
end

%% 完成
t = toc;
disp(['AdpCohEst 运算完成，耗时 ', num2str(t/60, '%.2f'), ' 分钟']);
disp('完成！');

end

%% ========================================================================
%  串行模式 - Stack 输出
%  ========================================================================
function [Coh, Ph] = AdpCohEst_Serial_Stack(mlistack, infstack, idx, SHP, BiasCorr, ...
                                             nlines, nwidths, npages, RadiusRow, RadiusCol)
    
    Coh = zeros(nlines, nwidths, npages, 'single');
    
    for ii = 1:npages
        m1 = mlistack(:,:,idx(ii,1));
        m2 = mlistack(:,:,idx(ii,2));
        Intf = sqrt(m1.*m2) .* infstack(:,:,ii);
        
        m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
        m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
        Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');
        
        nu = zeros(nlines, nwidths, 'single');
        de1 = nu; de2 = nu;
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
                nu(kk,jj) = mean(InterfValue);
                de1(kk,jj) = mean(MasterValue);
                de2(kk,jj) = mean(SlaveValue);
                num = num + 1;
            end
        end
        
        Coh(:,:,ii) = nu ./ sqrt(de1.*de2);
        
        fprintf('adp. coherence: %3d / %d is finished...\n', ii, npages);
    end
    
    Coh(isnan(Coh)) = 0;
    
    if nargout > 1
        Ph = angle(Coh);
    end
    
    Coh = abs(Coh);
    
    % Bias correction
    if strcmpi(BiasCorr, 'y')
        Coh = applyBiasCorrection_Stack(Coh, SHP, nlines, nwidths, npages, RadiusRow, RadiusCol);
    end
end

%% ========================================================================
%  串行模式 - Average 输出
%  ========================================================================
function Coh = AdpCohEst_Serial_Average(mlistack, infstack, idx, SHP, BiasCorr, ...
                                         nlines, nwidths, npages, RadiusRow, RadiusCol)
    
    CohSum = zeros(nlines, nwidths, 'single');
    
    for ii = 1:npages
        m1 = mlistack(:,:,idx(ii,1));
        m2 = mlistack(:,:,idx(ii,2));
        Intf = sqrt(m1.*m2) .* infstack(:,:,ii);
        
        m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
        m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
        Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');
        
        nu = zeros(nlines, nwidths, 'single');
        de1 = nu; de2 = nu;
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
                nu(kk,jj) = mean(InterfValue);
                de1(kk,jj) = mean(MasterValue);
                de2(kk,jj) = mean(SlaveValue);
                num = num + 1;
            end
        end
        
        CohSingle = nu ./ sqrt(de1.*de2);
        CohSingle(isnan(CohSingle)) = 0;
        CohSingle = abs(CohSingle);
        
        % Bias correction for current page
        if strcmpi(BiasCorr, 'y')
            CohSingle = applyBiasCorrection_Single(CohSingle, SHP, nlines, nwidths, RadiusRow, RadiusCol);
        end
        
        CohSum = CohSum + CohSingle;
        
        fprintf('adp. coherence: %3d / %d is finished...\n', ii, npages);
    end
    
    Coh = CohSum / npages;
end

%% ========================================================================
%  并行模式 - Stack 输出（空间分块）
%  ========================================================================
function [Coh, Ph] = AdpCohEst_Parallel_Stack(mlistack, infstack, idx, SHP, BiasCorr, ...
                                               nlines, nwidths, npages, RadiusRow, RadiusCol, NumWorkers)
    
    %% 计算分块参数（按行分块）
    blockInfo = computeBlockInfo(nlines, NumWorkers, RadiusRow);
    actualNumBlocks = length(blockInfo);
    disp(['实际分块数: ', num2str(actualNumBlocks)]);
    
    %% 准备各分块的数据
    blockData_mli = cell(actualNumBlocks, 1);
    blockData_inf = cell(actualNumBlocks, 1);
    
    for b = 1:actualNumBlocks
        startRow = blockInfo(b).startRow;
        endRow = blockInfo(b).endRow;
        blockData_mli{b} = mlistack(startRow:endRow, :, :);
        blockData_inf{b} = infstack(startRow:endRow, :, :);
    end
    
    % 提取各分块的SHP索引
    blockData_SHP = extractBlockSHP(SHP, blockInfo, nlines, nwidths);
    
    %% 并行处理各分块
    blockResults = cell(actualNumBlocks, 1);
    
    disp('开始并行计算...');
    parfor b = 1:actualNumBlocks
        fprintf('分块 %d/%d 开始处理...\n', b, actualNumBlocks);
        
        % 处理当前分块
        blockResults{b} = processBlock_Stack(blockData_mli{b}, blockData_inf{b}, ...
                                              blockData_SHP{b}, idx, npages, ...
                                              RadiusRow, RadiusCol, blockInfo(b), nwidths);
        
        fprintf('分块 %d/%d 处理完成\n', b, actualNumBlocks);
    end
    
    % 清理分块数据（parfor后不再需要）
    clear blockData_mli blockData_inf blockData_SHP;
    
    %% 合并结果
    disp('正在合并分块结果...');
    Coh = zeros(nlines, nwidths, npages, 'single');
    Ph = zeros(nlines, nwidths, npages, 'single');
    
    for b = 1:actualNumBlocks
        globalStartRow = blockInfo(b).globalStartRow;
        globalEndRow = blockInfo(b).globalEndRow;
        Coh(globalStartRow:globalEndRow, :, :) = blockResults{b}.Coh;
        Ph(globalStartRow:globalEndRow, :, :) = blockResults{b}.Ph;
    end
    
    %% Bias correction
    if strcmpi(BiasCorr, 'y')
        disp('正在进行偏差校正...');
        Coh = applyBiasCorrection_Stack_Parallel(Coh, SHP, nlines, nwidths, npages, ...
                                                  RadiusRow, RadiusCol, NumWorkers);
    end
end

%% ========================================================================
%  并行模式 - Average 输出（空间分块）
%  按页处理流程：分块计算相干性 -> 合并为整张 -> BiasCorr -> 累加到sum -> 除以页数
%  ========================================================================
function Coh = AdpCohEst_Parallel_Average(mlistack, infstack, idx, SHP, BiasCorr, ...
                                           nlines, nwidths, npages, RadiusRow, RadiusCol, NumWorkers)
    
    %% 计算分块参数
    blockInfo = computeBlockInfo(nlines, NumWorkers, RadiusRow);
    actualNumBlocks = length(blockInfo);
    disp(['实际分块数: ', num2str(actualNumBlocks)]);
    
    %% 预先准备各分块的 SHP 索引（只需计算一次）
    blockSHP_all = extractBlockSHP(SHP, blockInfo, nlines, nwidths);
    
    %% 初始化累加变量
    CohSum = zeros(nlines, nwidths, 'single');
    
    %% 按页处理
    for ii = 1:npages
        %% 1. 准备当前页的分块数据
        blockData_mli = cell(actualNumBlocks, 1);
        blockData_inf = cell(actualNumBlocks, 1);
        
        for b = 1:actualNumBlocks
            startRow = blockInfo(b).startRow;
            endRow = blockInfo(b).endRow;
            
            % 只读取当前页需要的两个 mli 和一个 inf
            blockData_mli{b}.m1 = mlistack(startRow:endRow, :, idx(ii,1));
            blockData_mli{b}.m2 = mlistack(startRow:endRow, :, idx(ii,2));
            blockData_inf{b} = infstack(startRow:endRow, :, ii);
        end
        
        %% 2. 并行计算各分块的相干性
        blockResults = cell(actualNumBlocks, 1);
        
        parfor b = 1:actualNumBlocks
            blockResults{b} = processBlock_Average_SinglePage(blockData_mli{b}, blockData_inf{b}, ...
                                                               blockSHP_all{b}, RadiusRow, RadiusCol, ...
                                                               blockInfo(b), nwidths);
        end
        
        %% 3. 合并为整张相干性
        CohSingle = zeros(nlines, nwidths, 'single');
        for b = 1:actualNumBlocks
            globalStartRow = blockInfo(b).globalStartRow;
            globalEndRow = blockInfo(b).globalEndRow;
            CohSingle(globalStartRow:globalEndRow, :) = blockResults{b};
        end
        
        %% 4. 进行偏差矫正（使用完整的相干性矩阵和完整的 SHP）
        if strcmpi(BiasCorr, 'y')
            CohSingle = applyBiasCorrection_Single(CohSingle, SHP, nlines, nwidths, RadiusRow, RadiusCol);
        end
        
        %% 5. 累加到 sum 中间变量
        CohSum = CohSum + CohSingle;
        
        fprintf('adp. coherence (parallel): %3d / %d is finished...\n', ii, npages);
    end
    
    Coh = CohSum / npages;
end

%% ========================================================================
%  计算分块信息
%  ========================================================================
function blockInfo = computeBlockInfo(nlines, NumWorkers, RadiusRow)
    validRowsPerBlock = ceil(nlines / NumWorkers);
    blockInfo = struct('startRow', {}, 'endRow', {}, ...
                       'validStartRow', {}, 'validEndRow', {}, ...
                       'globalStartRow', {}, 'globalEndRow', {});
    
    for b = 1:NumWorkers
        validStart = (b-1) * validRowsPerBlock + 1;
        validEnd = min(b * validRowsPerBlock, nlines);
        
        if validStart > nlines
            break;
        end
        
        % 实际读取区域（包含边界重叠）
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
end

%% ========================================================================
%  提取分块对应的SHP索引
%  ========================================================================
function blockSHP = extractBlockSHP(SHP, blockInfo, nlines, nwidths)
    numBlocks = length(blockInfo);
    blockSHP = cell(numBlocks, 1);
    
    for b = 1:numBlocks
        validStartRow = blockInfo(b).globalStartRow;
        validEndRow = blockInfo(b).globalEndRow;
        numValidRows = validEndRow - validStartRow + 1;
        
        % 计算全局像素索引范围
        pixelIndices = zeros(1, numValidRows * nwidths);
        count = 0;
        for col = 1:nwidths
            for row = validStartRow:validEndRow
                count = count + 1;
                pixelIndices(count) = (col-1)*nlines + row;
            end
        end
        blockSHP{b}.PixelInd = SHP.PixelInd(:, pixelIndices);
        blockSHP{b}.CalWin = SHP.CalWin;
    end
end

%% ========================================================================
%  处理单个分块 - Stack 模式
%  ========================================================================
function result = processBlock_Stack(blockMli, blockInf, blockSHP, idx, npages, ...
                                      RadiusRow, RadiusCol, blockInfo, nwidths)
    
    validStartRow = blockInfo.validStartRow;
    validEndRow = blockInfo.validEndRow;
    numValidRows = validEndRow - validStartRow + 1;
    
    % 预分配结果
    result.Coh = zeros(numValidRows, nwidths, npages, 'single');
    result.Ph = zeros(numValidRows, nwidths, npages, 'single');
    
    for ii = 1:npages
        m1 = blockMli(:,:,idx(ii,1));
        m2 = blockMli(:,:,idx(ii,2));
        Intf = sqrt(m1.*m2) .* blockInf(:,:,ii);
        
        m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
        m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
        Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');
        
        nu = zeros(numValidRows, nwidths, 'single');
        de1 = nu; de2 = nu;
        num = 1;
        
        for jj = 1:nwidths
            for kk = 1:numValidRows
                x_global = jj + RadiusCol;
                y_global = kk + validStartRow - 1 + RadiusRow;
                
                MasterValue = m1(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
                SlaveValue = m2(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
                InterfValue = Intf(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
                
                MasterValue = MasterValue(blockSHP.PixelInd(:,num));
                SlaveValue = SlaveValue(blockSHP.PixelInd(:,num));
                InterfValue = InterfValue(blockSHP.PixelInd(:,num));
                
                nu(kk,jj) = mean(InterfValue);
                de1(kk,jj) = mean(MasterValue);
                de2(kk,jj) = mean(SlaveValue);
                num = num + 1;
            end
        end
        
        CohSingle = nu ./ sqrt(de1.*de2);
        CohSingle(isnan(CohSingle)) = 0;
        
        result.Ph(:,:,ii) = angle(CohSingle);
        result.Coh(:,:,ii) = abs(CohSingle);
    end
end

%% ========================================================================
%  处理单个分块 - Average 模式（单页）
%  ========================================================================
function result = processBlock_Average_SinglePage(blockMli, blockInf, blockSHP, ...
                                                   RadiusRow, RadiusCol, blockInfo, nwidths)
    
    validStartRow = blockInfo.validStartRow;
    validEndRow = blockInfo.validEndRow;
    numValidRows = validEndRow - validStartRow + 1;
    
    % 获取当前页的 mli 数据
    m1 = blockMli.m1;
    m2 = blockMli.m2;
    Intf = sqrt(m1.*m2) .* blockInf;
    
    m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
    m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
    Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');
    
    nu = zeros(numValidRows, nwidths, 'single');
    de1 = nu; de2 = nu;
    num = 1;
    
    for jj = 1:nwidths
        for kk = 1:numValidRows
            x_global = jj + RadiusCol;
            y_global = kk + validStartRow - 1 + RadiusRow;
            
            MasterValue = m1(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            SlaveValue = m2(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            InterfValue = Intf(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            
            MasterValue = MasterValue(blockSHP.PixelInd(:,num));
            SlaveValue = SlaveValue(blockSHP.PixelInd(:,num));
            InterfValue = InterfValue(blockSHP.PixelInd(:,num));
            
            nu(kk,jj) = mean(InterfValue);
            de1(kk,jj) = mean(MasterValue);
            de2(kk,jj) = mean(SlaveValue);
            num = num + 1;
        end
    end
    
    CohSingle = nu ./ sqrt(de1.*de2);
    CohSingle(isnan(CohSingle)) = 0;
    result = abs(CohSingle);
end

%% ========================================================================
%  偏差校正 - Stack 模式（串行）
%  ========================================================================
function Coh = applyBiasCorrection_Stack(Coh, SHP, nlines, nwidths, npages, RadiusRow, RadiusCol)
    
    tmp = padarray(Coh, [RadiusRow RadiusCol 0], 'symmetric');
    
    for ii = 1:npages
        num = 1;
        for jj = 1:nwidths
            for kk = 1:nlines
                x_global = jj + RadiusCol;
                y_global = kk + RadiusRow;
                CohValue = tmp(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol, ii);
                CohValue = CohValue(SHP.PixelInd(:,num));
                Coh(kk,jj,ii) = mean(log(CohValue));
                num = num + 1;
            end
        end
        fprintf('bias mitigation: %3d / %d is finished...\n', ii, npages);
    end
    
    Coh = exp(Coh);
end

%% ========================================================================
%  偏差校正 - Stack 模式（并行）
%  ========================================================================
function Coh = applyBiasCorrection_Stack_Parallel(Coh, SHP, nlines, nwidths, npages, ...
                                                   RadiusRow, RadiusCol, NumWorkers)
    
    % 计算分块信息
    blockInfo = computeBlockInfo(nlines, NumWorkers, RadiusRow);
    actualNumBlocks = length(blockInfo);
    
    % 准备分块数据
    blockData_Coh = cell(actualNumBlocks, 1);
    for b = 1:actualNumBlocks
        startRow = blockInfo(b).startRow;
        endRow = blockInfo(b).endRow;
        blockData_Coh{b} = Coh(startRow:endRow, :, :);
    end
    
    % 提取各分块的SHP索引
    blockData_SHP = extractBlockSHP(SHP, blockInfo, nlines, nwidths);
    
    % 并行处理
    blockResults = cell(actualNumBlocks, 1);
    parfor b = 1:actualNumBlocks
        blockResults{b} = processBiasCorrection_Block(blockData_Coh{b}, blockData_SHP{b}, ...
                                                       npages, RadiusRow, RadiusCol, blockInfo(b));
    end
    
    % 合并结果
    for b = 1:actualNumBlocks
        globalStartRow = blockInfo(b).globalStartRow;
        globalEndRow = blockInfo(b).globalEndRow;
        Coh(globalStartRow:globalEndRow, :, :) = blockResults{b};
    end
end

%% ========================================================================
%  处理单个分块的偏差校正
%  ========================================================================
function result = processBiasCorrection_Block(blockCoh, blockSHP, npages, RadiusRow, RadiusCol, blockInfo)
    
    validStartRow = blockInfo.validStartRow;
    validEndRow = blockInfo.validEndRow;
    numValidRows = validEndRow - validStartRow + 1;
    [~, nwidths, ~] = size(blockCoh);
    
    result = zeros(numValidRows, nwidths, npages, 'single');
    
    tmp = padarray(blockCoh, [RadiusRow RadiusCol 0], 'symmetric');
    
    for ii = 1:npages
        num = 1;
        for jj = 1:nwidths
            for kk = 1:numValidRows
                y_global = kk + validStartRow - 1 + RadiusRow;
                x_global = jj + RadiusCol;
                CohValue = tmp(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol, ii);
                CohValue = CohValue(blockSHP.PixelInd(:,num));
                result(kk,jj,ii) = mean(log(CohValue));
                num = num + 1;
            end
        end
    end
    
    result = exp(result);
end

%% ========================================================================
%  偏差校正 - 单幅图像
%  ========================================================================
function CohSingle = applyBiasCorrection_Single(CohSingle, SHP, nlines, nwidths, RadiusRow, RadiusCol)
    
    tmp = padarray(CohSingle, [RadiusRow RadiusCol], 'symmetric');
    CohCorrected = zeros(nlines, nwidths, 'single');
    num = 1;
    
    for jj = 1:nwidths
        for kk = 1:nlines
            x_global = jj + RadiusCol;
            y_global = kk + RadiusRow;
            CohValue = tmp(y_global-RadiusRow:y_global+RadiusRow, x_global-RadiusCol:x_global+RadiusCol);
            CohValue = CohValue(SHP.PixelInd(:,num));
            CohCorrected(kk,jj) = mean(log(CohValue));
            num = num + 1;
        end
    end
    
    CohSingle = exp(CohCorrected);
end