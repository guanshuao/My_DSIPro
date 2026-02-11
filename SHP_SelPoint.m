function [SHP]=SHP_SelPoint(mlistack,CalWin,Alpha,EstAgr,NumBlocks)
%SHP_SelPoint 并行版本的同质像元(SHP)选取函数
%   本函数用于在SAR强度图堆栈上选取统计同质像元(Statistically Homogeneous Pixels)
%   针对BWS算法进行了并行优化，大幅提升计算速度
%
%   用法:
%       [SHP]=SHP_SelPoint(mlistack,CalWin,Alpha,EstAgr,NumBlocks)
%   
%   输入参数:
%   - mlistack: 高度 × 宽度 × 页数 的三维强度图堆栈
%   - CalWin:   固定的盒式窗口大小 [行数, 列数]，默认 [15, 15]
%   - Alpha:    显著性水平，取值范围0到1，默认0.05（即5%显著性水平）
%   - EstAgr:   像元选取算法: 'FaSHPS' 或 'BWS'（仅BWS支持Alpha=0.05和0.01）
%               默认 'FaSHPS'
%   - NumBlocks: 并行分块数量，默认为可用CPU核心数
%               仅对BWS算法有效，FaSHPS算法忽略此参数
%
%   输出参数:
%   - SHP.PixelInd: CalWin(1)*CalWin(2) × (高度*宽度) 的逻辑矩阵
%                   每列包含对应像素的同质像元集合
%   - SHP.BroNum:   每个像素的同质像元数量（不包含参考像素本身）
%   - SHP.CalWin:   固定盒式窗口大小
%
%   参考文献:
%   [1] Distributed scatterer interferometry with the refinement of spatiotemporal coherence
%        Mi Jiang, Andrea Monti-Guarnieri
%        IEEE Transactions on Geoscience and Remote Sensing vol. 58, no. 6, 
%        June 2020, pp. 3977-3987
%
%   [2] Fast Statistically Homogeneous Pixel Selection for Covariance Matrix Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding, Ramon F. Hanssen, Rakesh Malhotra and Ling Chang,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 53, no. 3,
%        March 2015, pp. 1213-1224.
% 
%   [3] The potential of more accurate InSAR covariance matrix estimation for land cover mapping
%        Mi Jiang, Bin Yong, Xin Tian, Rakesh Malhotra, Rui Hu, Zhiwei Li, Zhongbo Yu and Xinxin Zhang,
%        ISPRS Journal of Photogrammetry and Remote Sensing Vol. 126,
%        April 2017, pp. 120-128.
% 
%   提示: 当使用N视强度图像时，只能使用BWS算法！
%
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University
%   并行版本优化
%   ======================================================================
%   01/2026 并行版本：针对BWS算法进行并行优化
%   ======================================================================

%% 参数默认值设置
if nargin < 5 || isempty(NumBlocks)
    % 默认为1，即串行模式（与原始SHP_SelPoint兼容）
    NumBlocks = 1;
end

if nargin < 4 || isempty(EstAgr)
    EstAgr = 'BWS';
end

if nargin < 3 || isempty(Alpha)
    Alpha = 0.05;
end

if nargin < 2
    help SHP_SelPoint
    return;
end

%% 输入验证
if length(size(mlistack)) ~= 3
    error('请输入三维矩阵...');
end

tic;

[nlines, nwidths, npages] = size(mlistack);
mlistack = single(mlistack);

% 窗口参数计算
RadiusRow = (CalWin(1)-1)/2;
RadiusCol = (CalWin(2)-1)/2;
InitRow = (CalWin(1)+1)/2;  % 中心行位置
InitCol = (CalWin(2)+1)/2;  % 中心列位置

%% 根据算法类型选择处理方式
if strcmpi(EstAgr, 'FaSHPS')
    % FaSHPS算法：使用原始串行实现（已经很快）
    disp('使用FaSHPS算法（串行模式）...');
    SHP = SHP_SelPoint_FaSHPS(mlistack, CalWin, Alpha, nlines, nwidths, npages, ...
                               RadiusRow, RadiusCol, InitRow, InitCol);
else
    % BWS算法
    if NumBlocks == 1
        % 串行模式：与原始SHP_SelPoint完全一致
        disp('使用BWS算法（串行模式）...');
        SHP = SHP_SelPoint_BWS_Serial(mlistack, CalWin, Alpha, nlines, nwidths, npages, ...
                                       RadiusRow, RadiusCol, InitRow, InitCol);
    else
        % 并行模式
        disp(['使用BWS算法（并行模式，', num2str(NumBlocks), '个分块）...']);
        SHP = SHP_SelPoint_BWS_Parallel(mlistack, CalWin, Alpha, nlines, nwidths, npages, ...
                                         RadiusRow, RadiusCol, InitRow, InitCol, NumBlocks);
    end
end

%% 输出结果
t = toc;
figure;
imagesc(SHP.BroNum);
axis image off;
colormap jet;
ti = title('同质像元数量');
colorbar;
set(ti, 'fontweight', 'bold');
disp(['SHP_SelPoint 运算完成，耗时 ', num2str(t/60, '%.2f'), ' 分钟']);
disp('完成！');

end

%% ========================================================================
%  FaSHPS算法实现（串行版本，与原函数一致）
%  ========================================================================
function SHP = SHP_SelPoint_FaSHPS(mlistack, CalWin, Alpha, nlines, nwidths, npages, ...
                                    RadiusRow, RadiusCol, InitRow, InitCol)
    
    % 边缘镜像扩展
    mlistack = padarray(mlistack, [RadiusRow RadiusCol], 'symmetric');
    meanmli = mean(mlistack, 3);
    [nlines_EP, nwidths_EP] = size(meanmli);
    SHP.PixelInd = false(CalWin(1)*CalWin(2), nlines*nwidths);
    
    % 似然比检验的窗口半径
    LRT_nl = 3;
    LRT_nw = 3;
    if RadiusRow < LRT_nl
        LRT_nl = 1;
    end
    if RadiusCol < LRT_nw
        LRT_nw = 1;
    end
    
    % 临界区域计算
    CR_lo = finv(Alpha/2, 2*npages, 2*npages);
    CR_up = finv(1-Alpha/2, 2*npages, 2*npages);
    Galpha_L = gaminv(Alpha/2, npages, 1);
    Galpha_U = gaminv(1-Alpha/2, npages, 1);
    
    % 进度显示
    num = 1;
    p = 1;
    all_pixels = nlines * nwidths;
    all_step = floor(all_pixels/10);
    
    for kk = InitCol:nwidths_EP-RadiusCol
        for ll = InitRow:nlines_EP-RadiusRow
            % 初始估计（似然比检验）
            temp = meanmli(ll-LRT_nl:ll+LRT_nl, kk-LRT_nw:kk+LRT_nw);
            T = meanmli(ll, kk) ./ temp;
            T = T > CR_lo & T < CR_up;
            SeedPoint = mean(temp(T));
            
            % 迭代（Gamma置信区间）
            MeanMatrix = meanmli(ll-RadiusRow:ll+RadiusRow, kk-RadiusCol:kk+RadiusCol);
            SeedPoint = MeanMatrix > Galpha_L*SeedPoint/npages & MeanMatrix < Galpha_U*SeedPoint/npages;
            SeedPoint(InitRow, InitCol) = true;
            
            % 连通性分析
            LL = bwlabel(SeedPoint);
            SHP.PixelInd(:, num) = LL(:) == LL(InitRow, InitCol);
            num = num + 1;
            
            if num == all_step * p
                disp(['进度: ', num2str(10*p), '%']);
                p = p + 1;
            end
        end
    end
    
    % 计算同质像元数量
    SHP.BroNum = sum(SHP.PixelInd, 1);
    SHP.BroNum = uint16(reshape(SHP.BroNum(:), [nlines, nwidths]));
    SHP.CalWin = CalWin;
end

%% ========================================================================
%  BWS算法串行实现（与原始SHP_SelPoint完全一致）
%  ========================================================================
function SHP = SHP_SelPoint_BWS_Serial(mlistack, CalWin, Alpha, nlines, nwidths, npages, ...
                                        RadiusRow, RadiusCol, InitRow, InitCol)
    
    % 边缘镜像扩展
    mlistack = padarray(mlistack, [RadiusRow RadiusCol], 'symmetric');
    [nlines_EP, nwidths_EP, ~] = size(mlistack);
    SHP.PixelInd = false(CalWin(1)*CalWin(2), nlines*nwidths);
    
    % 进度显示
    num = 1;
    p = 1;
    all_pixels = nlines * nwidths;
    all_step = floor(all_pixels/10);
    
    for kk = InitCol:nwidths_EP-RadiusCol
        % 展示进度和当前时间
        disp(['Processing column ', num2str(kk-InitCol+1), ' of ', num2str(nwidths), ' at ', datestr(now,'HH:MM:SS')]);
        for ll = InitRow:nlines_EP-RadiusRow
            Matrix = mlistack(ll-RadiusRow:ll+RadiusRow, kk-RadiusCol:kk+RadiusCol, :);
            Ref = Matrix(InitRow, InitCol, :);
            T = BWStest(repmat(Ref(:), [1, CalWin(1)*CalWin(2)]), ...
                        reshape(Matrix, [CalWin(1)*CalWin(2), npages])', Alpha);
            temp = reshape(~T, [CalWin(1), CalWin(2)]);
            % 连通性分析
            LL = bwlabel(temp);
            SHP.PixelInd(:, num) = LL(:) == LL(InitRow, InitCol);
            num = num + 1;
            if num == all_step * p
                disp(['progress: ', num2str(10*p), '%']);
                p = p + 1;
            end
        end
    end
    
    % 计算同质像元数量
    SHP.BroNum = sum(SHP.PixelInd, 1);
    SHP.BroNum = uint16(reshape(SHP.BroNum(:), [nlines, nwidths]));
    SHP.CalWin = CalWin;
end

%% ========================================================================
%  BWS算法并行实现
%  ========================================================================
function SHP = SHP_SelPoint_BWS_Parallel(mlistack, CalWin, Alpha, nlines, nwidths, npages, ...
                                          RadiusRow, RadiusCol, InitRow, InitCol, NumBlocks)
    
    % 确保并行池已启动
    pool = gcp('nocreate');
    if isempty(pool)
        disp('正在启动并行池...');
        parpool('local', NumBlocks);
    else
        disp(['并行池已就绪，工作进程数: ', num2str(pool.NumWorkers)]);
    end
    
    %% 计算分块参数
    % 重要：分块时需要考虑边界重叠，重叠量等于窗口半径
    % 这样可以保证每个分块边缘的像素也能正确计算
    
    % 按行分块（更适合MATLAB的列优先存储）
    % 计算每个分块的有效行数（不含重叠区域）
    validRowsPerBlock = ceil(nlines / NumBlocks);
    
    % 预分配分块信息
    blockInfo = struct('startRow', {}, 'endRow', {}, ...
                       'validStartRow', {}, 'validEndRow', {}, ...
                       'globalStartRow', {});
    
    for b = 1:NumBlocks
        % 当前分块的有效区域（在原始图像中的位置）
        validStart = (b-1) * validRowsPerBlock + 1;
        validEnd = min(b * validRowsPerBlock, nlines);
        
        if validStart > nlines
            break;  % 如果分块数过多，可能有空分块
        end
        
        % 实际需要读取的区域（包含重叠边界）
        % 需要向上和向下各扩展RadiusRow行
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
    
    %% 并行处理各分块
    % 预分配结果存储
    blockResults = cell(actualNumBlocks, 1);
    
    % 准备各分块的数据
    blockData = cell(actualNumBlocks, 1);
    for b = 1:actualNumBlocks
        % 提取当前分块的数据（包含边界重叠）
        blockData{b} = mlistack(blockInfo(b).startRow:blockInfo(b).endRow, :, :);
    end
    
    % 清理原始大矩阵以释放内存（数据已转移至blockData）
    clear mlistack;

    % 并行处理
    disp('开始并行计算...');
    parfor b = 1:actualNumBlocks
        fprintf('分块 %d/%d 开始处理...\n', b, actualNumBlocks);
        
        % 获取当前分块数据
        currentBlock = blockData{b};
        currentBlockInfo = blockInfo(b);
        
        % 处理当前分块
        blockResults{b} = processSingleBlock_BWS(currentBlock, CalWin, Alpha, npages, ...
                                                   RadiusRow, RadiusCol, InitRow, InitCol, ...
                                                   currentBlockInfo, nwidths, b);
        
        fprintf('分块 %d/%d 处理完成\n', b, actualNumBlocks);
    end
    
    %% 合并结果
    % 清理未使用的变量以释放内存
    clear blockData;
    
    disp('正在合并分块结果...');
    SHP.PixelInd = false(CalWin(1)*CalWin(2), nlines*nwidths);
    
    for b = 1:actualNumBlocks
        % 获取当前分块的有效结果
        validResult = blockResults{b};
        blockResults{b} = []; % 立即清空，释放内存
        
        % 计算全局索引范围
        globalStartRow = blockInfo(b).globalStartRow;
        globalEndRow = blockInfo(b).globalEndRow;
        numValidRows = globalEndRow - globalStartRow + 1;
        
        % 将结果放入正确的位置 - 向量化优化以提升速度和内存效率
        for col = 1:nwidths
            % 计算源和目标索引范围
            % 源数据（validResult）是列优先排列的块数据
            srcStart = (col-1)*numValidRows + 1;
            srcEnd = srcStart + numValidRows - 1;
            
            % 目标数据（SHP.PixelInd）也是列优先排列的全图数据
            dstStart = (col-1)*nlines + globalStartRow;
            dstEnd = dstStart + numValidRows - 1;
            
            % 向量化赋值
            SHP.PixelInd(:, dstStart:dstEnd) = validResult(:, srcStart:srcEnd);
        end
        
        % 清理临时变量
        clear validResult;
    end
    clear blockResults; % 彻底清除容器
    
    %% 计算同质像元数量
    SHP.BroNum = sum(SHP.PixelInd, 1);
    SHP.BroNum = uint16(reshape(SHP.BroNum(:), [nlines, nwidths]));
    SHP.CalWin = CalWin;
end

%% ========================================================================
%  处理单个分块的BWS计算
%  ========================================================================
function result = processSingleBlock_BWS(blockData, CalWin, Alpha, npages, ...
                                          RadiusRow, RadiusCol, InitRow, InitCol, ...
                                          blockInfo, nwidths, blockIdx)
    
    % 获取分块尺寸
    [blockLines, ~, ~] = size(blockData);
    
    % 边缘镜像扩展（仅在列方向，行方向已经有足够的重叠）
    % 但为了保持算法一致性，仍然进行完整的扩展
    blockData = padarray(blockData, [RadiusRow RadiusCol], 'symmetric');
    [nlines_EP, nwidths_EP, ~] = size(blockData);
    
    % 计算有效处理区域（在扩展后坐标系中）
    % 原始分块的第一行在扩展后位于 RadiusRow+1
    % 有效区域的起始行
    validStartInPadded = RadiusRow + blockInfo.validStartRow;
    validEndInPadded = RadiusRow + blockInfo.validEndRow;
    
    % 计算有效行数和总像素数
    numValidRows = blockInfo.validEndRow - blockInfo.validStartRow + 1;
    numValidPixels = numValidRows * nwidths;
    
    % 预分配结果
    result = false(CalWin(1)*CalWin(2), numValidPixels);
    
    % 处理每个像素
    pixelCount = 0;
    for kk = InitCol:nwidths_EP-RadiusCol  % 列循环
        for ll = validStartInPadded:validEndInPadded  % 仅处理有效行
            pixelCount = pixelCount + 1;
            
            % 提取窗口数据
            Matrix = blockData(ll-RadiusRow:ll+RadiusRow, kk-RadiusCol:kk+RadiusCol, :);
            Ref = Matrix(InitRow, InitCol, :);
            
            % BWS检验
            T = BWStest(repmat(Ref(:), [1, CalWin(1)*CalWin(2)]), ...
                        reshape(Matrix, [CalWin(1)*CalWin(2), npages])', Alpha);
            temp = reshape(~T, [CalWin(1), CalWin(2)]);
            
            % 连通性分析
            LL = bwlabel(temp);
            result(:, pixelCount) = LL(:) == LL(InitRow, InitCol);
        end
    end
    
    % 显示进度
    fprintf('  分块 %d: 处理了 %d 个像素\n', blockIdx, pixelCount);
end