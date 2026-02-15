function [Coh, Ph] = AdpCohEst_New(slcPath, slcSuffix, slcFormat, slcEndian, ...
                                   intfPath, intfSuffix, intfFormat, intfEndian, ...
                                   nline, SHP, BiasCorr, OutputMode, NumWorkers)
%AdpCohEst_New 内存优化版自适应相干性估计
%   与 AdpCohEst 相比，本函数不再接收庞大的 SLC/干涉堆栈矩阵，
%   而是仅传入文件路径与读入参数，在计算每一对干涉图时按需读取：
%   - 2 张 SLC（转强度）
%   - 1 张干涉图（归一化到单位幅值）
%   从而显著降低峰值内存占用。
%
%   用法:
%   [Coh, Ph] = AdpCohEst_New(slcPath, slcSuffix, slcFormat, slcEndian, ...
%                             intfPath, intfSuffix, intfFormat, intfEndian, ...
%                             nline, SHP, BiasCorr, OutputMode, NumWorkers)
%
%   输入参数（全部必选）:
%   - slcPath:     SLC 文件目录
%   - slcSuffix:   SLC 文件后缀（例如 'slc'）
%   - slcFormat:   SLC 数据格式（例如 'cpxfloat32'）
%   - slcEndian:   SLC 字节序（例如 'b'）
%   - intfPath:    干涉图文件目录
%   - intfSuffix:  干涉图文件后缀（例如 'int'）
%   - intfFormat:  干涉图数据格式（例如 'cpxfloat32'）
%   - intfEndian:  干涉图字节序（例如 'b'）
%   - nline:       图像行数
%   - SHP:         同质像元结构体（至少含 PixelInd, CalWin）
%   - BiasCorr:    偏差矫正开关，'y' 或 'n'
%   - OutputMode:  输出模式，'stack' 或 'average'
%   - NumWorkers:  并行工作进程数，>=1
%
%   输出参数:
%   - Coh:         相干性幅值（stack: 3D；average: 2D）
%   - Ph:          相位（仅 stack 模式有效）
%
%   说明:
%   - 本函数默认使用列优先像元编号，与 SHP.PixelInd 的列组织方式一致。
%   - 若 OutputMode='average' 且调用方请求第二输出，返回空 Ph 并给出提示。

%% 参数检查（全部必选）
if nargin ~= 13
    error(['AdpCohEst_New 需要 13 个必选参数：', ...
           'slcPath, slcSuffix, slcFormat, slcEndian, intfPath, intfSuffix, ', ...
           'intfFormat, intfEndian, nline, SHP, BiasCorr, OutputMode, NumWorkers']);
end

if ~ischar(slcPath) && ~isstring(slcPath)
    error('slcPath 必须是字符串。');
end
if ~ischar(intfPath) && ~isstring(intfPath)
    error('intfPath 必须是字符串。');
end
if ~ischar(slcSuffix) && ~isstring(slcSuffix)
    error('slcSuffix 必须是字符串。');
end
if ~ischar(intfSuffix) && ~isstring(intfSuffix)
    error('intfSuffix 必须是字符串。');
end
if ~ischar(slcFormat) && ~isstring(slcFormat)
    error('slcFormat 必须是字符串。');
end
if ~ischar(intfFormat) && ~isstring(intfFormat)
    error('intfFormat 必须是字符串。');
end
if ~ischar(slcEndian) && ~isstring(slcEndian)
    error('slcEndian 必须是字符串。');
end
if ~ischar(intfEndian) && ~isstring(intfEndian)
    error('intfEndian 必须是字符串。');
end
if ~isscalar(nline) || nline <= 0 || floor(nline) ~= nline
    error('nline 必须是正整数。');
end
if ~isscalar(NumWorkers) || NumWorkers < 1 || floor(NumWorkers) ~= NumWorkers
    error('NumWorkers 必须是 >= 1 的整数。');
end
if ~isstruct(SHP) || ~isfield(SHP, 'PixelInd') || ~isfield(SHP, 'CalWin')
    error('SHP 结构体至少需要包含 PixelInd 和 CalWin 字段。');
end

BiasCorr = lower(char(BiasCorr));
OutputMode = lower(char(OutputMode));
if ~ismember(BiasCorr, {'y', 'n'})
    error('BiasCorr 仅支持 ''y'' 或 ''n''。');
end
if ~ismember(OutputMode, {'stack', 'average'})
    error('OutputMode 仅支持 ''stack'' 或 ''average''。');
end

slcPath = char(slcPath);
intfPath = char(intfPath);
slcSuffix = char(slcSuffix);
intfSuffix = char(intfSuffix);
slcFormat = char(slcFormat);
intfFormat = char(intfFormat);
slcEndian = char(slcEndian);
intfEndian = char(intfEndian);

if ~exist(slcPath, 'dir')
    error('SLC 路径不存在: %s', slcPath);
end
if ~exist(intfPath, 'dir')
    error('干涉图路径不存在: %s', intfPath);
end

%% 初始化
mainTimer = tic;
CalWin = SHP.CalWin;
if numel(CalWin) ~= 2
    error('SHP.CalWin 必须是长度为 2 的窗口大小 [行, 列]。');
end
RadiusRow = (CalWin(1)-1)/2;
RadiusCol = (CalWin(2)-1)/2;
if RadiusRow < 0 || RadiusCol < 0 || floor(RadiusRow) ~= RadiusRow || floor(RadiusCol) ~= RadiusCol
    error('SHP.CalWin 必须是奇数窗口大小，例如 [15 15]。');
end

%% 1) 仅生成文件列表，不读取堆栈
[slcFiles, slcList] = buildFileList(slcPath, slcSuffix, 1);
[intfFiles, intfList] = buildFileList(intfPath, intfSuffix, 2);

if isempty(slcFiles)
    error('在 %s 下未找到后缀为 %s 的 SLC 文件。', slcPath, slcSuffix);
end
if isempty(intfFiles)
    error('在 %s 下未找到后缀为 %s 的干涉图文件。', intfPath, intfSuffix);
end

npages = numel(intfFiles);

% 干涉对到 SLC 日期映射
[tf, idx] = ismember(intfList, slcList);
if any(~tf(:))
    badRows = find(~all(tf, 2));
    error('存在无法在 SLC 列表中匹配的干涉对，首个问题干涉对序号: %d', badRows(1));
end

% 读取第一张 SLC 以获取列数（后续按需读取）
firstSlc = SingleRead(slcFiles{1}, nline, slcFormat, slcEndian);
[~, nwidths] = size(firstSlc);
clear firstSlc;

expectedPixels = nline * nwidths;
if size(SHP.PixelInd, 2) ~= expectedPixels
    error(['SHP.PixelInd 列数与图像尺寸不一致：', ...
           'size(SHP.PixelInd,2)=%d，预期=%d(nline*nwidths)'], ...
           size(SHP.PixelInd,2), expectedPixels);
end

fprintf('SLC 文件数: %d\n', numel(slcFiles));
fprintf('干涉图文件数: %d\n', npages);
fprintf('图像尺寸: %d x %d\n', nline, nwidths);

%% 并行池准备
if NumWorkers > 1
    pool = gcp('nocreate');
    if isempty(pool)
        parpool(NumWorkers);
    elseif pool.NumWorkers ~= NumWorkers
        delete(pool);
        parpool(NumWorkers);
    end
end

%% 2) 预计算分块参数和分块 SHP
blockInfo = computeBlockInfo(nline, NumWorkers, RadiusRow);
actualNumBlocks = length(blockInfo);
blockSHP_all = extractBlockSHP(SHP, blockInfo, nline, nwidths);

%% 输出初始化
if strcmpi(OutputMode, 'stack')
    Coh = zeros(nline, nwidths, npages, 'single');
    if nargout > 1
        Ph = zeros(nline, nwidths, npages, 'single');
    else
        Ph = [];
    end
else
    CohSum = zeros(nline, nwidths, 'single');
    Coh = [];
    Ph = [];
end

%% SLC 两槽缓存（减少重复 I/O）
cache(1).idx = -1;
cache(1).power = [];
cache(2).idx = -1;
cache(2).power = [];

%% 3) 逐干涉对按需读取 + 分块并行
for ii = 1:npages
    m1 = readPowerFromCacheOrDisk(idx(ii,1));
    m2 = readPowerFromCacheOrDisk(idx(ii,2));

    intf = SingleRead(intfFiles{ii}, nline, intfFormat, intfEndian);
    intfAbs = abs(intf);
    intfValid = intfAbs > 0;
    intfNorm = zeros(size(intf), 'like', intf);
    intfNorm(intfValid) = intf(intfValid) ./ intfAbs(intfValid);
    clear intf intfAbs intfValid;

    blockData_mli = cell(actualNumBlocks, 1);
    blockData_inf = cell(actualNumBlocks, 1);

    for b = 1:actualNumBlocks
        startRow = blockInfo(b).startRow;
        endRow = blockInfo(b).endRow;

        blockData_mli{b}.m1 = m1(startRow:endRow, :);
        blockData_mli{b}.m2 = m2(startRow:endRow, :);
        blockData_inf{b} = intfNorm(startRow:endRow, :);
    end

    blockResults = cell(actualNumBlocks, 1);
    blockPhResults = cell(actualNumBlocks, 1);

    if NumWorkers == 1
        for b = 1:actualNumBlocks
            [blockResults{b}, blockPhResults{b}] = processBlock_Average_SinglePage_WithPh( ...
                blockData_mli{b}, blockData_inf{b}, blockSHP_all{b}, RadiusRow, RadiusCol, blockInfo(b), nwidths);
        end
    else
        parfor b = 1:actualNumBlocks
            [blockResults{b}, blockPhResults{b}] = processBlock_Average_SinglePage_WithPh( ...
                blockData_mli{b}, blockData_inf{b}, blockSHP_all{b}, RadiusRow, RadiusCol, blockInfo(b), nwidths);
        end
    end

    CohSingle = zeros(nline, nwidths, 'single');
    PhSingle = zeros(nline, nwidths, 'single');
    for b = 1:actualNumBlocks
        globalStartRow = blockInfo(b).globalStartRow;
        globalEndRow = blockInfo(b).globalEndRow;
        CohSingle(globalStartRow:globalEndRow, :) = blockResults{b};
        PhSingle(globalStartRow:globalEndRow, :) = blockPhResults{b};
    end

    if strcmpi(BiasCorr, 'y')
        CohSingle = applyBiasCorrection_Single(CohSingle, SHP, nline, nwidths, RadiusRow, RadiusCol);
    end

    if strcmpi(OutputMode, 'stack')
        Coh(:,:,ii) = CohSingle;
        if nargout > 1
            Ph(:,:,ii) = PhSingle;
        end
    else
        CohSum = CohSum + CohSingle;
    end

    fprintf('adp. coherence (on-demand): %3d / %d is finished...\n', ii, npages);

    clear m1 m2 intfNorm blockData_mli blockData_inf blockResults blockPhResults CohSingle PhSingle;
end

%% 4) 收尾
if strcmpi(OutputMode, 'average')
    Coh = CohSum / npages;
    if nargout > 1
        warning('average 模式下只输出 Coh，Ph 输出将被忽略。');
        Ph = [];
    end
elseif strcmpi(BiasCorr, 'y')
    % stack 模式的 BiasCorr 已逐页处理，无需额外处理
end

elapsedMin = toc(mainTimer) / 60;
disp(['AdpCohEst_New 运算完成，耗时 ', num2str(elapsedMin, '%.2f'), ' 分钟']);
disp('完成！');


    function powerImg = readPowerFromCacheOrDisk(slcIdx)
        % 先查缓存
        if cache(1).idx == slcIdx
            powerImg = cache(1).power;
            return;
        end
        if cache(2).idx == slcIdx
            powerImg = cache(2).power;
            return;
        end

        % 缓存未命中，读取并计算强度
        slcImg = SingleRead(slcFiles{slcIdx}, nline, slcFormat, slcEndian);
        powerImg = abs(slcImg).^2;
        clear slcImg;

        % 简单双槽替换：把 slot1 挪到 slot2，新数据放 slot1
        cache(2) = cache(1);
        cache(1).idx = slcIdx;
        cache(1).power = powerImg;
    end

end

%% ========================================================================
% 生成文件列表（仅解析文件名，不读影像）
% ========================================================================
function [filePaths, fileDates] = buildFileList(basePath, suffixname, dateCount)
    if ~strcmp(basePath(end), filesep)
        basePath = [basePath, filesep];
    end

    tagFiles = dir([basePath, '*', suffixname]);
    imgNum = length(tagFiles);

    filePaths = cell(imgNum, 1);
    fileDates = zeros(imgNum, dateCount);

    for ii = 1:imgNum
        filePaths{ii} = [basePath, tagFiles(ii).name];
        temp = regexp(tagFiles(ii).name, '\d+', 'match');

        if dateCount == 1
            if ~isscalar(temp)
                error('SLC 文件名格式应包含 1 个日期: %s', tagFiles(ii).name);
            end
            fileDates(ii,1) = str2double(temp{1});
        else
            if length(temp) ~= 2
                error('干涉图文件名格式应包含 2 个日期: %s', tagFiles(ii).name);
            end
            fileDates(ii,1) = str2double(temp{1});
            fileDates(ii,2) = str2double(temp{2});
        end
    end

    % 按日期排序，保证顺序稳定
    if dateCount == 1
        [fileDates, order] = sort(fileDates, 'ascend');
    else
        [fileDates, order] = sortrows(fileDates, [1 2]);
    end
    filePaths = filePaths(order);
end

%% ========================================================================
% 计算分块信息
% ========================================================================
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

        actualStart = max(1, validStart - RadiusRow);
        actualEnd = min(nlines, validEnd + RadiusRow);

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
% 提取分块对应的 SHP 索引
% ========================================================================
function blockSHP = extractBlockSHP(SHP, blockInfo, nlines, nwidths)
    numBlocks = length(blockInfo);
    blockSHP = cell(numBlocks, 1);

    for b = 1:numBlocks
        validStartRow = blockInfo(b).globalStartRow;
        validEndRow = blockInfo(b).globalEndRow;
        numValidRows = validEndRow - validStartRow + 1;

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
% 处理单个分块 - 单页（返回 Coh 与 Ph）
% ========================================================================
function [CohResult, PhResult] = processBlock_Average_SinglePage_WithPh(blockMli, blockInf, blockSHP, ...
                                                                         RadiusRow, RadiusCol, blockInfo, nwidths)

    validStartRow = blockInfo.validStartRow;
    validEndRow = blockInfo.validEndRow;
    numValidRows = validEndRow - validStartRow + 1;

    m1 = blockMli.m1;
    m2 = blockMli.m2;
    Intf = sqrt(m1.*m2) .* blockInf;

    m1 = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
    m2 = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
    Intf = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');

    nu = zeros(numValidRows, nwidths, 'single');
    de1 = nu;
    de2 = nu;
    num = 1;

    for jj = 1:nwidths
        for kk = 1:numValidRows
            xGlobal = jj + RadiusCol;
            yGlobal = kk + validStartRow - 1 + RadiusRow;

            masterValue = m1(yGlobal-RadiusRow:yGlobal+RadiusRow, xGlobal-RadiusCol:xGlobal+RadiusCol);
            slaveValue = m2(yGlobal-RadiusRow:yGlobal+RadiusRow, xGlobal-RadiusCol:xGlobal+RadiusCol);
            interfValue = Intf(yGlobal-RadiusRow:yGlobal+RadiusRow, xGlobal-RadiusCol:xGlobal+RadiusCol);

            masterValue = masterValue(blockSHP.PixelInd(:,num));
            slaveValue = slaveValue(blockSHP.PixelInd(:,num));
            interfValue = interfValue(blockSHP.PixelInd(:,num));

            nu(kk,jj) = mean(interfValue);
            de1(kk,jj) = mean(masterValue);
            de2(kk,jj) = mean(slaveValue);
            num = num + 1;
        end
    end

    CohComplex = nu ./ sqrt(de1.*de2);
    CohComplex(isnan(CohComplex)) = 0;

    PhResult = angle(CohComplex);
    CohResult = abs(CohComplex);
end

%% ========================================================================
% 偏差校正 - 单幅图像
% ========================================================================
function CohSingle = applyBiasCorrection_Single(CohSingle, SHP, nlines, nwidths, RadiusRow, RadiusCol)

    tmp = padarray(CohSingle, [RadiusRow RadiusCol], 'symmetric');
    CohCorrected = zeros(nlines, nwidths, 'single');
    num = 1;

    for jj = 1:nwidths
        for kk = 1:nlines
            xGlobal = jj + RadiusCol;
            yGlobal = kk + RadiusRow;
            CohValue = tmp(yGlobal-RadiusRow:yGlobal+RadiusRow, xGlobal-RadiusCol:xGlobal+RadiusCol);
            CohValue = CohValue(SHP.PixelInd(:,num));
            CohCorrected(kk,jj) = mean(log(CohValue));
            num = num + 1;
        end
    end

    CohSingle = exp(CohCorrected);
end
