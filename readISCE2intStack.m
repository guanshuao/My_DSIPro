function [intfstack, metaList] = readISCE2intStack(baseDir, varargin)
%READISCE2INTSTACK 按主从日期批量读取 ISCE2 干涉图并堆栈。
%   intfstack = readISCE2intStack(baseDir) 会在 baseDir 下查找形如
%   yymmdd-yymmdd 的文件夹，进入其中的 insar 子目录，寻找 diff*.int
%   干涉图（及同名 .xml/.vrt），调用 readISCE2SLC 逐个读取并堆栈。
%   返回结构体：
%       intfstack.datastack : 行 x 列 x 干涉对数 的复数三维矩阵
%       intfstack.filename  : 每行一个 [master slave] 数字日期
%   [intfstack, metaList] 同时返回每个干涉图的元数据 cell 数组。
%
%   可选 Name-Value 参数（适用时透传给 readISCE2SLC）：
%       'PairFolders' : 指定要读取的 yymmdd-yymmdd 文件夹列表；默认自动检测
%       'InnerPath'   : 干涉图所在子路径，默认 'insar'
%       'Pattern'     : 干涉图文件通配符，默认 'diff*.int'
%       'Subset'      : [startRow, endRow, startCol, endCol] 子区间
%       'Width'       : 覆盖宽度
%       'Length'      : 覆盖长度
%       'ByteOrder'   : 'l' 小端或 'b' 大端（默认 'l'）
%       'DataType'    : 'cfloat'（默认）或 'cfloat64'
%       'Verbose'     : 是否打印进度，默认 true
%
%   示例：
%       baseDir = fullfile(pwd, '..', 'process', 'pais');
%       [intfstack, meta] = readISCE2intStack(baseDir, 'Subset', [1 2000 1 3000]);
%
%   输出字段与 ImgRead 读取的干涉栈一致，可直接用于后续处理。

    p = inputParser;
    addRequired(p, 'baseDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'PairFolders', {}, @(x) iscellstr(x) || isstring(x));
    addParameter(p, 'InnerPath', 'insar', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Pattern', 'diff*.int', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Subset', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 4));
    addParameter(p, 'Width', [], @isnumeric);
    addParameter(p, 'Length', [], @isnumeric);
    addParameter(p, 'ByteOrder', 'l', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DataType', 'cfloat', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Verbose', true, @islogical);
    parse(p, baseDir, varargin{:});

    opts = p.Results;
    baseDir = char(opts.baseDir);
    innerPath = char(opts.InnerPath);
    pattern = char(opts.Pattern);

    if ~isfolder(baseDir)
        error('Base directory does not exist: %s', baseDir);
    end

    % 确定要处理的干涉对文件夹 yymmdd-yymmdd
    if isempty(opts.PairFolders)
        entries = dir(baseDir);
        isPairDir = [entries.isdir];
        names = {entries(isPairDir).name};
        names = names(~ismember(names, {'.', '..'}));
        pairNames = names(~cellfun(@isempty, regexp(names, '^\d{6}-\d{6}$', 'once')));
    else
        pairNames = cellstr(opts.PairFolders);
    end

    if isempty(pairNames)
        error('No pair folders found under %s', baseDir);
    end

    % 在每个干涉对文件夹中定位 diff*.int
    intFiles = {};
    pairDates = [];
    validPairs = {};
    for i = 1:numel(pairNames)
        pairName = pairNames{i};
        parts = regexp(pairName, '^(\d{6})-(\d{6})$', 'tokens', 'once');
        if isempty(parts)
            continue; % 已经过滤，保护性处理
        end
        master = str2double(parts{1});
        slave  = str2double(parts{2});

        searchDir = fullfile(baseDir, pairName, innerPath);
        files = dir(fullfile(searchDir, pattern));
        if isempty(files)
            warning('跳过 %s：未找到匹配的 %s', pairName, pattern);
            continue;
        end

        % 若有多个 diff*.int，默认取排序后的第一个，避免混淆
        [~, sortIdxLocal] = sort({files.name});
        fileChosen = files(sortIdxLocal(1)).name;
        intPath = fullfile(searchDir, fileChosen);

        intFiles{end+1} = intPath; %#ok<AGROW>
        pairDates(end+1, :) = [master, slave]; %#ok<AGROW>
        validPairs{end+1} = pairName; %#ok<AGROW>
    end

    if isempty(intFiles)
        error('未找到任何干涉图文件，请检查路径与命名。');
    end

    % 按主从日期排序
    [~, sortIdx] = sortrows(pairDates, [1 2]);
    intFiles = intFiles(sortIdx);
    pairDates = pairDates(sortIdx, :);
    validPairs = validPairs(sortIdx);

    % 准备透传给 readISCE2SLC 的参数
    forwardArgs = {};
    if ~isempty(opts.Width),     forwardArgs = [forwardArgs, {'Width', opts.Width}]; end
    if ~isempty(opts.Length),    forwardArgs = [forwardArgs, {'Length', opts.Length}]; end
    if ~isempty(opts.Subset),    forwardArgs = [forwardArgs, {'Subset', opts.Subset}]; end
    if ~isempty(opts.ByteOrder), forwardArgs = [forwardArgs, {'ByteOrder', opts.ByteOrder}]; end
    if ~isempty(opts.DataType),  forwardArgs = [forwardArgs, {'DataType', opts.DataType}]; end

    intfstack = struct('datastack', [], 'filename', []);
    metaList = cell(numel(intFiles), 1);

    for k = 1:numel(intFiles)
        intPath = intFiles{k};
        if opts.Verbose
            fprintf('Reading %3d / %3d : %s\n', k, numel(intFiles), intPath);
        end

        [intData, meta] = readISCE2SLC(intPath, forwardArgs{:});

        if k == 1
            [nRows, nCols] = size(intData);
            intfstack.datastack = zeros(nRows, nCols, numel(intFiles), 'like', intData);
            intfstack.filename = zeros(numel(intFiles), 2);
        else
            if ~isequal(size(intData), [size(intfstack.datastack, 1), size(intfstack.datastack, 2)])
                error('尺寸不一致：%s 为 %dx%d，期望 %dx%d', ...
                      intPath, size(intData,1), size(intData,2), ...
                      size(intfstack.datastack,1), size(intfstack.datastack,2));
            end
        end

        intfstack.datastack(:, :, k) = intData;
        intfstack.filename(k, :) = pairDates(k, :);
        metaList{k} = meta;
    end

    if opts.Verbose
        fprintf('已加载 %d 个干涉图，尺寸 %d 行 x %d 列。\n', ...
                numel(intFiles), size(intfstack.datastack,1), size(intfstack.datastack,2));
    end
end
