function [slcstack, metaList] = readISCE2SLCStack(baseDir, varargin)
%READISCE2SLCSTACK 按日期批量读取 ISCE2 生成的 SLC 并堆栈。
%   slcstack = readISCE2SLCStack(baseDir) 会在 baseDir 下自动查找形如
%   yymmdd 的日期文件夹，读取 <date>/f1_2960/s1/<date>.slc（调用
%   readISCE2SLC），返回包含以下字段的结构体：
%       slcstack.datastack : 行 x 列 x 日期数 的复数三维矩阵
%       slcstack.filename  : 来自文件夹名的数字日期
%   [slcstack, metaList] 同时返回每个日期的元数据 cell 数组。
%
%   可选 Name-Value 参数（适用时透传给 readISCE2SLC）：
%       'DateFolders' : 指定要读取的日期文件夹列表；默认自动检测 yymmdd
%       'InnerPath'   : 日期文件夹内通往 SLC 的子路径，默认 f1_2960/s1
%       'Subset'      : [startRow, endRow, startCol, endCol] 子区间
%       'Width'       : 覆盖宽度
%       'Length'      : 覆盖长度
%       'ByteOrder'   : 'l' 小端或 'b' 大端（默认 'l'）
%       'DataType'    : 'cfloat'（默认）或 'cfloat64'
%       'Verbose'     : 是否打印进度，默认 true
%
%   示例：
%       baseDir = fullfile(pwd, '..', 'process', 'dates_resampled');
%       [slcstack, meta] = readISCE2SLCStack(baseDir, 'Subset', [1 2000 1 3000]);
%
%   输出与 ImgRead 保持一致：datastack 和 filename 可直接用于
%   Kumamoto_SHP 等后续脚本。

    p = inputParser;
    addRequired(p, 'baseDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DateFolders', {}, @(x) iscellstr(x) || isstring(x));
    addParameter(p, 'InnerPath', fullfile('f1_2960', 's1'), @(x) ischar(x) || isstring(x));
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

    if ~isfolder(baseDir)
        error('Base directory does not exist: %s', baseDir);
    end

    % 确定要处理的日期文件夹
    if isempty(opts.DateFolders)
        entries = dir(baseDir);
        isDateDir = [entries.isdir];
        names = {entries(isDateDir).name};
        names = names(~ismember(names, {'.', '..'}));
        dateNames = names(~cellfun(@isempty, regexp(names, '^\d{6}$', 'once')));
    else
        dateNames = cellstr(opts.DateFolders);
    end

    if isempty(dateNames)
        error('No date folders found under %s', baseDir);
    end

    % 拼接各日期的 SLC 完整路径
    slcFiles = {};
    validDates = {};
    for i = 1:numel(dateNames)
        d = dateNames{i};
        slcPath = fullfile(baseDir, d, innerPath, [d, '.slc']);
        if exist(slcPath, 'file')
            slcFiles{end+1} = slcPath; %#ok<AGROW>
            validDates{end+1} = d; %#ok<AGROW>
        else
            warning('Skip %s: SLC not found at %s', d, slcPath);
        end
    end

    if isempty(slcFiles)
        error('No SLC files were found. Check path and folder names.');
    end

    % 按日期排序（yymmdd）
    [~, sortIdx] = sort(cellfun(@str2double, validDates));
    slcFiles = slcFiles(sortIdx);
    validDates = validDates(sortIdx);

    % 准备透传给 readISCE2SLC 的参数
    forwardArgs = {};
    if ~isempty(opts.Width),     forwardArgs = [forwardArgs, {'Width', opts.Width}]; end
    if ~isempty(opts.Length),    forwardArgs = [forwardArgs, {'Length', opts.Length}]; end
    if ~isempty(opts.Subset),    forwardArgs = [forwardArgs, {'Subset', opts.Subset}]; end
    if ~isempty(opts.ByteOrder), forwardArgs = [forwardArgs, {'ByteOrder', opts.ByteOrder}]; end
    if ~isempty(opts.DataType),  forwardArgs = [forwardArgs, {'DataType', opts.DataType}]; end

    slcstack = struct('datastack', [], 'filename', []);
    metaList = cell(numel(slcFiles), 1);

    for k = 1:numel(slcFiles)
        slcPath = slcFiles{k};
        if opts.Verbose
            fprintf('Reading %3d / %3d : %s\n', k, numel(slcFiles), slcPath);
        end

        [slcData, meta] = readISCE2SLC(slcPath, forwardArgs{:});

        if k == 1
            [nRows, nCols] = size(slcData);
            slcstack.datastack = zeros(nRows, nCols, numel(slcFiles), 'like', slcData);
            slcstack.filename = zeros(numel(slcFiles), 1);
        else
            if ~isequal(size(slcData), [size(slcstack.datastack, 1), size(slcstack.datastack, 2)])
                error('Size mismatch at %s (got %dx%d, expected %dx%d)', ...
                      slcPath, size(slcData,1), size(slcData,2), ...
                      size(slcstack.datastack,1), size(slcstack.datastack,2));
            end
        end

        slcstack.datastack(:, :, k) = slcData;
        slcstack.filename(k) = str2double(validDates{k});
        metaList{k} = meta;
    end

    if opts.Verbose
        fprintf('已加载 %d 个 SLC，尺寸 %d 行 x %d 列。\n', ...
                numel(slcFiles), size(slcstack.datastack,1), size(slcstack.datastack,2));
    end
end
