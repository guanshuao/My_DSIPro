function BoxCohEst_StreamByPair(datesDir, pairsDir, varargin)
% BoxCohEst_StreamByPair 逐对读取主从SLC与干涉图，估计并保存Boxcar相干性。

    p = inputParser;
    addRequired(p, 'datesDir', @(x) ischar(x) || isstring(x));
    addRequired(p, 'pairsDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Subset', [], @(x) isnumeric(x) && (isempty(x) || numel(x)==4));
    addParameter(p, 'SLCInnerPath', fullfile('f1_2960','s1'), @(x) ischar(x) || isstring(x));
    addParameter(p, 'PairsInnerPath', 'insar', @(x) ischar(x) || isstring(x));
    addParameter(p, 'IntPattern', 'diff*.int', @(x) ischar(x) || isstring(x));
    addParameter(p, 'CalWin', [7 7], @(x) (isnumeric(x) && (numel(x)==2 || size(x,2)==2)) || iscell(x));
    addParameter(p, 'CohEstAgr', 'fast', @(x) ischar(x) || isstring(x));
    addParameter(p, 'OutDir', fullfile(pwd,'BoxCoh'), @(x) ischar(x) || isstring(x) || iscell(x));
    addParameter(p, 'Verbose', true, @islogical);
    parse(p, datesDir, pairsDir, varargin{:});

    opts = p.Results;
    datesDir = char(opts.datesDir);
    pairsDir = char(opts.pairsDir);
    slcInner = char(opts.SLCInnerPath);
    pairInner = char(opts.PairsInnerPath);
    intPattern = char(opts.IntPattern);
    subset = opts.Subset;
    method = lower(char(opts.CohEstAgr));

    % 处理多窗口与输出目录
    calWinList = normalizeCalWin(opts.CalWin);
    outDirList = normalizeOutDir(opts.OutDir, calWinList);

    % 构建日期->SLC路径映射
    dateEntries = dir(datesDir);
    dateNames = {dateEntries([dateEntries.isdir]).name};
    dateNames = dateNames(~ismember(dateNames,{'.','..'}));
    dateNames = dateNames(~cellfun(@isempty, regexp(dateNames,'^\d{6}$','once')));
    slcMapKeys = zeros(0,1);
    slcMapVals = {};
    for i = 1:numel(dateNames)
        d = dateNames{i};
        slcPath = fullfile(datesDir, d, slcInner, [d '.slc']);
        if exist(slcPath,'file')
            slcMapKeys(end+1,1) = str2double(d); %#ok<AGROW>
            slcMapVals{end+1,1} = slcPath; %#ok<AGROW>
        end
    end
    if isempty(slcMapKeys)
        error('No SLC files found under %s', datesDir);
    end

    % 构建干涉对列表及路径
    pairEntries = dir(pairsDir);
    pairNames = {pairEntries([pairEntries.isdir]).name};
    pairNames = pairNames(~ismember(pairNames,{'.','..'}));
    pairNames = pairNames(~cellfun(@isempty, regexp(pairNames,'^\d{6}-\d{6}$','once')));
    pairDates = zeros(0,2);
    intPaths = {};
    for i = 1:numel(pairNames)
        pn = pairNames{i};
        toks = regexp(pn,'^(\d{6})-(\d{6})$','tokens','once');
        if isempty(toks), continue; end
        m = str2double(toks{1});
        s = str2double(toks{2});
        searchDir = fullfile(pairsDir, pn, pairInner);
        files = dir(fullfile(searchDir, intPattern));
        if isempty(files),
            if opts.Verbose, warning('Skip %s: no %s', pn, intPattern); end
            continue;
        end
        [~, idx] = sort({files.name});
        intPath = fullfile(searchDir, files(idx(1)).name);
        pairDates(end+1,:) = [m s]; %#ok<AGROW>
        intPaths{end+1,1} = intPath; %#ok<AGROW>
    end
    if isempty(intPaths)
        error('No interferograms found under %s', pairsDir);
    end

    % 按主从日期排序
    [~, sortIdx] = sortrows(pairDates,[1 2]);
    pairDates = pairDates(sortIdx,:);
    intPaths = intPaths(sortIdx);

    hList = cellfun(@(w) fspecial('average', w), calWinList, 'UniformOutput', false);

    % 缓存当前参考（主）图像，减少重复读取
    currentMasterDate = NaN;
    currentMasterSLC = [];
    for k = 1:size(pairDates,1)
        m = pairDates(k,1);
        s = pairDates(k,2);
        % 定位主从SLC
        mIdx = find(slcMapKeys==m,1);
        sIdx = find(slcMapKeys==s,1);
        if isempty(mIdx) || isempty(sIdx)
            warning('Missing SLC for pair %06d-%06d, skip.', m, s);
            continue;
        end
        mPath = slcMapVals{mIdx};
        sPath = slcMapVals{sIdx};

        % 读取参考主图像（若与上一次不同则更新缓存）
        if ~isequal(currentMasterDate, m) || isempty(currentMasterSLC)
            [currentMasterSLC, ~] = readISCE2SLC(mPath, 'Subset', subset);
            currentMasterDate = m;
        end
        % 仅读取本次的从图与干涉图
        [sslc, ~] = readISCE2SLC(sPath, 'Subset', subset);
        [intData, ~] = readISCE2SLC(intPaths{k}, 'Subset', subset);

        % 归一化干涉图至单位振幅
        amp = abs(intData);
        mask = amp>0;
        intData(~mask) = 0;
        intData(mask) = intData(mask)./amp(mask);

        m1 = abs(currentMasterSLC).^2;
        m2 = abs(sslc).^2;

        for wIdx = 1:numel(calWinList)
            win = calWinList{wIdx};
            coh = computeCoh(m1, m2, intData, win, method, hList{wIdx});

            % coh = single(abs(coh));
            pairName = sprintf('%06d_%06d.mat', m, s);
            save(fullfile(outDirList{wIdx}, pairName), 'coh', '-v7.3');

            if opts.Verbose
                fprintf('Saved BoxCoh for %06d-%06d (win %dx%d) to %s\n', m, s, win(1), win(2), outDirList{wIdx});
            end
        end
    end

    if opts.Verbose
        disp('BoxCohEst_StreamByPair completed.');
    end
end

function calWinList = normalizeCalWin(calWin)
    if iscell(calWin)
        calWinList = cellfun(@(w) reshape(w, 1, []), calWin, 'UniformOutput', false);
    elseif isnumeric(calWin) && numel(calWin)==2
        calWinList = {reshape(calWin, 1, [])};
    elseif isnumeric(calWin) && size(calWin,2)==2
        calWinList = mat2cell(calWin, ones(size(calWin,1),1), 2);
    else
        error('CalWin must be [r c], Nx2, or a cell array of [r c].');
    end
end

function outDirList = normalizeOutDir(outDir, calWinList)
    if iscell(outDir)
        outDirList = outDir(:);
    elseif isstring(outDir) && numel(outDir) > 1
        outDirList = cellstr(outDir(:));
    else
        outDirList = {char(outDir)};
    end

    if numel(outDirList) == 1 && numel(calWinList) > 1
        baseDir = outDirList{1};
        outDirList = cell(numel(calWinList),1);
        for i = 1:numel(calWinList)
            w = calWinList{i};
            outDirList{i} = sprintf('%s_%dx%d', baseDir, w(1), w(2));
        end
    elseif numel(outDirList) ~= numel(calWinList)
        error('OutDir count must match CalWin count when providing multiple windows.');
    end

    for i = 1:numel(outDirList)
        if ~exist(outDirList{i},'dir'), mkdir(outDirList{i}); end
    end
end

function coh = computeCoh(m1, m2, intData, win, method, h)
    switch method
        case 'fast'
            nu = imfilter(sqrt(m1.*m2).*intData, h, 'symmetric');
            de1 = imfilter(m1, h, 'symmetric');
            de2 = imfilter(m2, h, 'symmetric');
            coh = nu ./ sqrt(de1 .* de2);
        case 'normal'
            [nlines, nwidths] = size(m1);
            RadiusRow = (win(1) - 1) / 2;
            RadiusCol = (win(2) - 1) / 2;
            m1p = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
            m2p = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
            intp = padarray(intData, [RadiusRow RadiusCol], 'symmetric');
            coh = zeros(nlines, nwidths, 'single');
            for jj = 1:nwidths
                for kk = 1:nlines
                    x = jj + RadiusCol;
                    y = kk + RadiusRow;
                    m1w = m1p(y-RadiusRow:y+RadiusRow, x-RadiusCol:x+RadiusCol);
                    m2w = m2p(y-RadiusRow:y+RadiusRow, x-RadiusCol:x+RadiusCol);
                    nuw = sqrt(m1w.*m2w) .* intp(y-RadiusRow:y+RadiusRow, x-RadiusCol:x+RadiusCol);
                    coh(kk,jj) = mean(nuw(:)) ./ sqrt(mean(m1w(:)) .* mean(m2w(:)));
                end
            end
        otherwise
            error('Unsupported CohEstAgr: %s', method);
    end
end
