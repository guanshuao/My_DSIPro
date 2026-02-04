function AdpCohEst_StreamByPair(datesDir, pairsDir, SHP, BiasCorr, outDir, varargin)
% AdpCohEst_StreamByPair 逐对读取主从SLC与干涉图，采用同质像素估计相干性并保存。

    p = inputParser;
    addRequired(p, 'datesDir', @(x) ischar(x) || isstring(x));
    addRequired(p, 'pairsDir', @(x) ischar(x) || isstring(x));
    addRequired(p, 'SHP', @isstruct);
    addRequired(p, 'BiasCorr', @(x) ischar(x) || isstring(x));
    addRequired(p, 'outDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Subset', [], @(x) isnumeric(x) && (isempty(x) || numel(x)==4));
    addParameter(p, 'SLCInnerPath', fullfile('f1_2960','s1'), @(x) ischar(x) || isstring(x));
    addParameter(p, 'PairsInnerPath', 'insar', @(x) ischar(x) || isstring(x));
    addParameter(p, 'IntPattern', 'diff*.int', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Verbose', true, @islogical);
    parse(p, datesDir, pairsDir, SHP, BiasCorr, outDir, varargin{:});

    opts = p.Results;
    datesDir = char(opts.datesDir);
    pairsDir = char(opts.pairsDir);
    slcInner = char(opts.SLCInnerPath);
    pairInner = char(opts.PairsInnerPath);
    intPattern = char(opts.IntPattern);
    subset = opts.Subset;
    BiasCorr = char(opts.BiasCorr);
    outDir = char(opts.outDir);

    if ~exist(outDir,'dir'), mkdir(outDir); end

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
        if isempty(files)
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

    CalWin = SHP.CalWin;
    RadiusRow = (CalWin(1) - 1) / 2;
    RadiusCol = (CalWin(2) - 1) / 2;

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
        Intf = sqrt(m1 .* m2) .* intData;

        [nlines, nwidths] = size(m1);
        m1p = padarray(m1, [RadiusRow RadiusCol], 'symmetric');
        m2p = padarray(m2, [RadiusRow RadiusCol], 'symmetric');
        Intfp = padarray(Intf, [RadiusRow RadiusCol], 'symmetric');

        nu = zeros(nlines, nwidths, 'single');
        de1 = nu;
        de2 = nu;
        num = 1;

        for jj = 1:nwidths
            for kk = 1:nlines
                x_global = jj + RadiusCol;
                y_global = kk + RadiusRow;

                MasterValue = m1p(y_global - RadiusRow:y_global + RadiusRow, x_global - RadiusCol:x_global + RadiusCol);
                SlaveValue = m2p(y_global - RadiusRow:y_global + RadiusRow, x_global - RadiusCol:x_global + RadiusCol);
                InterfValue = Intfp(y_global - RadiusRow:y_global + RadiusRow, x_global - RadiusCol:x_global + RadiusCol);

                MasterValue = MasterValue(SHP.PixelInd(:, num));
                SlaveValue = SlaveValue(SHP.PixelInd(:, num));
                InterfValue = InterfValue(SHP.PixelInd(:, num));

                nu(kk, jj) = mean(InterfValue);
                de1(kk, jj) = mean(MasterValue);
                de2(kk, jj) = mean(SlaveValue);
                num = num + 1;
            end
        end

        coh = nu ./ sqrt(de1 .* de2);
        coh(isnan(coh)) = 0;

        if strcmpi(BiasCorr, 'y')
            tmp = padarray(coh, [RadiusRow RadiusCol], 'symmetric');
            num = 1;
            for jj = 1:nwidths
                for kk = 1:nlines
                    x_global = jj + RadiusCol;
                    y_global = kk + RadiusRow;

                    CohValue = tmp(y_global - RadiusRow:y_global + RadiusRow, x_global - RadiusCol:x_global + RadiusCol);
                    CohValue = CohValue(SHP.PixelInd(:, num));

                    coh(kk, jj) = exp(mean(log(CohValue)));
                    num = num + 1;
                end
            end
        end

        % coh = single(coh);
        pairName = sprintf('%06d_%06d.mat', m, s);
        save(fullfile(outDir, pairName), 'coh', '-v7.3');

        if opts.Verbose
            fprintf('Saved AdpCoh for %06d-%06d to %s\n', m, s, outDir);
        end
    end

    if opts.Verbose
        disp('AdpCohEst_StreamByPair completed.');
    end
end
