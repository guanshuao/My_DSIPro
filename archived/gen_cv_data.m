function X = gen_cv_data(sigma, m, n, p, varargin)
%GEN_CV_DATA 生成正数、双峰且指定变异系数的数据。
%
%   X = GEN_CV_DATA(sigma, m, n, p) 返回一个 m×n×p 的数组 X。
%   生成的数据满足：
%   - 所有取值严格为正数
%   - 分布近似双峰（两组对数正态混合）
%   - 样本均值约为 1（通过归一化实现，浮点误差下可视为 1）
%   - 样本标准差约为 sigma（使用 std(.,0)），因此变异系数 CV=std/mean≈sigma
%
%   可选参数（Name-Value）：
%   - 'Weight'   : 第一峰所占比例 w∈(0,1)，默认 0.5
%   - 'LogSigma' : 对数域标准差（峰宽），默认 0.25
%   - 'Delta'    : 两峰在对数域的间隔（越大越双峰），默认 0.8
%   - 'Seed'     : 随机种子（空则不设置），默认 []
%   - 'MaxIter'  : 二分迭代次数上限，默认 60
%   - 'Tol'      : 目标 sigma 的相对/绝对容差，默认 1e-10
%
%   说明：
%   1) 先生成正数双峰基样本 x（两组对数正态混合）；
%   2) 再做幂变换 x^k，通过二分搜索选取 k，使得 CV(x^k)=sigma；
%   3) 最后归一化到均值 1。

    validateattributes(sigma, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'sigma');
    validateattributes(m, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'm');
    validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'n');
    validateattributes(p, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'p');

    ip = inputParser;
    ip.FunctionName = mfilename;
    ip.CaseSensitive = false;
    addParameter(ip, 'Weight', 0.5, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0 && v < 1);
    addParameter(ip, 'LogSigma', 0.25, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    addParameter(ip, 'Delta', 0.8, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 0);
    addParameter(ip, 'Seed', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && isfinite(v)));
    addParameter(ip, 'MaxIter', 60, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1);
    addParameter(ip, 'Tol', 1e-10, @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v > 0);
    parse(ip, varargin{:});
    opts = ip.Results;

    if ~isempty(opts.Seed)
        rng(opts.Seed);
    end

    N = double(m) * double(n) * double(p);

    if sigma == 0
        X = ones(m, n, p, 'double');
        return;
    end

    % ---- 第 1 步：生成正数双峰基样本（两组对数正态混合）
    % 可通过 opts.Delta / opts.LogSigma / opts.Weight 控制双峰形态
    w = opts.Weight;         % 混合权重
    delta = opts.Delta;      % 对数域两峰间隔

    N1 = max(1, round(w * N));
    N2 = max(1, N - N1);

    x1 = exp(-delta + opts.LogSigma * randn(N1, 1));
    x2 = exp( delta + opts.LogSigma * randn(N2, 1));

    x = [x1; x2];
    if numel(x) ~= N
        x = x(1:N);
    end

    % 打乱顺序，使两类样本充分混合
    x = x(randperm(N));

    % ---- 第 2 步：选择幂指数 k，使得 CV(x^k) == sigma
    lx = log(x);

    % 对于非退化的正数样本，CV(k) 通常随 k 增大而增大，可用二分法求解
    kLow = 0.0;
    kHigh = 1.0;

    fHigh = local_cv_from_log(lx, kHigh) - sigma;

    % 扩展上界，直到可以达到目标 sigma
    maxK = 1e6;
    while fHigh < 0 && kHigh < maxK
        kHigh = kHigh * 2;
        fHigh = local_cv_from_log(lx, kHigh) - sigma;
    end

    if fHigh < 0
        error('gen_cv_data:UnachievableSigma', ...
            '无法达到 sigma=%g；建议增大样本量 m*n*p，或减小 sigma，或调整 ''Delta''/''LogSigma''。', sigma);
    end

    % 二分法求解 k
    k = 0.5 * (kLow + kHigh);
    for iter = 1:round(opts.MaxIter)
        k = 0.5 * (kLow + kHigh);
        fMid = local_cv_from_log(lx, k) - sigma;

        if abs(fMid) <= max(1e-12, opts.Tol * max(1, sigma))
            break;
        end

        if fMid > 0
            kHigh = k;
        else
            kLow = k;
        end
    end

    % ---- 第 3 步：应用变换并归一化到均值 1
    % 数值稳定处理：y = exp(k*lx)/mean(exp(k*lx))
    % 先减去最大值避免溢出；比例不变
    z = k * lx;
    a = max(z);
    v = exp(z - a);
    mu = mean(v);
    if mu <= 0 || ~isfinite(mu)
        error('gen_cv_data:NumericalIssue', '归一化时出现数值问题（mean=%g）。', mu);
    end
    y = v ./ mu;

    % Reshape to requested output size
    X = reshape(y, m, n, p);
end

function cv = local_cv_from_log(lx, k)
% 在对数域计算 CV：CV(exp(k*lx)) = std/mean。
% 为了数值稳定，使用 exp(z-max(z)) 缩放，CV 不受缩放影响。
    z = k * lx;
    a = max(z);
    v = exp(z - a);
    mu = mean(v);
    if mu == 0 || ~isfinite(mu)
        cv = Inf;
        return;
    end
    cv = std(v, 0) ./ mu;
end
