function [varargout] = LUT_lookup(direction, n, value, lut_path)
% LUT_LOOKUP  基于预生成LUT进行正向或反向查找
%
% 正向查找 (forward): 已知样本数n和真实相干系数rho，查找rho_hat期望和标准差
%   rho_hat = LUT_lookup('forward', n, rho)
%   [rho_hat, std] = LUT_LOOKUP('forward', n, rho)
%   [rho_hat, std] = LUT_LOOKUP('forward', n, rho, lut_path)
%
% 反向查找 (inverse): 已知样本数n和估计相干系数rho_hat，反推真实相干系数rho
%   rho = LUT_LOOKUP('inverse', n, rho_hat)
%   rho = LUT_LOOKUP('inverse', n, rho_hat, lut_path)
%
% 输入:
%   direction - 查找方向: 'forward' 或 'inverse' (必需)
%   n         - 样本数，必须为整数且在LUT范围内 (必需)
%   value     - 正向时为rho，反向时为rho_hat (必需)
%   lut_path  - LUT文件路径 (可选，默认 'LUT.mat')
%
% 输出:
%   正向: rho_hat (期望), std (标准差，可选)
%   反向: rho (真实相干系数)

    %% 参数验证
    if nargin < 3
        error('LUT_lookup:NotEnoughInputs', '至少需要3个输入参数: direction, n, value。');
    end
    if nargin < 4 || isempty(lut_path)
        lut_path = 'LUT.mat';
    end
    
    direction = validateDirection(direction);
    
    %% 加载LUT（带缓存）
    LUT = loadLUT(lut_path);
    log_path = resolveLogPath(lut_path);
    logMessage(log_path, sprintf('LUT_lookup start: direction=%s, n_size=%s, value_size=%s', ...
        direction, mat2str(size(n)), mat2str(size(value))));
    
    %% 验证并预处理输入
    n = double(n);
    value = double(value);
    
    validateN(n, LUT.n_values);
    [n_idx, value_vec, output_size] = prepareInputs(n, value, LUT.n_values);
    
    %% 执行查找
    if strcmp(direction, 'forward')
        validateRho(value, LUT.rho_values);
        [result1, result2] = forwardLookup(LUT, n_idx, value_vec, output_size);
        varargout{1} = result1;
        if nargout >= 2
            varargout{2} = result2;
        end
    else
        result = inverseLookup(LUT, n_idx, value_vec, output_size, log_path);
        varargout{1} = result;
        if nargout > 1
            warning('LUT_lookup:ExtraOutput', '反向查找仅返回一个输出值rho。');
        end
    end
    logMessage(log_path, 'LUT_lookup end.');
end

%% ========================================================================
%  输入验证函数
%  ========================================================================

function direction = validateDirection(direction)
    if ~ischar(direction) && ~isstring(direction)
        error('LUT_lookup:InvalidDirection', 'direction必须为字符串。');
    end
    direction = lower(char(direction));
    if ~ismember(direction, {'forward', 'inverse'})
        error('LUT_lookup:InvalidDirection', 'direction必须为 ''forward'' 或 ''inverse''。');
    end
end

function validateN(n, n_values)
    if any(n(:) ~= round(n(:)))
        error('LUT_lookup:InvalidN', 'n必须为整数。');
    end
    if any(n(:) < n_values(1)) || any(n(:) > n_values(end))
        error('LUT_lookup:OutOfRange', 'n超出LUT支持范围[%d, %d]。', n_values(1), n_values(end));
    end
end

function validateRho(rho, rho_values)
    if any(rho(:) < rho_values(1)) || any(rho(:) > rho_values(end))
        error('LUT_lookup:OutOfRange', 'rho超出LUT支持范围[%g, %g]。', rho_values(1), rho_values(end));
    end
end

%% ========================================================================
%  LUT加载函数（带缓存）
%  ========================================================================

function LUT = loadLUT(lut_path)
    persistent cached_path cached_LUT
    
    if isempty(cached_LUT) || ~strcmp(cached_path, lut_path)
        data = load(lut_path, 'LUT');
        if ~isfield(data, 'LUT')
            error('LUT_lookup:InvalidFile', '文件%s中缺少LUT结构体。', lut_path);
        end
        
        required = {'n_values', 'rho_values', 'mean_rho_hat', 'std_rho_hat'};
        missing = required(~isfield(data.LUT, required));
        if ~isempty(missing)
            error('LUT_lookup:InvalidStruct', 'LUT结构缺少字段: %s。', strjoin(missing, ', '));
        end
        
        cached_LUT = data.LUT;
        cached_path = lut_path;
    end
    
    LUT = cached_LUT;
end

%% ========================================================================
%  输入预处理函数
%  ========================================================================

function [n_idx, value_vec, output_size] = prepareInputs(n, value, n_values)
    % 确定输出尺寸
    if isscalar(n) && isscalar(value)
        output_size = [1, 1];
    elseif isscalar(n)
        output_size = size(value);
    elseif isscalar(value)
        output_size = size(n);
    else
        if ~isequal(size(n), size(value))
            error('LUT_lookup:SizeMismatch', 'n与value必须同尺寸或有一个为标量。');
        end
        output_size = size(n);
    end
    
    % 向量化并扩展标量
    n_vec = n(:);
    value_vec = value(:);
    num_elements = max(numel(n_vec), numel(value_vec));
    
    if isscalar(n_vec)
        n_vec = repmat(n_vec, num_elements, 1);
    end
    if isscalar(value_vec)
        value_vec = repmat(value_vec, num_elements, 1);
    end
    
    % 计算n的索引（LUT行号）
    n_idx = n_vec - n_values(1) + 1;
end

%% ========================================================================
%  查找函数
%  ========================================================================

function [rho_hat, rho_std] = forwardLookup(LUT, n_idx, rho_vec, output_size)
% 正向查找: 给定n和rho，返回rho_hat期望和标准差

    num_pts = numel(n_idx);
    rho_hat_vec = zeros(num_pts, 1);
    rho_std_vec = zeros(num_pts, 1);
    
    % 按唯一n值分组处理，减少重复提取
    [unique_idx, ~, ic] = unique(n_idx);
    
    for k = 1:numel(unique_idx)
        idx = unique_idx(k);
        mask = (ic == k);
        rho_hat_vec(mask) = interp1(LUT.rho_values, LUT.mean_rho_hat(idx, :), rho_vec(mask), 'linear');
        rho_std_vec(mask) = interp1(LUT.rho_values, LUT.std_rho_hat(idx, :), rho_vec(mask), 'linear');
    end
    
    rho_hat = reshape(rho_hat_vec, output_size);
    rho_std = reshape(rho_std_vec, output_size);
end

function rho = inverseLookup(LUT, n_idx, rho_hat_vec, output_size, log_path)
% 反向查找: 给定n和rho_hat，返回真实rho

    num_pts = numel(n_idx);
    rho_vec = zeros(num_pts, 1);
    
    % 按唯一n值分组处理
    [unique_idx, ~, ic] = unique(n_idx);
    
    for k = 1:numel(unique_idx)
        idx = unique_idx(k);
        mask = (ic == k);
        actual_n = idx + LUT.n_values(1) - 1;
        
        rho_hat_k = rho_hat_vec(mask);
        num_total = numel(rho_hat_k);
        mean_curve = LUT.mean_rho_hat(idx, :);
        valid_min = mean_curve(1);
        valid_max = mean_curve(end);
        
        % n=1时不做任何限制，全部保留原值
        if actual_n == 1
            logMessage(log_path, sprintf('n=1: keep all %d values unchanged.', num_total));
            rho_vec(mask) = rho_hat_k;
            continue;
        end
        
        % n>1时：超出范围的值进行截断
        out_low = rho_hat_k < valid_min;
        out_high = rho_hat_k > valid_max;
        num_low = sum(out_low);
        num_high = sum(out_high);
        if num_low > 0 || num_high > 0
            logMessage(log_path, sprintf('n=%d: %d below min->0, %d above max->1, total=%d (range=[%.4f, %.4f]).', ...
                actual_n, num_low, num_high, num_total, valid_min, valid_max));
        end
        
        % 反向插值：范围内的值插值，超出范围的置0或1
        rho_result_k = rho_hat_k;
        rho_result_k(out_low) = 0;
        rho_result_k(out_high) = 1;
        in_range = ~(out_low | out_high);
        if any(in_range)
            rho_result_k(in_range) = interp1(mean_curve, LUT.rho_values, rho_hat_k(in_range), 'linear');
        end
        rho_vec(mask) = rho_result_k;
    end
    
    rho = reshape(rho_vec, output_size);
end

%% ========================================================================
%  日志函数
%  ========================================================================

function log_path = resolveLogPath(lut_path)
    [folder, ~, ~] = fileparts(lut_path);
    if isempty(folder)
        folder = pwd;
    end
    log_path = fullfile(folder, 'LUT_lookup.log');
end

function logMessage(log_path, msg)
    try
        fid = fopen(log_path, 'a');
        if fid == -1
            return;
        end
        timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        fprintf(fid, '[%s] %s\n', timestamp, msg);
        fclose(fid);
    catch
        % 忽略日志写入错误
    end
end