% GEN_LUT  生成样本相干系数统计特性的查找表 (LUT)
%
% 功能:
%   基于广义超几何函数理论，计算不同样本数 (n) 和真实相干系数 (rho) 下，样本相干系数 (rho_hat) 的期望值和标准差。
%
% 理论公式:
%   E[rho_hat]   = C1 * (1-rho²)^n * 3F2(1.5, n, n; 1, n+0.5; rho²)
%   E[rho_hat²]  = C2 * (1-rho²)^n * 3F2(2, n, n; 1, n+1; rho²)
%   其中 C1 = Gamma(n)*Gamma(1.5)/Gamma(n+0.5), C2 = Gamma(n)*Gamma(2)/Gamma(n+1) = 1/n
%
% 数值方法:
%   使用对数空间级数展开，将 (1-rho²)^n 与超几何级数合并计算，
%   避免大 n 或 rho 接近 1 时的数值溢出问题。
%
% 输出:
%   LUT.mat - 包含结构体 LUT:
%     n_values     - 窗口尺寸向量
%     rho_values   - 真实相干系数向量
%     mean_rho_hat - 期望矩阵 [num_n × num_rho]
%     std_rho_hat  - 标准差矩阵 [num_n × num_rho]
%
% 用法: 直接运行 gen_LUT
% -------------------------------------------------------------------------

%% 参数设置
clear; clc; close all;

N = 841;  % 最大样本数
RhoStep = 0.01;  % 相干系数步长

n_values   = 1:N;
rho_values = 0:RhoStep:1;
num_n      = numel(n_values);
num_rho    = numel(rho_values);

mean_rho_hat = zeros(num_n, num_rho);
std_rho_hat  = zeros(num_n, num_rho);

%% 主计算循环

for idx_n = 1:num_n
    if mod(idx_n, 50) == 0
        fprintf('Processing n = %d\n', idx_n);
    end
    
    n = n_values(idx_n);
    log_C1 = gammaln(n) + gammaln(1.5) - gammaln(n + 0.5);
    log_C2 = gammaln(n) + gammaln(2)   - gammaln(n + 1);
    
    mean_row = zeros(1, num_rho);
    std_row  = zeros(1, num_rho);
    
    for idx_rho = 1:num_rho
        rho = rho_values(idx_rho);
        rho_sq = rho^2;
        
        % 边界情况: rho = 1
        if rho >= 1 - 1e-12
            mean_row(idx_rho) = 1;
            std_row(idx_rho)  = 0;
            continue;
        end
        
        % 计算一阶矩和二阶矩
        mean_val      = compute_hypergeom_product(n, rho_sq, 1.5, n+0.5, log_C1);
        second_moment = compute_hypergeom_product(n, rho_sq, 2,   n+1,   log_C2);
        
        % 计算标准差
        variance = max(second_moment - mean_val^2, 0);
        
        mean_row(idx_rho) = mean_val;
        std_row(idx_rho)  = sqrt(variance);
    end
    
    mean_rho_hat(idx_n, :) = mean_row;
    std_rho_hat(idx_n, :)  = std_row;
end

%% 保存结果
% 将N和RhoStep信息包含在info字段中，先生成字符串以避免解析问题
info_str = sprintf('The range of n is 1 to %d, RhoStep is %.4f', N, RhoStep);
LUT = struct( ...
    'info',         info_str, ...
    'n_values',     n_values, ...
    'rho_values',   rho_values, ...
    'mean_rho_hat', mean_rho_hat, ...
    'std_rho_hat',  std_rho_hat);

save('LUT.mat', 'LUT', '-v7.3');
fprintf('LUT生成完成: %d x %d, 已保存到 LUT.mat\n', num_n, num_rho);

%% ========================================================================
%  辅助函数
%  ========================================================================

function result = compute_hypergeom_product(n, z, a1, b2, log_C)
% COMPUTE_HYPERGEOM_PRODUCT  计算 C * (1-z)^n * 3F2(a1,n,n; 1,b2; z)
%
% 使用对数空间在线累加级数，处理 z 接近 1 时的数值稳定性问题。
% 级数展开: 3F2 = sum_{k=0}^{inf} [(a1)_k (n)_k^2 / (1)_k (b2)_k] * z^k / k!

    % 特殊情况: z ≈ 0
    if z < 1e-15
        result = exp(log_C + n * log(1 - z));
        return;
    end
    
    % 预计算
    log_z = log(z);
    log_one_minus_z = log(1 - z);
    
    % 初始化: k=0 项
    log_term = log_C + n * log_one_minus_z;
    log_sum  = log_term;
    
    % 迭代参数
    max_iter = (z > 0.9) * 80000 + 20000;  % z > 0.9 时用 100000, 否则 20000
    tol = 1e-14;
    converge_count = 0;
    
    for k = 0:max_iter-1
        % 递推: term_{k+1} = term_k * (a1+k)(n+k)^2 / [(1+k)(b2+k)(k+1)] * z
        log_ratio = log(a1 + k) + 2*log(n + k) ...
                  - log(1 + k) - log(b2 + k) - log(k + 1) + log_z;
        log_term = log_term + log_ratio;
        
        % 对数空间稳定累加
        diff = log_term - log_sum;
        if diff > 0
            log_sum = log_term + log1p(exp(-diff));
        else
            log_sum = log_sum + log1p(exp(diff));
        end
        
        % 收敛判断: 连续 5 次相对贡献 < tol
        if diff < log(tol)
            converge_count = converge_count + 1;
            if converge_count >= 5, break; end
        else
            converge_count = 0;
        end
    end
    
    result = exp(log_sum);
    
    % 边界约束
    if ~isfinite(result), result = 0; end
    result = min(max(result, 0), 1);
end