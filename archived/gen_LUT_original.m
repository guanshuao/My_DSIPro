% 该脚本运行速度差，数值稳定性差，已停用

clear;clc;close all;
n_values = 1:900;
rho_values = 0:0.001:1;
rho_squared = rho_values.^2;
tail_base = max(1 - rho_squared, 0);  % 避免浮点误差导致负数

num_n = numel(n_values);
num_rho = numel(rho_values);

mean_rho_hat = zeros(num_n, num_rho);
std_rho_hat = zeros(num_n, num_rho);

log_gamma_three_half = gammaln(1.5);
log_gamma_two = gammaln(2);  % 等于0，但保留写法以便阅读

parfor idx_n = 1:num_n
    disp(idx_n);
	n = n_values(idx_n);
	log_gamma_n = gammaln(n);
	log_gamma_n_plus_half = gammaln(n + 0.5);
	log_gamma_n_plus_one = gammaln(n + 1);

	% tail_term = tail_base.^n;
	log_tail_term = n * log(tail_base);

	mean_hyper = hypergeom([1.5, n, n], [1, n + 0.5], rho_squared);
	mean_row = exp(log_gamma_n + log_gamma_three_half - log_gamma_n_plus_half) .* mean_hyper .* exp(log_tail_term);
	mean_row = real(mean_row);  % 去除可能出现的微小虚部

	second_hyper = hypergeom([2, n, n], [1, n + 1], rho_squared);
	second_moment_row = exp(log_gamma_n + log_gamma_two - log_gamma_n_plus_one) .* second_hyper .* exp(log_tail_term);
	second_moment_row = real(second_moment_row);

	variance_row = max(second_moment_row - mean_row.^2, 0);
	std_row = sqrt(variance_row);

	mean_rho_hat(idx_n, :) = mean_row;
	std_rho_hat(idx_n, :) = std_row;
end

LUT = struct( ...
	'n_values', n_values, ...
	'rho_values', rho_values, ...
	'mean_rho_hat', mean_rho_hat, ...
	'std_rho_hat', std_rho_hat);

save('LUT.mat', 'LUT', '-v7.3');

fprintf('LUT生成完成：%d个n点 × %d个rho点，结果已保存到LUT.mat\n', num_n, num_rho);