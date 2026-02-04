% 本程序通过蒙特卡洛模拟相干性估计的偏差和方差，通过选取一定比例的“被替换点”来模拟非均匀场景，对比均匀场景与非均匀场景下的估计性能。
clear;
clc;

window_size = 3;number_of_pixels = window_size^2; % 相干性估计的窗口大小、窗口内像素总数
rho = linspace(0, 1, 101); % 真实相关系数范围
number_of_MC_trials = 100000;% 蒙特卡洛试验次数，一般取10万，需要精确绘制时可取100万
cv_list = [0.7 1.0 1.3]; % 非均一场景的 cv 取值
% 为了避免占用巨大的三维内存（M1/M2/X 等），这里改为按 rho(k) 逐点生成样本并即时统计
n_rho = length(rho);
n_cv = length(cv_list);

rho_hat_abs_mean_homo = zeros(1, n_rho);
rho_hat_abs_std_homo = zeros(1, n_rho);
rho_rmse_homo = zeros(1, n_rho);

rho_hat_abs_mean_hetero = zeros(n_cv, n_rho);
rho_hat_abs_std_hetero = zeros(n_cv, n_rho);
rho_rmse_hetero = zeros(n_cv, n_rho);

fprintf('Running Monte Carlo (homo + hetero for cv = [%s])...\n', num2str(cv_list));
for k = 1:n_rho
    fprintf('rho progress: %d/%d\n', k, n_rho);

    r = rho(k);
    A = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
    B = r * A + sqrt(1 - r^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));

    % 均一场景
    abs_rho_hat = sample_abs_corr(A, B);
    rho_hat_abs_mean_homo(k) = mean(abs_rho_hat);
    rho_hat_abs_std_homo(k) = std(abs_rho_hat, 0);
    rho_rmse_homo(k) = sqrt(mean((abs_rho_hat - r).^2));

    % 非均一场景：对 A/B 相同位置进行替换/置零（取决于 gen_cv_data 的定义）
    for icv = 1:n_cv
        X = gen_cv_data(cv_list(icv), number_of_MC_trials, number_of_pixels, 1);
        if ndims(X) == 3
            X = X(:, :, 1);
        end

        abs_rho_hat_r = sample_abs_corr(X .* A, X .* B);
        rho_hat_abs_mean_hetero(icv, k) = mean(abs_rho_hat_r);
        rho_hat_abs_std_hetero(icv, k) = std(abs_rho_hat_r, 0);
        rho_rmse_hetero(icv, k) = sqrt(mean((abs_rho_hat_r - r).^2));
    end
end

rho_bias_homo = rho_hat_abs_mean_homo - rho;
rho_bias_hetero = rho_hat_abs_mean_hetero - rho;

% 绘制平均估计值与真值对比（均一 + 不同 cv）
figure;
hold on;
plot(rho, rho, '--k', 'LineWidth', 1.2, 'DisplayName', '$|\gamma|$');
colors = lines(n_cv + 1);
plot(rho, rho_hat_abs_mean_homo, 'Color', colors(1, :), 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{homo}$');
for icv = 1:n_cv
    lbl = sprintf('$|\\hat{\\gamma}|_{cv=%.1f}$', cv_list(icv));
    plot(rho, rho_hat_abs_mean_hetero(icv, :), 'Color', colors(icv + 1, :), 'LineWidth', 2, 'DisplayName', lbl);
end


grid on;
xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]); % 保证横纵坐标比例一致，不受窗口拉伸影响
legend('Interpreter', 'latex', 'Location', 'southeast');
hold off;

% 绘制标准差对比（均一 + 不同 cv）
figure;
hold on;
plot(rho, rho_hat_abs_std_homo, 'Color', colors(1, :), 'LineWidth', 2, 'DisplayName', '$\sigma_{homo}$');
for icv = 1:n_cv
    lbl = sprintf('$\\sigma_{cv=%.1f}$', cv_list(icv));
    plot(rho, rho_hat_abs_std_hetero(icv, :), 'Color', colors(icv + 1, :), 'LineWidth', 2, 'DisplayName', lbl);
end


grid on;
xlim([0 1]);
pbaspect([1 1 1]);
legend('Interpreter', 'latex', 'Location', 'northeast');
hold off;


function abs_rho_hat = sample_abs_corr(A, B)
% 计算每一行（一个 Monte Carlo 试验）在像素维度上的复相关系数幅度 |rho_hat|
    A_centered = A - mean(A, 2);
    B_centered = B - mean(B, 2);

    numerator = sum(A_centered .* conj(B_centered), 2);
    denom = sqrt(sum(abs(A_centered).^2, 2) .* sum(abs(B_centered).^2, 2));
    denom(denom == 0) = eps;

    rho_hat = numerator ./ denom;
    abs_rho_hat = abs(rho_hat);
end