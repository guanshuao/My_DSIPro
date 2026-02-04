clear;
clc;

window_sizes = [3, 11, 21]; % 尝试不同的窗口尺寸
rho = linspace(0, 1, 101); 
number_of_Monte_Carlo = 100000; % 蒙特卡洛试验次数
colors = {'r', 'g', 'b'}; % 曲线对应的颜色

rho_hat_abs_mean = zeros(numel(window_sizes), numel(rho));
rho_hat_abs_std = zeros(numel(window_sizes), numel(rho));
rho_bias = zeros(numel(window_sizes), numel(rho));
rho_rmse = zeros(numel(window_sizes), numel(rho));

for w_idx = 1:numel(window_sizes)
    window_size = window_sizes(w_idx);
    number_of_pixels = window_size^2;

    for r_idx = 1:numel(rho)
        current_rho = rho(r_idx);

        A = randn(number_of_Monte_Carlo, number_of_pixels) + 1i * randn(number_of_Monte_Carlo, number_of_pixels);
        B = current_rho * A + sqrt(1 - current_rho^2) * ...
            (randn(number_of_Monte_Carlo, number_of_pixels) + 1i * randn(number_of_Monte_Carlo, number_of_pixels));

        A_mean = mean(A, 2);
        B_mean = mean(B, 2);
        A_centered = A - A_mean;
        B_centered = B - B_mean;

        numerator = sum(A_centered .* conj(B_centered), 2);
        denom = sqrt(sum(abs(A_centered).^2, 2) .* sum(abs(B_centered).^2, 2));
        denom(denom == 0) = eps;

        rho_vals = numerator ./ denom;
        abs_vals = abs(rho_vals);

        rho_hat_abs_mean(w_idx, r_idx) = mean(abs_vals);
        rho_hat_abs_std(w_idx, r_idx) = std(abs_vals, 0, 1);
        rho_bias(w_idx, r_idx) = rho_hat_abs_mean(w_idx, r_idx) - current_rho;
        rho_rmse(w_idx, r_idx) = sqrt(mean((abs_vals - current_rho).^2));
    end
end

% 绘制不同窗口尺寸下的平均估计曲线
figure;
hold on;
plot(rho, rho, '--k', 'LineWidth', 1.2, 'DisplayName', '$|\gamma|$');

for w_idx = 1:numel(window_sizes)
    mean_curve = rho_hat_abs_mean(w_idx, :);
    std_curve = rho_hat_abs_std(w_idx, :);
    lower = max(mean_curve - std_curve, 0);
    upper = min(mean_curve + std_curve, 1);

    fill([rho, fliplr(rho)], [lower, fliplr(upper)], colors{w_idx}, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(rho, mean_curve, 'LineWidth', 2, 'Color', colors{w_idx}, ...
        'DisplayName', sprintf('%d x %d', window_sizes(w_idx), window_sizes(w_idx)));
end

xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('$|\hat{\gamma}|$', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'southeast');
title('Monte Carlo Calibration of $|\hat{\gamma}|$ with $\pm1\sigma$ Envelopes', 'Interpreter', 'latex');
grid on;
xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]);
hold off;