clear;clc;close all;

%% 参数设置
coh_window_size = 3; % 相干性估计的窗口大小
number_of_pixels = coh_window_size^2; % 窗口内像素总数
lee_window_size = 7; % Lee滤波的窗口大小
rho = linspace(0, 1, 101); % 真实相关系数取值范围
number_of_MC_trials = 10000;% 蒙特卡洛试验次数
cv_list = [0.4 0.7 1.0 2.0 4];  % CV值列表
num_cv = length(cv_list);         % CV值数量
cv_tolerance = 0.02;      % CV浮动范围
cv_window_size = 21;      % 计算CV的窗口大小

img = SingleRead('F:\Coh_Est\SL\MLI\20220518.rslc',1201,'cpxfloat32');
img = flip(img);
power = (abs(img)).^2;

colors = [
    0.0 0.0 1.0;   % 蓝色 - homo
    1.0 0.0 0.0;   % 红色 - cv1
    0.0 0.8 0.0;   % 绿色 - cv2
    1.0 0.5 0.0;   % 橙色 - cv3
    0.5 0.0 0.5;   % 紫色 - cv4
    0.0 0.8 0.8;   % 青色 - cv5
];

%% 初始化矩阵
M1 = zeros(number_of_MC_trials, number_of_pixels, length(rho));
M2 = zeros(number_of_MC_trials, number_of_pixels, length(rho));
rho_MC = zeros(number_of_MC_trials, length(rho));
abs_rho_MC = zeros(number_of_MC_trials, length(rho));
% 为每个CV值初始化存储矩阵
rho_MC_replace = zeros(number_of_MC_trials, length(rho), num_cv);
abs_rho_MC_replace = zeros(number_of_MC_trials, length(rho), num_cv);

%% 生成蒙特卡洛样本
fprintf('Generating Monte Carlo samples...\n');
parfor k = 1:length(rho)
    fprintf('Progress: %d/%d\n', k, length(rho));
    M1(:, :, k) = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
    M2(:, :, k) = rho(k) * M1(:, :, k) + sqrt(1 - (rho(k))^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));
end
fprintf('Monte Carlo samples generated.\n');

%% 对原始样本计算样本相关系数
fprintf('Calculating sample correlation coefficients for original samples...\n');
parfor k = 1:length(rho)
    fprintf('Progress: %d/%d\n', k, length(rho));
    A = M1(:, :, k);
    B = M2(:, :, k);
    A_mean = mean(A, 2);
    B_mean = mean(B, 2);
    A_centered = A - A_mean;
    B_centered = B - B_mean;
    numerator = sum(A_centered .* conj(B_centered), 2);
    denom = sqrt(sum(abs(A_centered).^2, 2) .* sum(abs(B_centered).^2, 2));
    denom(denom == 0) = eps;
    rho_MC(:, k) = numerator ./ denom;
    abs_rho_MC(:, k) = abs(rho_MC(:, k));
end

%% 预先计算Lee滤波
fprintf('Applying Lee filter...\n');
img_lee = LeeFilter(power, lee_window_size, 1);

%% 对每个CV值生成替换矩阵并计算样本相关系数
for cv_idx = 1:num_cv
    target_cv = cv_list(cv_idx);
    fprintf('\n=== Processing CV = %.2f (%d/%d) ===\n', target_cv, cv_idx, num_cv);
    
    % 生成替换矩阵
    fprintf('Generating replacement matrix based on CV filtering...\n');
    sigma = find_pixels_by_cv(img_lee, target_cv, cv_tolerance, cv_window_size);
    texture = gen_iid_matrix(sqrt(sigma), number_of_MC_trials, number_of_pixels, length(rho));
    M1_replace = texture.*M1;
    M2_replace = texture.*M2;
    
    % 对替换后的样本计算样本相关系数
    fprintf('Calculating sample correlation coefficients for replaced samples...\n');
    temp_rho_MC_replace = zeros(number_of_MC_trials, length(rho));
    temp_abs_rho_MC_replace = zeros(number_of_MC_trials, length(rho));
    
    parfor k = 1:length(rho)
        A_replace = M1_replace(:, :, k);
        B_replace = M2_replace(:, :, k);
        A_replace_mean = mean(A_replace, 2);
        B_replace_mean = mean(B_replace, 2);
        A_replace_centered = A_replace - A_replace_mean;
        B_replace_centered = B_replace - B_replace_mean;
        numerator_replace = sum(A_replace_centered .* conj(B_replace_centered), 2);
        denom_replace = sqrt(sum(abs(A_replace_centered).^2, 2) .* sum(abs(B_replace_centered).^2, 2));
        denom_replace(denom_replace == 0) = eps;
        temp_rho_MC_replace(:, k) = numerator_replace ./ denom_replace;
        temp_abs_rho_MC_replace(:, k) = abs(temp_rho_MC_replace(:, k));
    end
    
    rho_MC_replace(:, :, cv_idx) = temp_rho_MC_replace;
    abs_rho_MC_replace(:, :, cv_idx) = temp_abs_rho_MC_replace;
end
clear M1 M2 M1_replace M2_replace texture;

%% 计算样本相关系数幅度的统计量
rho_hat_abs_mean = mean(abs_rho_MC, 1); 
rho_hat_abs_std = std(abs_rho_MC, 0, 1);

% 为每个CV值计算统计量
rho_hat_abs_mean_replace = zeros(num_cv, length(rho));
rho_hat_abs_std_replace = zeros(num_cv, length(rho));
for cv_idx = 1:num_cv
    rho_hat_abs_mean_replace(cv_idx, :) = mean(abs_rho_MC_replace(:, :, cv_idx), 1);
    rho_hat_abs_std_replace(cv_idx, :) = std(abs_rho_MC_replace(:, :, cv_idx), 0, 1);
end

%% 绘制平均估计值与真值对比（不含标准差范围带）
figure;
hold on;
plot(rho, rho, '--k', 'LineWidth', 1.2, 'DisplayName', '$|\gamma|$');
plot(rho, rho_hat_abs_mean, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{homo}$');
for cv_idx = 1:num_cv
    plot(rho, rho_hat_abs_mean_replace(cv_idx, :), 'Color', colors(cv_idx+1,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('$|\\hat{\\gamma}|_{CV=%.1f}$', cv_list(cv_idx)));
end
xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('$|\hat{\gamma}|$', 'Interpreter', 'latex');
title(sprintf('Mean Coherence Estimate (Window Size: %d)', coh_window_size), 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'southeast');
grid on;
xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]);
hold off;

%%  绘制标准差对比
figure;
hold on;
plot(rho, rho_hat_abs_std, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '$\sigma_{homo}$');
for cv_idx = 1:num_cv
    plot(rho, rho_hat_abs_std_replace(cv_idx, :), 'Color', colors(cv_idx+1,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('$\\sigma_{CV=%.1f}$', cv_list(cv_idx)));
end
xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('Standard Deviation', 'Interpreter', 'latex');
title(sprintf('Standard Deviation of Coherence Estimate (Window Size: %d)', coh_window_size), 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northeast');
grid on;
xlim([0 1]);
pbaspect([1 1 1]);
hold off;

%%  绘制平均估计值与真值对比（含标准差范围带）
figure;
hold on;
plot(rho, rho, '--k', 'LineWidth', 1.2, 'DisplayName', '$|\gamma|$');
plot(rho, rho_hat_abs_mean, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{homo}$');

% 绘制homo的标准差范围
fill([rho, fliplr(rho)], [rho_hat_abs_mean - rho_hat_abs_std, fliplr(rho_hat_abs_mean + rho_hat_abs_std)], ...
    colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '$\pm 1 \sigma_{homo}$');

for cv_idx = 1:num_cv
    plot(rho, rho_hat_abs_mean_replace(cv_idx, :), 'Color', colors(cv_idx+1,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('$|\\hat{\\gamma}|_{CV=%.1f}$', cv_list(cv_idx)));
    % 绘制每个CV的标准差范围
    fill([rho, fliplr(rho)], ...
        [rho_hat_abs_mean_replace(cv_idx, :) - rho_hat_abs_std_replace(cv_idx, :), ...
         fliplr(rho_hat_abs_mean_replace(cv_idx, :) + rho_hat_abs_std_replace(cv_idx, :))], ...
        colors(cv_idx+1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'DisplayName', sprintf('$\\pm 1 \\sigma_{CV=%.1f}$', cv_list(cv_idx)));
end
xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('$|\hat{\gamma}|$', 'Interpreter', 'latex');
title(sprintf('Mean Coherence Estimate with $\\pm 1\\sigma$ (Window Size: %d)', coh_window_size), 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'southeast');
grid on;
xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]);
hold off;
