% 蒙特卡洛范例程序
% 本程序通过蒙特卡洛模拟相干性估计的偏差和方差，通过选取一定比例的“被替换点”来模拟非均匀场景，对比均匀场景与非均匀场景下的估计性能。
clear;clc;close all;

%% 参数设置
window_size = 3; % 相干性估计的窗口大小
number_of_pixels = window_size^2; % 窗口内像素总数
rho = linspace(0, 1, 101); % 真实相关系数取值范围
number_of_MC_trials = 100000;% 蒙特卡洛试验次数，一般取10万，需要精确绘制时可取100万
contamination_ratio = 0.5; % 被替换点占像素的比例，称为污染率
amplification_scale = 2; % 被替换点幅度与原始幅度的比值，称为amplification scale

%% 初始化矩阵
M1 = zeros(number_of_MC_trials, number_of_pixels, length(rho));
M2 = zeros(number_of_MC_trials, number_of_pixels, length(rho));
rho_MC = zeros(number_of_MC_trials, length(rho));
abs_rho_MC = zeros(number_of_MC_trials, length(rho));
rho_MC_replace = zeros(number_of_MC_trials, length(rho));
abs_rho_MC_replace = zeros(number_of_MC_trials, length(rho));

%% 生成蒙特卡洛样本
fprintf('Generating Monte Carlo samples...\n');
parfor k = 1:length(rho)
    % r = rho(k);
    % A = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
    % B = rho(k) * A + sqrt(1 - (rho(k))^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));
    % M1(:, :, k) = A;
    % M2(:, :, k) = B;
    % AB都是number_of_MC_trials行，number_of_pixels列的矩阵，每一行代表一次蒙特卡洛试验中生成的两个矩阵对应位置的number_of_pixels个像素值
    % A中的每个元素 randn + 1i*randn
    % B中的每个元素 rho*对应A中的元素 + sqrt(1-(rho)^2)*(randn + 1i*randn)
    fprintf('Progress: %d/%d\n', k, length(rho));
    M1(:, :, k) = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
    M2(:, :, k) = rho(k) * M1(:, :, k) + sqrt(1 - (rho(k))^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));
    % M1和M2是三维矩阵，共有101个切片，每个切片是number_of_MC_trials行，number_of_pixels列的矩阵
end
fprintf('Monte Carlo samples generated.\n');

%% 在整个三维样本中随机挑选 contamination_ratio 比例的元素并放大
M1_replace = M1;
M2_replace = M2;
fprintf('Selecting contaminated points...\n');
selected_idx = randperm(numel(M1), round(contamination_ratio * numel(M1))); % 随机选择索引
fprintf('Selecting contaminated points done.\n');
fprintf('Replacing contaminated points...\n');
M1_replace(selected_idx) = amplification_scale * M1_replace(selected_idx);
M2_replace(selected_idx) = amplification_scale * M2_replace(selected_idx);
clear selected_idx;
fprintf('Replacement of contaminated points done.\n');

% 每次计算 number_of_pixels 对样本的相关系数，生成一个 number_of_MC_trials 行，length(rho) 列的矩阵 RHO_MC
% RHO_MC(i, j) 表示第 i 次蒙特卡洛试验中，第 j 个 rho 值对应的两个矩阵的相关系数，每个元素都是复数
% n对样本的样本相关系数公式：\hat{\rho}=\frac{\sum_{i=1}^{n}(X_i-\bar{X})(Y_i-\bar{Y})^*}{\sqrt{\sum_{i=1}^{n}(X_i-\bar{X})^2\sum_{i=1}^{n}(Y_i-\bar{Y})^2}}


%% 对原始样本计算样本相关系数
fprintf('Calculating sample correlation coefficients for original samples...\n');
parfor k = 1:length(rho)
    fprintf('Progress: %d/%d\n', k, length(rho));
    % 提取当前rho对应的Monte Carlo样本
    A = M1(:, :, k);
    B = M2(:, :, k);

    % 沿像素维度计算样本均值并去中心化
    A_mean = mean(A, 2);
    B_mean = mean(B, 2);
    A_centered = A - A_mean;
    B_centered = B - B_mean;

    % 根据复相关系数定义计算分子和分母
    numerator = sum(A_centered .* conj(B_centered), 2);
    denom = sqrt(sum(abs(A_centered).^2, 2) .* sum(abs(B_centered).^2, 2));
    denom(denom == 0) = eps;

    % 计算样本相关系数及其幅度
    rho_MC(:, k) = numerator ./ denom;
    abs_rho_MC(:, k) = abs(rho_MC(:, k));
end
clear M1 M2 A B A_centered B_centered;

%% 对替换后的样本计算样本相关系数
fprintf('Calculating sample correlation coefficients for replaced samples...\n');
parfor k = 1:length(rho)
    fprintf('Progress: %d/%d\n', k, length(rho));
    A_replace = M1_replace(:, :, k);
    B_replace = M2_replace(:, :, k);
    A_replace_mean = mean(A_replace, 2);
    B_replace_mean = mean(B_replace, 2);
    A_replace_centered = A_replace - A_replace_mean;
    B_replace_centered = B_replace - B_replace_mean;
    numerator_replace = sum(A_replace_centered .* conj(B_replace_centered), 2);
    denom_replace = sqrt(sum(abs(A_replace_centered).^2, 2) .* sum(abs(B_replace_centered).^2, 2));
    denom_replace(denom_replace == 0) = eps;
    rho_MC_replace(:, k) = numerator_replace ./ denom_replace;
    abs_rho_MC_replace(:, k) = abs(rho_MC_replace(:, k));
end
clear M1_replace M2_replace A_replace B_replace A_replace_centered B_replace_centered;

%% 计算样本相关系数幅度的统计量，用于窗口标定
rho_hat_abs_mean = mean(abs_rho_MC, 1); 
rho_hat_abs_std = std(abs_rho_MC, 0, 1);
rho_hat_abs_mean_replace = mean(abs_rho_MC_replace, 1);
rho_hat_abs_std_replace = std(abs_rho_MC_replace, 0, 1);

rho_matrix = repmat(rho, number_of_MC_trials, 1);
rho_bias = rho_hat_abs_mean - rho;
rho_rmse = sqrt(mean((abs(rho_MC) - rho_matrix).^2, 1));
rho_bias_replace = rho_hat_abs_mean_replace - rho;
rho_rmse_replace = sqrt(mean((abs(rho_MC_replace) - rho_matrix).^2, 1));


%% 绘制平均估计值与真值对比（不含标准差范围带）
figure;
hold on;
plot(rho, rho, '--k', 'LineWidth', 1.2, 'DisplayName', '$|\gamma|$');
plot(rho, rho_hat_abs_mean, 'b', 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{homo}$');
plot(rho, rho_hat_abs_mean_replace, 'r', 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{hetero}$');
xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('$|\hat{\gamma}|$', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'southeast');

title(sprintf('Homo vs Hetero for Window Size: %d, Contamination Ratio: %.2f, Amplification Scale: %.2f', window_size, contamination_ratio, amplification_scale));

grid on;
xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]); % 保证横纵坐标比例一致，不受窗口拉伸影响
hold off;

%%  绘制标准差对比
figure;
hold on;
plot(rho, rho_hat_abs_std, 'b', 'LineWidth', 2, 'DisplayName', '$\sigma_{homo}$');
plot(rho, rho_hat_abs_std_replace, 'r', 'LineWidth', 2, 'DisplayName', '$\sigma_{hetero}$');

xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('Standard Deviation', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northeast');

title(sprintf('Std Dev Comparison: Homo vs Hetero (Win: %d, Ratio: %.2f, Amp: %.2f)', window_size, contamination_ratio, amplification_scale));

grid on;
xlim([0 1]);
pbaspect([1 1 1]);
hold off;


%%  绘制平均估计值与真值对比（含标准差范围带）
figure;
hold on;
plot(rho, rho, '--k', 'LineWidth', 1.2, 'DisplayName', '$|\gamma|$');
plot(rho, rho_hat_abs_mean, 'b', 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{homo}$');
plot(rho, rho_hat_abs_mean_replace, 'r', 'LineWidth', 2, 'DisplayName', '$|\hat{\gamma}|_{hetero}$');

% 可视化一个标准差范围（原始与替换）
fill([rho, fliplr(rho)], [rho_hat_abs_mean - rho_hat_abs_std, fliplr(rho_hat_abs_mean + rho_hat_abs_std)], ...
    [0.8 0.85 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', '$\pm 1 \sigma_{homo}$');
fill([rho, fliplr(rho)], [rho_hat_abs_mean_replace - rho_hat_abs_std_replace, fliplr(rho_hat_abs_mean_replace + rho_hat_abs_std_replace)], ...
    [1 0.85 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', '$\pm 1 \sigma_{hetero}$');

xlabel('$|\gamma|$', 'Interpreter', 'latex');
ylabel('$|\hat{\gamma}|$', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'southeast');

title(sprintf('Homo vs Hetero for Window Size: %d, Contamination Ratio: %.2f, Amplification Scale: %.2f', window_size, contamination_ratio, amplification_scale));

grid on;
xlim([0 1]);
ylim([0 1]);
pbaspect([1 1 1]); % 保证横纵坐标比例一致，不受窗口拉伸影响
hold off;

