% gen_bias_data.m
% 蒙特卡洛生成用来绘制不均匀场景下的相干性估计图的数据
clear; clc;

rho = linspace(0, 1, 101)';
number_of_MC_trials = 10000; % 为了节省时间使用1万次，如果需要更高精度可以改为10万

%% 1. Varying w (r = 0.5, h = 3)
w_list = [3, 7];
contamination_ratio = 0.5;
amplification_scale = 3;
bias_w_data = rho;

fprintf('Generating Data for Figure 1 (varying w)...\n');
for w = w_list
    number_of_pixels = w^2;
    rho_hat_abs_mean = zeros(length(rho), 1);
    rho_hat_abs_mean_replace = zeros(length(rho), 1);
    
    parfor k = 1:length(rho)
        M1 = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
        M2 = rho(k) * M1 + sqrt(1 - rho(k)^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));
        
        M1_replace = M1; M2_replace = M2;
        % 针对每个trial分别随机选取点
        for t = 1:number_of_MC_trials
            idx = randperm(number_of_pixels, round(contamination_ratio * number_of_pixels));
            M1_replace(t, idx) = amplification_scale * M1_replace(t, idx);
            M2_replace(t, idx) = amplification_scale * M2_replace(t, idx);
        end
        
        % Homo
        A_c = M1 - mean(M1, 2); B_c = M2 - mean(M2, 2);
        num = sum(A_c .* conj(B_c), 2);
        den = sqrt(sum(abs(A_c).^2, 2) .* sum(abs(B_c).^2, 2));
        rho_hat_abs_mean(k) = mean(abs(num ./ max(den, eps)));
        
        % Hetero
        Ar_c = M1_replace - mean(M1_replace, 2); Br_c = M2_replace - mean(M2_replace, 2);
        num_r = sum(Ar_c .* conj(Br_c), 2);
        den_r = sqrt(sum(abs(Ar_c).^2, 2) .* sum(abs(Br_c).^2, 2));
        rho_hat_abs_mean_replace(k) = mean(abs(num_r ./ max(den_r, eps)));
    end
    bias_w_data = [bias_w_data, rho_hat_abs_mean, rho_hat_abs_mean_replace];
end
writematrix(bias_w_data, './bias_w_data.txt', 'Delimiter', '\t');


%% 2. Varying h (w = 3, r = 0.5)
h_list = [1, 2, 4];
w = 3; number_of_pixels = w^2;
contamination_ratio = 0.5;
bias_h_data = rho;

fprintf('Generating Data for Figure 2 (varying h)...\n');
for h = h_list
    rho_hat_abs_mean_replace = zeros(length(rho), 1);
    parfor k = 1:length(rho)
        M1 = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
        M2 = rho(k) * M1 + sqrt(1 - rho(k)^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));
        
        M1_replace = M1; M2_replace = M2;
        for t = 1:number_of_MC_trials
            idx = randperm(number_of_pixels, round(contamination_ratio * number_of_pixels));
            M1_replace(t, idx) = h * M1_replace(t, idx);
            M2_replace(t, idx) = h * M2_replace(t, idx);
        end
        
        Ar_c = M1_replace - mean(M1_replace, 2); Br_c = M2_replace - mean(M2_replace, 2);
        num_r = sum(Ar_c .* conj(Br_c), 2);
        den_r = sqrt(sum(abs(Ar_c).^2, 2) .* sum(abs(Br_c).^2, 2));
        rho_hat_abs_mean_replace(k) = mean(abs(num_r ./ max(den_r, eps)));
    end
    bias_h_data = [bias_h_data, rho_hat_abs_mean_replace];
end
writematrix(bias_h_data, './bias_h_data.txt', 'Delimiter', '\t');


%% 3. Varying r (w = 3, h = 3)
r_list = [0, 0.1, 0.3, 0.5];
w = 3; number_of_pixels = w^2;
amplification_scale = 3;
bias_r_data = rho;

fprintf('Generating Data for Figure 3 (varying r)...\n');
for r = r_list
    rho_hat_abs_mean_replace = zeros(length(rho), 1);
    parfor k = 1:length(rho)
        M1 = randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels);
        M2 = rho(k) * M1 + sqrt(1 - rho(k)^2) * (randn(number_of_MC_trials, number_of_pixels) + 1i * randn(number_of_MC_trials, number_of_pixels));
        
        M1_replace = M1; M2_replace = M2;
        for t = 1:number_of_MC_trials
            idx = randperm(number_of_pixels, round(r * number_of_pixels));
            M1_replace(t, idx) = amplification_scale * M1_replace(t, idx);
            M2_replace(t, idx) = amplification_scale * M2_replace(t, idx);
        end
        
        Ar_c = M1_replace - mean(M1_replace, 2); Br_c = M2_replace - mean(M2_replace, 2);
        num_r = sum(Ar_c .* conj(Br_c), 2);
        den_r = sqrt(sum(abs(Ar_c).^2, 2) .* sum(abs(Br_c).^2, 2));
        rho_hat_abs_mean_replace(k) = mean(abs(num_r ./ max(den_r, eps)));
    end
    bias_r_data = [bias_r_data, rho_hat_abs_mean_replace];
end
writematrix(bias_r_data, './bias_r_data.txt', 'Delimiter', '\t');
fprintf('Done.\n');
