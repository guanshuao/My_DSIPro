clear;clc;close all;

temp_loader = load('/titan/guanshuao/Kumamoto/DSIPro/AdpCoh/average.mat');
temp_fields = fieldnames(temp_loader);
adpcoh = temp_loader.(temp_fields{1});
clear temp_loader temp_fields;
disp('Adaptive coherence loaded.');

temp_loader = load('/titan/guanshuao/Kumamoto/DSIPro/BoxCoh/average.mat');
temp_fields = fieldnames(temp_loader);
boxcoh = temp_loader.(temp_fields{1});
clear temp_loader temp_fields;
disp('Boxcar coherence loaded.');

load('SHP_FAST.mat');
disp('SHP loaded.');

load('LUT.mat');
disp('LUT loaded.');

temp_loader = load('avg_amp.mat');
temp_fields = fieldnames(temp_loader);
mla_image = temp_loader.(temp_fields{1});
clear temp_loader temp_fields;
disp('Average amplitude image loaded.');

temp_loader = load('avg_intensity.mat');
temp_fields = fieldnames(temp_loader);
mli_image = temp_loader.(temp_fields{1});
clear temp_loader temp_fields;
disp('Average intensity image loaded.');


%%%%%%%%%%%%%%%%
CV_mli_map = CVEst(mli_image,21);
CV_mla_map = CVEst(mla_image,21);



% 首先进行同质像元选取，开一个21x21的窗口，共441个像元，在该窗口内选取与之同类的地物，记录其数量即为SHP Number(Statistically Homogeneous Pixels Number)
% Boxcar方法估计相关系数，计算的是两个窗口内所有像元的相关系数
% 自适应方法估计相关系数，只会选取窗口内的同质像元参与计算
% Modified Adaptive方法，通过查找表，把使用SHPNumber个像元计算得到的相关系数，转换为使用441个像元计算得到的相关系数，从而消除由于像元数量变化带来的误差


SHP_num = SHP_FAST.BroNum; % SHP数量图，例如某个像元的值为X，则证明该像元周围有X个像元，是它的同类像元
% boxcoh = mean(BoxCoh_ml,3,'omitnan'); % 用boxcar方法估计的相干系数图
% adpcoh = mean(UbCoh_ml,3,'omitnan'); % 用自适应方法估计的相干系数图

ratio = zeros(3,441); % 1: Box/Adp, 2: Box/Mod, 3: Adp/Mod
% 1是boxcar与adaptive的比值，2反映了不均一带来的误差，3反映了像素数变化带来的误差
one = ones(1,441); % 用以绘制y=1的参考线
coh = zeros(4,441); % 用以存储boxcar、自适应以及修正后的相关系数、真实相关系数的均值
stds = zeros(3,441); % 用以存储boxcar、自适应以及修正后的相关系数的标准差
CV = zeros(2,441); % 第一行存储振幅图CV，第二行存储强度图CV

for n_shp = 1:441
    pixel_indices = find(SHP_num == n_shp); % 找出SHP Number等于n_shp的像元索引
    CV_mla_vals = CV_mla_map(pixel_indices); % SHP Number为n_shp的像元对应的CV值
    CV_mli_vals = CV_mli_map(pixel_indices); % SHP Number为n_shp的像元对应的CV值
    CV(1,n_shp) = mean(CV_mla_vals(:)); % 计算振幅图CV的均值
    CV(2,n_shp) = mean(CV_mli_vals(:)); % 计算强度图CV的均值
end

figure('Visible', 'off');
plot(1:441,CV(1,:),'r', 'LineWidth', 3);
hold on;
plot(1:441,CV(2,:),'b', 'LineWidth', 3);
legend('Amplitude CV','Intensity CV');
xlim([1,441]);
xlabel('SHP Number');
ylabel('Coefficient of Variation (CV)');
title('Mean Coefficient of Variation (CV)');
hold off;
saveas(gcf, 'Mean_CV.png');

for n_shp = 1:441
    pixel_indices = find(SHP_num == n_shp); % 找出SHP Number等于n_shp的像元索引
    box_coh_vals = boxcoh(pixel_indices); % SHP Number为n_shp的像元对应的boxcar相关系数
    adp_coh_vals = adpcoh(pixel_indices); % SHP Number为n_shp的像元对应的自适应相关系数
    coh(1,n_shp) = mean(box_coh_vals(:),'omitnan'); % 计算boxcar相关系数的均值
    coh(2,n_shp) = mean(adp_coh_vals(:),'omitnan'); % 计算自适应相关系数的均值
    stds(1,n_shp) = std(box_coh_vals(:),'omitnan'); % 计算boxcar相关系数的标准差
    stds(2,n_shp) = std(adp_coh_vals(:),'omitnan'); % 计算自适应相关系数的标准差

    % 自适应修正计算 (Adaptive Modified Calculation)
    idx_n = find(LUT.n_values == n_shp); % 在查找表中找到当前SHP数量对应的索引
    if ~isempty(idx_n)
        % 反向查找：给定估计的相关系数 rho_hat (即 coh(2,n_shp)) 和样本数 n (即 n_shp)，查找真实的 rho
        expected_curve = LUT.mean_rho_hat(idx_n, :); % 获取该样本数下的期望相关系数曲线
        
        % 过滤掉期望曲线中的非有限值 (NaN 或 Inf)
        valid_mask = isfinite(expected_curve);
        expected_curve_clean = expected_curve(valid_mask);
        rho_values_clean = LUT.rho_values(valid_mask);

        % 确保插值点唯一，去除重复值
        [unique_curve, unique_idx] = unique(expected_curve_clean);
        
        if length(unique_curve) > 1 && isfinite(coh(2,n_shp))
             % 通过插值反推真实的相干系数 true_rho
             true_rho = interp1(unique_curve, rho_values_clean(unique_idx), coh(2,n_shp), 'linear', 'extrap');
             true_rho = max(0, min(1, true_rho)); % 将结果限制在 [0, 1] 范围内
             
             % 正向查找：给定真实的 rho 和 n=441 (假设满窗口)，查找对应的估计相关系数 rho_hat
             idx_441 = find(LUT.n_values == 441); % 找到 n=441 在查找表中的索引
             if ~isempty(idx_441)
                 expected_curve_441 = LUT.mean_rho_hat(idx_441, :); % 获取 n=441 时的期望曲线
                 coh(3,n_shp) = interp1(LUT.rho_values, expected_curve_441, true_rho, 'linear'); % 插值得到修正后的相关系数
                 coh(4,n_shp) = true_rho; % 记录真实的相干系数

                 % 计算修正后的标准差 (Calculate Modified Std)
                 std_curve_441 = LUT.std_rho_hat(idx_441, :);
                 stds(3,n_shp) = interp1(LUT.rho_values, std_curve_441, true_rho, 'linear');
             else
                 coh(3,n_shp) = NaN;
                 coh(4,n_shp) = NaN;
                 stds(3,n_shp) = NaN;
             end
        else
             coh(3,n_shp) = NaN;
             coh(4,n_shp) = NaN;
             stds(3,n_shp) = NaN;
        end
    else
        coh(3,n_shp) = NaN;
        coh(4,n_shp) = NaN;
        stds(3,n_shp) = NaN;
    end

    % 计算比率 (Calculate ratios)
    ratio(1,n_shp) = coh(1,n_shp)/coh(2,n_shp); % Boxcar / Adaptive (Boxcar与自适应的比值)
    if ~isnan(coh(3,n_shp)) && coh(3,n_shp) ~= 0
        ratio(2,n_shp) = coh(1,n_shp)/coh(3,n_shp); % Boxcar / Modified (Boxcar与修正后的比值)
        ratio(3,n_shp) = coh(2,n_shp)/coh(3,n_shp); % Adaptive / Modified (自适应与修正后的比值)
    else
        ratio(2,n_shp) = NaN;
        ratio(3,n_shp) = NaN;
    end
end


figure('Visible', 'off');
plot(1:441,coh(1,:),'r', 'LineWidth', 3);
hold on;
plot(1:441,coh(2,:),'b', 'LineWidth', 3);
plot(1:441,coh(3,:),'g', 'LineWidth', 3);
plot(1:441,coh(4,:),'k--', 'LineWidth', 2);
legend('Boxcar','Adaptive','Adaptive\_Modified','''True'' Coherence');
xlim([1,441]);
ylim([0 1]);
xlabel('SHP Number');
ylabel('$|\hat{\rho}|$', 'Interpreter', 'latex');
title('Mean coherence');
hold off;
saveas(gcf, 'Mean_Coherence.png');

figure('Visible', 'off');
plot(1:441,ratio(1,:),'r', 'LineWidth', 3);
hold on;
plot(1:441,ratio(2,:),'b', 'LineWidth', 3);
plot(1:441,ratio(3,:),'g', 'LineWidth', 3);
xlim([1,441]);
ylim([0 2]);
plot(1:441,one,'k:', 'LineWidth', 1);
legend('Boxcar/Adaptive','Boxcar/Modified','Adaptive/Modified');
xlabel('SHP Number');
ylabel('Ratio', 'Interpreter', 'latex');
title('Ratio of mean coherence');
hold off;
saveas(gcf, 'Ratio_Mean_Coherence.png');

figure('Visible', 'off');
subplot(1,2,1);
imagesc(boxcoh);
axis image;axis off;colorbar;
title('Coherence map estimated by the Boxcar method');
subplot(1,2,2);
imagesc(adpcoh);
axis image;axis off;colorbar;
title('Coherence map estimated by the adaptive method');
saveas(gcf, 'Coherence_Map_Comparison.png');

figure('Visible', 'off');
plot(1:441,stds(1,:),'r', 'LineWidth', 3);
hold on;
plot(1:441,stds(2,:),'b', 'LineWidth', 3);
plot(1:441,stds(3,:),'g', 'LineWidth', 3);
legend('Boxcar','Adaptive','Adaptive\_Modified');
xlim([1,441]);
xlabel('SHP Number');
ylabel('Standard Deviation');
title('Standard deviation of coherence');
hold off;
saveas(gcf, 'Std_Dev_Coherence.png');

figure('Visible', 'off');
x = 1:441;
hold on;
% Dummy patch for legend
h_dummy = fill([-10 -5 -5], [-10 -5 -10], [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Boxcar area
fill([x, fliplr(x)], [coh(1,:)-stds(1,:), fliplr(coh(1,:)+stds(1,:))], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Adaptive area
fill([x, fliplr(x)], [coh(2,:)-stds(2,:), fliplr(coh(2,:)+stds(2,:))], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot means
h1 = plot(x, coh(1,:), 'r', 'LineWidth', 3);
h2 = plot(x, coh(2,:), 'b', 'LineWidth', 3);

legend([h1, h2, h_dummy], {'Boxcar', 'Adaptive', '\pm\sigma'}, 'Location', 'Best');
xlim([1, 441]);
ylim([0 1]);
xlabel('SHP Number');
ylabel('Coherence');
title('Mean coherence with \pm\sigma range');
hold off;
saveas(gcf, 'Mean_Coherence_Sigma.png');