clc;clear;close all;
img = SingleRead('E:\Coh_Est\SL\MLI\20220518.rslc',1201,'cpxfloat32');
img = flip(img);

amp = abs(img);
power = abs(img).^2;
sigma = LeeFilter(power);

CV_amp = CVEst(amp, [21 21]);
CV_power = CVEst(power, [21 21]);
CV_sigma = CVEst(sigma, [21 21]);


%% 展示滤波前后图像对比
figure;
ax1 = subplot(3,1,1);
mlishow(amp);
title('Original Amplitude Image');
axis image;
ax2 = subplot(3,1,2);
mlishow(power);
title('Original Power Image');
axis image;
ax3 = subplot(3,1,3);
mlishow(sigma);
title('Sigma Image');
axis image;
linkaxes([ax1, ax2, ax3], 'xy');


%% 展示滤波前后直方图对比
% 只统计中间部分像素值的直方图，防止拖尾太长
figure;
subplot(3,1,1);
histogram(amp(:), 'Normalization', 'pdf', 'BinLimits', prctile(amp(:), [0, 97]));
title('Original Histogram of Amplitude');
subplot(3,1,2);
histogram(power(:), 'Normalization', 'pdf', 'BinLimits', prctile(power(:), [0, 97]));
title('Original Histogram of Power');
subplot(3,1,3);
histogram(sigma(:), 'Normalization', 'pdf', 'BinLimits', prctile(sigma(:), [0, 97]));
title('Histogram of Sigma');


%% 展示滤波前后CV（变异系数）图像对比
figure;
ax1 = subplot(3,1,1);
imagesc(CV_amp);
colorbar;
title('CV of Original Amplitude');
axis image;
ax2 = subplot(3,1,2);
imagesc(CV_power);
colorbar;
title('CV of Original Power');
axis image;
ax3 = subplot(3,1,3);
imagesc(CV_sigma);
colorbar;
title('CV of Sigma');
axis image;
linkaxes([ax1, ax2, ax3], 'xy');


%% 展示滤波前后CV（变异系数）直方图对比
figure;
ax1 = subplot(3,1,1);
histogram(CV_amp(:), 'Normalization', 'pdf', 'BinLimits', [0 prctile(CV_amp(:), 99)]);
title('CV Histogram of Original Amplitude');
ax2 = subplot(3,1,2);
histogram(CV_power(:), 'Normalization', 'pdf', 'BinLimits', [0 prctile(CV_power(:), 99)]);
title('CV Histogram of Original Power');
ax3 = subplot(3,1,3);
histogram(CV_sigma(:), 'Normalization', 'pdf', 'BinLimits', [0 prctile(CV_sigma(:), 99)]);
title('CV Histogram of Sigma');