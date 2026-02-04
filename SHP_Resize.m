function [SHP_small] = SHP_Resize(SHP_large, CalWin_small)
% SHP_Resize - 从大窗口的SHP结构体中提取小窗口的SHP结构体
%
%   Usage:
%       [SHP_small] = SHP_Resize(SHP_large, CalWin_small)
%
%   Inputs:
%       SHP_large    - 大窗口的SHP结构体（由SHP_SelPoint生成）
%       CalWin_small - 小窗口的大小，如[11, 11]，必须为奇数且不大于大窗口
%
%   Outputs:
%       SHP_small - 小窗口的SHP结构体，包含：
%           .PixelInd - 小窗口内的同质像元索引
%           .BroNum   - 每个像素的同质像元数量
%           .CalWin   - 小窗口大小
%
%   原理说明：
%       由于大窗口和小窗口的中心重合，且小窗口是大窗口的子集：
%       1. 从大窗口的PixelInd中提取小窗口对应位置的像素
%       2. 重新进行连通域标记（bwlabel），因为裁剪后连通性可能改变
%       3. 只保留与中心像素连通的同质像元
%
%   Example:
%       % 假设已有13x13窗口的SHP结构体
%       SHP_11 = SHP_Resize(SHP_13, [11, 11]);
%
%   See also: SHP_SelPoint
%
%   Mi JIANG, Sun Yat-sen University
%   ======================================================================

if nargin < 2
    help SHP_Resize
    return;
end

tic;

% 获取大窗口参数
CalWin_large = SHP_large.CalWin;

% 输入检查
if mod(CalWin_small(1), 2) == 0 || mod(CalWin_small(2), 2) == 0
    error('小窗口尺寸必须为奇数');
end

if CalWin_small(1) > CalWin_large(1) || CalWin_small(2) > CalWin_large(2)
    error('小窗口尺寸必须小于等于大窗口尺寸');
end

if CalWin_small(1) == CalWin_large(1) && CalWin_small(2) == CalWin_large(2)
    warning('小窗口与大窗口尺寸相同，直接返回原结构体');
    SHP_small = SHP_large;
    return;
end

% 获取原始图像尺寸
[nlines, nwidths] = size(SHP_large.BroNum);
nPixels = nlines * nwidths;

% 计算大窗口的中心位置（1-indexed）
CenterRow_large = (CalWin_large(1) + 1) / 2;
CenterCol_large = (CalWin_large(2) + 1) / 2;

% 计算小窗口的中心位置（1-indexed）
CenterRow_small = (CalWin_small(1) + 1) / 2;
CenterCol_small = (CalWin_small(2) + 1) / 2;

% 计算小窗口的半径
RadiusRow_small = (CalWin_small(1) - 1) / 2;
RadiusCol_small = (CalWin_small(2) - 1) / 2;

% 小窗口在大窗口中的起始和结束行列（1-indexed）
StartRow = CenterRow_large - RadiusRow_small;
EndRow = CenterRow_large + RadiusRow_small;
StartCol = CenterCol_large - RadiusCol_small;
EndCol = CenterCol_large + RadiusCol_small;

% 计算小窗口内像素在大窗口PixelInd中对应的行索引
% 大窗口中，PixelInd的行按照列优先（MATLAB默认）排列
% 即大窗口中位置(row, col)的索引为：(col-1)*CalWin_large(1) + row
idx_small = zeros(CalWin_small(1) * CalWin_small(2), 1);
count = 0;
for col = StartCol:EndCol
    for row = StartRow:EndRow
        count = count + 1;
        idx_large = (col - 1) * CalWin_large(1) + row;
        idx_small(count) = idx_large;
    end
end

% 初始化小窗口的PixelInd
SHP_small.PixelInd = false(CalWin_small(1) * CalWin_small(2), nPixels);

% 显示进度
disp('正在重新计算连通性...');
progress_step = floor(nPixels / 10);

% 对每个像素重新计算连通性
for p = 1:nPixels
    % 从大窗口提取小窗口区域的同质像元掩膜
    mask_large = SHP_large.PixelInd(:, p);
    mask_small_vec = mask_large(idx_small);
    
    % 重塑为小窗口形状
    mask_small_2d = reshape(mask_small_vec, [CalWin_small(1), CalWin_small(2)]);
    
    % 重新进行连通域标记
    LL = bwlabel(mask_small_2d);
    
    % 只保留与中心像素连通的区域
    center_label = LL(CenterRow_small, CenterCol_small);
    if center_label > 0
        connected_mask = (LL == center_label);
    else
        % 如果中心像素不在任何连通域中，只保留中心像素自身
        connected_mask = false(CalWin_small);
        connected_mask(CenterRow_small, CenterCol_small) = true;
    end
    
    % 存储结果（按列优先展开）
    SHP_small.PixelInd(:, p) = connected_mask(:);
    
    % 显示进度
    if mod(p, progress_step) == 0
        disp(['进度: ', num2str(round(100 * p / nPixels)), '%']);
    end
end

% 计算BroNum
BroNum_vec = sum(SHP_small.PixelInd, 1);
SHP_small.BroNum = uint16(reshape(BroNum_vec(:), [nlines, nwidths]));

% 记录小窗口大小
SHP_small.CalWin = CalWin_small;

t = toc;

% 可视化对比
figure;
subplot(1, 2, 1);
imagesc(SHP_large.BroNum); axis image off; colormap jet;
title(['大窗口 (', num2str(CalWin_large(1)), 'x', num2str(CalWin_large(2)), ') SHP数量']);
colorbar;

subplot(1, 2, 2);
imagesc(SHP_small.BroNum); axis image off; colormap jet;
title(['小窗口 (', num2str(CalWin_small(1)), 'x', num2str(CalWin_small(2)), ') SHP数量']);
colorbar;

disp(['SHP_Resize: ', num2str(CalWin_large(1)), 'x', num2str(CalWin_large(2)), ...
      ' -> ', num2str(CalWin_small(1)), 'x', num2str(CalWin_small(2))]);
disp(['操作完成，耗时 ', num2str(t, '%.2f'), ' 秒']);
disp('Done!');

end
