function amp_values = find_pixels_by_cv(img, target_cv, cv_tolerance, win_size)
% find_pixels_by_cv 依据CV(变异系数)筛选像元，并扩展为窗口覆盖区域后返回像元值。
%
% 用途：
%   1) 先对输入图像计算局部CV图(CVEst)；
%   2) 找到CV落在 [target_cv-cv_tolerance, target_cv+cv_tolerance] 的像元位置；
%   3) 以这些像元为中心，取 win_size×win_size 的邻域窗口并做并集(重叠去重)；
%   4) 返回该并集区域内的所有像元值向量。
%
% 输入：
%   img          : 待筛选的图像矩阵(通常为幅度/功率/滤波后的sigma等)。
%   target_cv    : 目标CV值。
%   cv_tolerance : CV容差，实际筛选范围为 target_cv±cv_tolerance。
%   win_size     : 以候选像元为中心扩展的窗口大小(标量，建议奇数)。
%
% 输出：
%   amp_values   : 被选中区域(所有窗口的并集)内的像元值(列向量形式)。
%
% 备注：
%   - 该函数只做“筛选+取值”，不改变像元值本身。
%   - 若 img 中含 NaN/Inf，CVEst 的卷积统计可能会传播 NaN；建议在上游先处理无效值。

% 目标CV范围
cv_min = target_cv - cv_tolerance;
cv_max = target_cv + cv_tolerance;

% 扩展窗口的半径(例如 win_size=21 => half_win=10)
half_win = floor(win_size / 2);

% 计算局部CV图：每个像元用 win_size×win_size 邻域估计CV
CV_img = CVEst(img, [win_size win_size]);
[nrows, ncols] = size(img);

% 在CV图上做阈值筛选，得到候选像元(中心点)
cv_mask = (CV_img >= cv_min) & (CV_img <= cv_max);
[row_idx, col_idx] = find(cv_mask);
num_points = length(row_idx);
fprintf('找到 %d 个CV在 [%.2f, %.2f] 范围内的点\n', num_points, cv_min, cv_max);

% selected_mask 用于累积所有中心点对应的窗口并集；逻辑矩阵自动完成“重叠去重”
selected_mask = false(nrows, ncols);

% 对每个候选中心点，取其邻域窗口并加入并集掩膜
for i = 1:num_points
    r = row_idx(i);
    c = col_idx(i);
    
    % 边界保护：窗口落在图像范围内
    r_start = max(1, r - half_win);
    r_end = min(nrows, r + half_win);
    c_start = max(1, c - half_win);
    c_end = min(ncols, c + half_win);
    
    selected_mask(r_start:r_end, c_start:c_end) = true;
end
% 按并集掩膜取值；返回为向量(按列主序展开)
amp_values = img(selected_mask);
fprintf('共选中 %d 个像素（重叠区域已去重）\n', length(amp_values));
end
