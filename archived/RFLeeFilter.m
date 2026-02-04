function img_filtered = RFLeeFilter(img, win_size, ENL)
%RFLeeFilter Refined Lee speckle filter (Refined Filtering Using Local Statistics)
%
%   精致Lee滤波器的优化实现。
%   该函数实现了 J.-S. Lee (1980) 描述的算法。
%   它利用局部统计量（均值和方差）和边缘检测来选择均匀的邻域进行滤波。
%
%   输入:
%       img      - 输入图像（强度/功率）。NaN值将被忽略（不参与统计计算）。
%       win_size - 窗口大小。必须为 7（由算法定义固定）。
%                  如果未提供，默认为 7。
%       ENL      - 等效视数 (Equivalent Number of Looks)。用于在局部估计不稳定时
%                  估计噪声方差。默认为 1。
%
%   输出:
%       img_filtered - 滤波后的图像。
%
%   Reference:
%       J.-S. Lee, “Refined Filtering of Image Noise Using Local Statistics”, 
%       IEEE Trans. Pattern Anal. Mach. Intell., vol. PAMI-3, no. 2, 1981.

    if nargin < 2 || isempty(win_size)
        win_size = 7;
    end
    if win_size ~= 7
        error('RFLeeFilter:winSize', 'Refined Lee algorithm requires 7x7 window.');
    end
    if nargin < 3 || isempty(ENL)
        ENL = 1;
    end

    % 1. 预处理：处理 NaN 和填充
    % 记录 NaN 的位置，并在计算统计量时将其排除
    nan_mask = isnan(img);
    img(nan_mask) = 0; % 将 NaN 替换为 0 以便进行卷积求和，但在计数时会排除
    valid_mask = ~nan_mask;
    
    [h, w] = size(img);
    
    % 填充图像以处理边界。'symmetric' 填充可减少边缘伪影。
    pad_width = 3;
    img_pad = padarray(img, [pad_width, pad_width], 'symmetric');
    valid_mask_pad = padarray(valid_mask, [pad_width, pad_width], 'symmetric');
    
    % 2. 计算整个填充图像的 3x3 均值和方差图
    % 这些将用于构建 9 个子区域的统计数据。
    kernel3 = ones(3);
    
    % 计算每个 3x3 窗口内的有效像素数
    count3 = conv2(double(valid_mask_pad), kernel3, 'valid');
    count3(count3 == 0) = 1; % 避免除以零
    
    % 'valid' 卷积使尺寸减小 2 (3x3 核)。
    % 结果对应于以 img_pad 的 (r+1, c+1) 为中心的 3x3 块均值。
    sum3 = conv2(img_pad, kernel3, 'valid'); 
    mu3 = sum3 ./ count3;
    
    sum3_sq = conv2(img_pad.^2, kernel3, 'valid');
    var3 = (sum3_sq - (sum3.^2)./count3) ./ (count3 - 1);
    var3(count3 <= 1) = 0;
    var3(var3 < 0) = 0; % 确保方差非负
    
    % 3. 提取每个像素的 9 个子区域统计量
    % 我们提取 mu3/var3 的 9 个移位版本，以与 7x7 窗口结构对齐。
    % 对于像素 (r,c)，7x7 窗口有 9 个子块。
    % 这些子块在填充图像坐标系中的中心允许我们正确地切片 mu3。
    
    M = cell(3,3);
    V = cell(3,3);
    
    for r_idx = 1:3
        for c_idx = 1:3
            % 根据 7x7 窗口内的 3x3 网格计算偏移量
            r_offset = 2*(r_idx-1);
            c_offset = 2*(c_idx-1);
            % 提取相应的切片
            M{r_idx, c_idx} = mu3(1+r_offset : h+r_offset, 1+c_offset : w+c_offset);
            V{r_idx, c_idx} = var3(1+r_offset : h+r_offset, 1+c_offset : w+c_offset);
        end
    end
    
    % 4. 计算梯度和方向
    % 计算 4 个主要方向的梯度
    G1 = abs(M{2,1} - M{2,3}); % 水平
    G2 = abs(M{1,3} - M{3,1}); % 反对角线
    G3 = abs(M{1,2} - M{3,2}); % 垂直
    G4 = abs(M{1,1} - M{3,3}); % 主对角线
    
    % 堆叠梯度以找到最大值
    G_stack = cat(3, G1, G2, G3, G4);
    [~, coarse_dir] = max(G_stack, [], 3); % 1..4
    
    % 将方向细化为 8 个方向 (0..7)
    D = zeros(h, w);
    
    % 细化条件
    C1 = abs(M{2,1} - M{2,2}) < abs(M{2,2} - M{2,3});
    C2 = abs(M{1,3} - M{2,2}) < abs(M{2,2} - M{3,1});
    C3 = abs(M{1,2} - M{2,2}) < abs(M{2,2} - M{3,2});
    C4 = abs(M{1,1} - M{2,2}) < abs(M{2,2} - M{3,3});
    
    % 应用细化逻辑
    D(coarse_dir==1 & ~C1) = 0; % 右
    D(coarse_dir==1 & C1)  = 4; % 左
    D(coarse_dir==2 & C2)  = 1; % 右上
    D(coarse_dir==2 & ~C2) = 5; % 左下
    D(coarse_dir==3 & C3)  = 2; % 上
    D(coarse_dir==3 & ~C3) = 6; % 下
    D(coarse_dir==4 & C4)  = 3; % 左上
    D(coarse_dir==4 & ~C4) = 7; % 右下
    
    % 5. 计算选定窗口（28 像素）的局部统计量
    % 与 8 个特定掩模进行卷积，以获得每个潜在方向的均值和方差
    masks = getRefinedLeeMasks();
    
    means_stack = zeros(h, w, 8);
    vars_stack = zeros(h, w, 8);
    
    img_sq_pad = img_pad.^2;
    
    for k = 1:8
        mask = masks{k};
        
        % 计算当前掩模下的有效像素数
        count_k = conv2(double(valid_mask_pad), mask, 'valid');
        count_k(count_k == 0) = 1; % 避免除以零
        
        s = conv2(img_pad, mask, 'valid');
        ss = conv2(img_sq_pad, mask, 'valid');
        
        m = s ./ count_k;
        v = (ss - (s.^2)./count_k) ./ (count_k - 1);
        v(count_k <= 1) = 0;
        v(v < 0) = 0;
        
        means_stack(:,:,k) = m;
        vars_stack(:,:,k) = v;
    end
    
    % 选择与确定方向 D 对应的统计量
    % 将 D (0..7) 转换为线性索引 (1..8)
    indices = sub2ind([h, w, 8], repmat((1:h)', 1, w), repmat(1:w, h, 1), D + 1);
    
    meanY = means_stack(indices);
    varY = vars_stack(indices);
    
    % 6. 估计局部噪声方差 (SigmaV)
    % 计算 9 个子区域的归一化方差
    sigmaV_prior = 1 / ENL;
    
    NV_stack = zeros(h, w, 9);
    for i = 1:9
        m_sub = M{i};
        v_sub = V{i};
        
        nv_sub = v_sub ./ (m_sub.^2);
        % 将无效计算（例如 mean=0）标记为 NaN
        nv_sub(m_sub <= 0) = NaN; 
        
        NV_stack(:,:,i) = nv_sub;
    end
    
    % 排序并取最小的 5 个的平均值
    NV_sorted = sort(NV_stack, 3);
    NV_top5 = NV_sorted(:,:,1:5);
    sigmaV_est = mean(NV_top5, 3, 'omitnan');
    
    % 如果估计失败（全为 NaN），则使用先验值
    sigmaV_est(isnan(sigmaV_est)) = sigmaV_prior;
    
    % 7. 应用 Lee 滤波器公式
    % 权重 W = (varY - meanY^2 * sigmaV) / (varY * (1 + sigmaV))
    
    center = img;
    
    num = varY - (meanY.^2) .* sigmaV_est;
    den = varY .* (1 + sigmaV_est);
    
    % 避免除以零
    W = zeros(h, w);
    valid_den = den > 0;
    W(valid_den) = num(valid_den) ./ den(valid_den);
    
    % 将 W 裁剪到 [0, 1]（因为方差不能为负）
    W(W < 0) = 0;
    W(W > 1) = 1;
    
    img_filtered = meanY + W .* (center - meanY);
    
    % 恢复原始的 NaN 值
    img_filtered(nan_mask) = NaN;

end

function masks = getRefinedLeeMasks()
    % 定义 7x7 窗口的 8 个掩模
    masks = cell(1,8);
    base = zeros(7);
    [c, r] = meshgrid(1:7, 1:7);
    
    % d=0: 右 (列 4-7)
    m = base; m(:, 4:7) = 1; masks{1} = m;
    
    % d=1: 右上 (c >= r)
    m = base; m(c >= r) = 1; masks{2} = m;
    
    % d=2: 上 (行 1-4)
    m = base; m(1:4, :) = 1; masks{3} = m;
    
    % d=3: 左上 (r+c <= 8)
    m = base; m(r+c <= 8) = 1; masks{4} = m;
    
    % d=4: 左 (列 1-4)
    m = base; m(:, 1:4) = 1; masks{5} = m;
    
    % d=5: 左下 (c <= r)
    m = base; m(c <= r) = 1; masks{6} = m;
    
    % d=6: 下 (行 4-7)
    m = base; m(4:7, :) = 1; masks{7} = m;
    
    % d=7: 右下 (r+c >= 8)
    m = base; m(r+c >= 8) = 1; masks{8} = m;
end