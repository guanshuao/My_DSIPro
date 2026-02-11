function img_filtered = LeeFilter(img, win_size, looks)

    % 强制输入三个参数
    if nargin < 3
        error('需要三个输入参数：img, win_size, looks');
    end

    % 在计算局部均值/方差时忽略无效像元，输出中保持这些像元为 NaN
    valid_mask = isfinite(img);
    img0 = img;
    img0(~valid_mask) = 0;
    % 定义均值滤波器核
    h = ones(win_size, win_size) / (win_size * win_size);

    % 计算局部均值/二阶矩
    valid_frac = imfilter(double(valid_mask), h, 'symmetric');
    denom = max(valid_frac, eps);
    I_mean = imfilter(img0, h, 'symmetric') ./ denom;
    I_sq_mean = imfilter(img0.^2, h, 'symmetric') ./ denom;
    % 若窗口内没有有效像元，则局部统计量置 NaN
    I_mean(valid_frac == 0) = NaN;
    I_sq_mean(valid_frac == 0) = NaN;

    % 计算局部方差
    I_var = I_sq_mean - I_mean.^2;

    % 计算VI（CV的平方，归一化方差）
    epsilon = 1e-10;
    VI = I_var ./ (I_mean.^2 + epsilon);

    % 计算Lee滤波权重
    k = 1 - (1 / looks) ./ (VI + epsilon);

    % 权重限制在 [0, 1] 之间
    k(k < 0) = 0; 
    k(k > 1) = 1;
    
    img_filtered = I_mean + k .* (img0 - I_mean);
    img_filtered(~valid_mask) = NaN;
end