function CV_map = CVEst(img, calwin)
    % CVEst 计算图像的变异系数(CV)地图。
    %   CV_map = CVEst(img, calwin)
    %   img    : 图像数据矩阵（可为灰度或RGB），支持包含NaN的数据。
    %   calwin : 计算窗口大小，[height width]，如 [21 21]。若为标量则使用方窗。
    %
    % 计算步骤：局部均值 mu、局部平方均值 E[X^2]，得到方差 sigma^2 = E[X^2]-mu^2，
    % 再计算 CV = sigma / mu，并对均值为 0 的位置置 0，边界通过零填充处理。
    % NaN值会被自动忽略，仅使用窗口内有效像素进行计算。

    if isscalar(calwin)
        calwin = [calwin calwin];
    end

    if numel(calwin) ~= 2 || any(calwin <= 0) || any(mod(calwin,2) == 0)
        error('calwin 必须是正奇数的 [height width] 或标量。');
    end

    if ndims(img) == 3
        img = rgb2gray(img);
    end

    img = double(img);

    % ------- 处理NaN：创建有效像素掩码 -------
    valid_mask = ~isnan(img);      % 有效像素掩码（非NaN为1，NaN为0）
    img_clean = img;
    img_clean(~valid_mask) = 0;    % 将NaN替换为0，便于卷积计算

    % ------- 核与有效像元计数 -------
    kernel = ones(calwin); % 窗口全 1 卷积核
    % 只计算有效像素的数量（NaN位置不计入）
    N_count = conv2(double(valid_mask), kernel, 'same');

    % ------- 局部统计量 -------
    sum_img = conv2(img_clean, kernel, 'same');
    mean_img = sum_img ./ N_count;

    sum_img2 = conv2(img_clean.^2, kernel, 'same');
    mean_img2 = sum_img2 ./ N_count;

    % 方差与标准差
    var_img = mean_img2 - mean_img.^2;
    var_img(var_img < 0) = 0; % 修正数值误差
    std_img = sqrt(var_img);

    % ------- 变异系数 -------
    CV_map = std_img ./ mean_img;
    CV_map(mean_img == 0) = 0;      % 避免除以 0
    CV_map(N_count == 0) = NaN;     % 窗口内全为NaN时，结果为NaN
end
