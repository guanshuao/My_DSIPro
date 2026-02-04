function [Coh,Ph] = BoxCohEst(mlistack,mlilist,infstack,inflist,CalWin,OutputMode)
%   BoxCohEst 使用 Boxcar 滤波估计相干性和相位
%   支持多视 (N-looks) 图像。
%   用法:
%       [Coh,Ph] = BoxCohEst(mlistack,mlilist,infstack,inflist,CalWin,CohEstAgr)
%   或者:
%       [Coh] = BoxCohEst(mlistack,mlilist,infstack,inflist,CalWin,CohEstAgr)
%   输入:
%   - mlistack: 高度 x 宽度 x 页数 的 (实数) 矩阵，例如 SAR 单视强度序列
%   - mlilist:  n*1 的文件列表，包含以 <yyyymmdd> 命名的强度图像
%   - infstack: 高度 x 宽度 x 页数 的 (复数) 矩阵，例如去除相位梯度 (地形、形变) 后且幅度归一化的 InSAR 单视干涉图
%   - inflist:  n*2 的文件列表，包含以 <yyyymmdd yyyymmdd> 命名的复数干涉图
%   - CalWin:   用于收集样本进行相干性估计的窗口大小 [方位向 距离向]
%   - OutputMode: 输出模式，'stack' 或 'average'
%
%   输出:
%   - Coh:      相干性幅度堆栈 (实数)，如果只有一个输出参数，则返回相干性幅度
%   - Ph:       干涉相位 (实数)，如果有两个输出参数，则返回相位    
% 
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University,  
%   ======================================================================
%   11/2021 MJ 将高斯低通滤波器替换为均值滤波器
%   ======================================================================

if nargin < 6 || isempty(OutputMode) % 如果未提供 OutputMode 参数，则使用默认值 'stack'
   OutputMode = 'stack'; % 'stack' or 'average'
end

if nargin < 5 % 如果参数还有缺失，显示帮助信息并返回
    help BoxCohEst
    return
end

tic;
[nlines,nwidths,npages]=size(infstack); % 获取干涉图堆栈的尺寸
Coh=zeros(nlines,nwidths,npages,'single'); % 初始化相干性矩阵

[~,idx]=ismember(inflist,mlilist);
% 匹配干涉图列表和强度图列表，找到对应的主辅影像索引。
% 将干涉对列表 inflist 中的主/辅影像日期，在单视强度列表 mlilist 中查找其对应的位置索引：
%  - mlilist: n×1 字符串/元胞数组，每一行是一个强度图对应的获取日期（如 'yyyymmdd'）
%  - inflist: m×2 字符串/元胞数组，每一行是一个干涉对的主/辅影像日期（如 ['yyyymmdd' 'yyyymmdd']）
%  - ismember 对每一列分别在 mlilist 中查找，返回与 inflist 中每个日期匹配的行号
%  - idx 为 m×2 的整型矩阵，第 i 行 [p q] 表示第 i 个干涉对的主影像来自 mlistack(:,:,p)，辅影像来自 mlistack(:,:,q)

h = fspecial('average',CalWin);
% 生成一个二维均值滤波核（Boxcar 滤波器）：
%  - CalWin 为 [方位向 距离向]
%  - h 的尺寸为 CalWin(1)×CalWin(2)，元素值相等且和为 1，即一个均值滤波核
%  - 在后续使用 filter2(h, ...) 时，相当于在 CalWin 大小的窗口内对数据做局部平均
% 注意：多视数据只能使用 'average' (均值滤波)



if strcmpi(OutputMode,'stack')
    for ii=1:npages
        disp(['Processing coherence estimation for interferogram ',num2str(ii),' / ',num2str(npages),' ...']);
        m1 = mlistack(:,:,idx(ii,1));% 主影像强度
        m2 = mlistack(:,:,idx(ii,2));% 辅影像强度
        
        % 创建有效像素掩膜（非NaN位置）
        validMask = ~isnan(m1) & ~isnan(m2) & ~isnan(infstack(:,:,ii));
        
        % 将NaN替换为0进行滤波，同时记录有效样本数
        m1_clean = m1; m1_clean(~validMask) = 0;
        m2_clean = m2; m2_clean(~validMask) = 0;
        intf_clean = infstack(:,:,ii); intf_clean(~validMask) = 0;
        
        % 计算每个窗口内的有效样本数
        validCount = imfilter(double(validMask), h, 'symmetric') * numel(h);
        validCount(validCount == 0) = NaN;  % 避免除以0
        
        % 相干性估计器的分子项（使用sum代替mean，然后手动除以有效样本数）
        nu_sum = imfilter(sqrt(m1_clean.*m2_clean).*intf_clean, ones(CalWin), 'symmetric');
        nu = nu_sum ./ validCount;
        
        % 分母项
        de1_sum = imfilter(m1_clean, ones(CalWin), 'symmetric');
        de2_sum = imfilter(m2_clean, ones(CalWin), 'symmetric');
        de1 = de1_sum ./ validCount;
        de2 = de2_sum ./ validCount;
        
        % 计算复数相干性
        denominator = sqrt(de1.*de2);
        denominator(denominator == 0) = NaN;  % 避免除以0
        Coh(:,:,ii) = nu./denominator;
    end   

    % 处理NaN和Inf值
    Coh(isnan(Coh) | isinf(Coh)) = 0;
    if nargout > 1 % stack 模式下，如果请求两个输出，则计算相位
        Ph = angle(Coh);
    end
    Coh = abs(Coh); % 取相干性幅度
end

if strcmpi(OutputMode,'average')
    CohSum = zeros(nlines,nwidths,'single');
    for ii=1:npages
        disp(['Processing coherence estimation for interferogram ',num2str(ii),' / ',num2str(npages),' ...']);
        m1 = mlistack(:,:,idx(ii,1));% 主影像强度
        m2 = mlistack(:,:,idx(ii,2));% 辅影像强度
        
        % 创建有效像素掩膜（非NaN位置）
        validMask = ~isnan(m1) & ~isnan(m2) & ~isnan(infstack(:,:,ii));
        
        % 将NaN替换为0进行滤波，同时记录有效样本数
        m1_clean = m1; m1_clean(~validMask) = 0;
        m2_clean = m2; m2_clean(~validMask) = 0;
        intf_clean = infstack(:,:,ii); intf_clean(~validMask) = 0;
        
        % 计算每个窗口内的有效样本数
        validCount = imfilter(double(validMask), h, 'symmetric') * numel(h);
        validCount(validCount == 0) = NaN;  % 避免除以0
        
        % 相干性估计器的分子项（使用sum代替mean，然后手动除以有效样本数）
        nu_sum = imfilter(sqrt(m1_clean.*m2_clean).*intf_clean, ones(CalWin), 'symmetric');
        nu = nu_sum ./ validCount;
        
        % 分母项
        de1_sum = imfilter(m1_clean, ones(CalWin), 'symmetric');
        de2_sum = imfilter(m2_clean, ones(CalWin), 'symmetric');
        de1 = de1_sum ./ validCount;
        de2 = de2_sum ./ validCount;
        
        % 计算复数相干性
        denominator = sqrt(de1.*de2);
        denominator(denominator == 0) = NaN;  % 避免除以0
        CohSum = CohSum + abs(nu./denominator);
    end
    Coh = CohSum / npages;

    % 处理NaN和Inf值
    Coh(isnan(Coh) | isinf(Coh)) = 0;
    if nargout > 1
        warning('average 模式下只输出 Coh，Ph 输出将被忽略。');
    end
end

t=toc; % 计时结束
disp(['BoxCohEst operation completed in ',num2str(t/60),' min(s).']);
disp('Done!');    