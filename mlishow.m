function mlishow(mli,sc,exp)
%mlishow 用于显示强度图像
%   该函数对强度图像进行指数变换和缩放，以便更好地可视化。
%
%   输入:
%   - mli: 强度图像 (二维矩阵)
%   - sc:  display scale factor，显示缩放因子 (默认=0.35，注意：原注释写的是 1.，但代码中默认值是 0.35)
%   - exp: display exponent，显示指数 (默认=1，注意：原注释写的是 .35，但代码中默认值是 1)
%
%   Mi JIANG, Sun Yat-sen Univeristy,

if nargin < 3 || isempty(exp)
    exp=1; % 默认指数为 1 (线性显示)
end

if nargin < 2 || isempty(sc)
    sc=.35; % 默认缩放因子为 0.35
end

if nargin < 1
    help mlishow;
    return;
end

mli=abs(mli); % 取绝对值 (防止复数输入)
mli(isnan(mli))=0; % 将 NaN 替换为 0

% 计算非零像素的平均强度 (经过指数变换后)
P = mli.^exp;
nv = numel(nonzeros(P)); % 非零像素数量
P = sum(P(:)); % 所有像素强度之和

if P==0
   av=1; % 避免除以零
else
   av=P/nv; % 平均强度
end

% 计算缩放系数
% 这里的 120 是一个经验常数，用于将平均强度映射到灰度值 120 左右 (范围 0-255)
scale = sc*120/av;

% 应用指数变换和缩放
IMG = mli.^exp*scale;

% 截断超过 255 的值 (饱和处理)
IMG(IMG>255) = 255;
% 将小于 1 但非 0 的值设为 1 (避免过暗)
IMG(IMG<1&IMG~=0)=1;

% 显示图像
imagesc(IMG);
colormap gray % 使用灰度色图