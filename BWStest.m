function [H]= BWStest(Xarray,Yarray,Alpha)
%BWStest BWS 检验的简化加速版
%   计算双样本情况下的 Baumgartner-Weiß-Schindler 检验统计量。
%   用于判断两个样本是否来自同一分布 (即是否为同质像素)。
%
%   输入:
%   - Xarray: 高度 x 宽度 矩阵，每一列对应一个观测样本 (通常是参考像素的时间序列)
%   - Yarray: 高度 x 宽度 矩阵，每一列对应一个观测样本 (通常是待检测像素的时间序列)
%   - Alpha:  显著性水平，仅支持 0.05 (默认) 和 0.01。
%
%   输出:
%   - H:      0 => 在显著性水平 ALPHA 下接受原假设 (分布相同，即同质)
%             1 => 在显著性水平 ALPHA 下拒绝原假设 (分布不同，即异质)
%
%   参考文献:
%   [1] A Nonparametric Test for the General Two-Sample Problem, W. Baumgartner, P. Weiß and H. Schindler
%   Biometrics Vol. 54, No. 3 (Sep., 1998), pp. 1129-1135 
%
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University  
%   ======================================================================
%   06/2019 MJ: 修复了接受 nan 向量的 bug
%   ======================================================================

if nargin < 3 || isempty(Alpha)
    Alpha = 0.05;
end

if nargin < 2
    help BWStest
end

[n,m] = size(Xarray);
% 计算合并数据的秩 (tiedrank 处理并列秩)
ranks = tiedrank(cat(1,Xarray, Yarray));
% 分离并排序秩
xrank = sort(ranks(1:n,:));
yrank = sort(ranks((n+1):end,:));

% 计算统计量 B 的中间变量
% 对应论文中的公式计算权重和积分项的离散近似
temp  = (1:n)'*ones(1,m);
tempx = (xrank - 2.*temp).^2;
tempy = (yrank - 2.*temp).^2;
temp  = temp/(n+1).*(1-temp/(n+1))*2*n;
BX    = 1/n*sum(tempx./temp,1);
BY    = 1/n*sum(tempy./temp,1);

% 检验统计量 B
B = 1/2*(BX + BY);

% 确定临界值 b
if Alpha ==.05
    b = 2.493; % 参考文献 [1] 表 1
else %(Alpha =.01)
    b = 3.880; % 参考文献 [1] 表 1
end

% 假设检验
% 如果 B < b，则统计量在临界值内，接受原假设 (H0: 分布相同)
H = B < b; %to aviod nan vector
% 输出约定: 0 表示接受 (同质), 1 表示拒绝 (异质)
% 所以取反: 如果 B < b (True), H 变为 False (0) -> 接受
%          如果 B >= b (False), H 变为 True (1) -> 拒绝
H = ~H;
