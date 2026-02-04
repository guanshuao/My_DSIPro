function [Footprint]=SHP_FootPrint(imagename,SHP,PointIdx,imagesize_radius)
%SHP_FootPrint 显示给定参考像素的同质像素(SHP)的分布范围(足迹)
%
%   用法:
%       [Footprint]=SHP_FootPrint(imagename,SHP,PointIdx,imagesize_radius)
%   
%
%   输入:
%   - imagename:        底图路径，ras 格式 (通常是平均强度图)
%   - SHP:              SHP 数据结构，详见脚本 "SHP_SelPoint.m"
%   - PointIdx:         参考像素在底图中的坐标 [行 列]
%   - imagesize_radius: 显示背景的半径 [行 列]，默认值为 [25 25]
%
%   输出:
%   - Footprint:        返回参考点的坐标 [行 列]
%
%   示例:
%   从原始底图中手动选择一个参考点，显示半径为 25 像素的区域:
%    [Footprint]=SHP_FootPrint('/home/user/INSAR/COHEST/MLI/mli_ave.ras',SHP);
%   比较不同 SHP 集合的结果 (例如 BWS 算法的结果):
%    SHP_FootPrint('/home/user/INSAR/COHEST/MLI/mli_ave.ras',BWSSHP,Footprint);   
%
%
%   
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen Univeristy,

% 读取底图 (通常是平均强度图)
[mean_amp,amp_color]=imread(imagename);

% 如果未提供显示半径，默认为 [25 25]
if nargin < 4
    imagesize_radius=[25 25]; %[row col]
end

% 如果未提供参考点坐标，允许用户在图上交互式选择
if nargin < 3
    disp('select a point from basemap...')
    clf;
    figure(1);image(mean_amp);colormap(amp_color);axis image
    PointIdx=fliplr(round(ginput(1))); % ginput 返回 [x y] 即 [列 行]，fliplr 转换为 [行 列]
end    

% 如果未提供 SHP 数据，显示帮助信息
if nargin < 2
    help SHP_FootPrint
    return;
end

% 获取计算窗口大小
CalWin=SHP.CalWin;
[nlines,nwidths]=size(mean_amp);

% 计算参考点在图像中的线性索引
IND = sub2ind([nlines,nwidths],PointIdx(1),PointIdx(2));

% 从 SHP 结构中提取该点的 SHP 索引掩码，并重塑为窗口大小
% SHP.PixelInd 的每一列对应图像中的一个像素，行对应窗口内的像素
SHPS= reshape(SHP.PixelInd(:,IND), [CalWin(1),CalWin(2)]);

% 计算窗口半径
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;

% 生成以参考点为中心的窗口内的所有像素坐标网格 (BoxCar 矩形窗口)
[X,Y] = meshgrid(PointIdx(2)-RadiusCol:PointIdx(2)+RadiusCol,PointIdx(1)-RadiusRow:PointIdx(1)+RadiusRow);

% 保存矩形窗口(BoxCar)的所有坐标用于显示
footprintColREC = X(:);%footprint for REC
footprintRowREC = Y(:);

% 根据 SHP 掩码筛选出同质像素的坐标
X=X(SHPS);
Y=Y(SHPS);

% 移除超出图像边界的坐标点
negvalue = X<=0|Y<=0|X>nwidths|Y>nlines;
X=X(~negvalue);
Y=Y(~negvalue);

% 定义用于显示的裁剪窗口范围 (以参考点为中心，imagesize_radius 为半径)
ShowRow = imagesize_radius(1);%image size
ShowCol = imagesize_radius(2);
CropWin = [PointIdx(1)-ShowRow,PointIdx(1)+ShowRow,PointIdx(2)-ShowCol,PointIdx(2)+ShowCol];

% 确保裁剪窗口不超出图像边界
if CropWin(1) < 1
        CropWin(1) = 1;
end
if CropWin(2) > nlines
        CropWin(2) = nlines;
end
if CropWin(3) < 1;
        CropWin(3) = 1;
end
if CropWin(4) > nwidths
        CropWin(4) = nwidths;
end
Footprint = [PointIdx(1),PointIdx(2)];

% 绘图输出
figure; colormap(amp_color);

% 右图: 显示 SHP 分布
subplot(1,2,2);
% 显示裁剪后的底图区域
image(CropWin(3):CropWin(4),CropWin(1):CropWin(2),mean_amp(CropWin(1):CropWin(2),CropWin(3):CropWin(4)));
hold on;
plot(X,Y,'g.','MarkerSize',10); % 绘制 SHP 像素 (绿色点)
plot(PointIdx(2),PointIdx(1),'r.','MarkerSize',10); % 绘制参考中心点 (红色点)
title('SHP'); axis image

% 左图: 显示 BoxCar (矩形窗口) 分布
subplot(1,2,1);
% 显示裁剪后的底图区域
image(CropWin(3):CropWin(4),CropWin(1):CropWin(2),mean_amp(CropWin(1):CropWin(2),CropWin(3):CropWin(4)));
hold on;
plot(footprintColREC,footprintRowREC,'g.','MarkerSize',10); % 绘制整个矩形窗口内的像素 (绿色点)
plot(X,Y,'b.','MarkerSize',10); % 叠加绘制 SHP 像素 (蓝色点)
plot(PointIdx(2),PointIdx(1),'r.','MarkerSize',10); % 绘制参考中心点 (红色点)
title('BoxCar');axis image
