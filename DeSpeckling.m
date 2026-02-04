function Mliimg = DeSpeckling(mlistack,SHP)
%DeSpeckling 对强度图像堆栈进行滤波 (去斑)
%   如果没有提供 SHP 文件，则返回 Boxcar 滤波后的强度图像。
%   支持多视 (N-looks) 强度图像。
%
%   用法:
%      Mliimg = DeSpeckling(mlistack,SHP)
%   
%   输入:
%   - mlistack: 高度 x 宽度 x 页数 的 (实数) 矩阵，例如 SAR 单视强度序列
%   - SHP:      (可选) 同质像素结构体，详见 "SHP_SelPoint.m" 脚本。
%               如果未提供，则默认使用 Boxcar 滤波。
%
%   输出:
%   - Mliimg:   滤波后的强度图像 
%
%   参考文献:
%   [1] Fast Statistically Homogeneous Pixel Selection for Covariance Matrix Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding, Ramon F. Hanssen, Rakesh Malhotra and Ling Chang,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 53, no. 3,
%        March 2015, pp. 1213-1224.
% 
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University,  
%   ===========================================
%   MJ 11/2021 移除 Acc 和高斯核
%   ===========================================

BoxFilt=false;

% 如果输入参数少于 2 个 (即没有提供 SHP)，则启用 Boxcar 滤波
if nargin < 2 || isempty(SHP)
    BoxFilt=true;
end

if nargin < 1
    help DeSpeckling
    return
end

tic;
[nlines,nwidths,npages]=size(mlistack);
Mliimg = mlistack; % 初始化输出

if BoxFilt
    % Boxcar 滤波 (均值滤波)
    h = fspecial('average',[7 7]); % 使用 7x7 的均值滤波器
    for ii=1:npages
        Mliimg(:,:,ii) = filter2(h,mlistack(:,:,ii));
        fprintf('BOXCAR DESPECKLING: %d / %d is finished...\n',ii,npages);
    end     
else
    % 自适应滤波 (基于 SHP)
    CalWin =SHP.CalWin;
    RadiusRow=(CalWin(1)-1)/2;
    RadiusCol=(CalWin(2)-1)/2;  
    % 边缘填充
    mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');
    
    %Despeckling (去斑)
    for ii=1:npages
        temp = mlistack(:,:,ii);
        num=1;
        for jj = 1:nwidths
            for kk= 1:nlines
                x_global  = jj+RadiusCol;
                y_global  = kk+RadiusRow;
                % 提取窗口数据
                MliValue  = temp(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                % 使用 SHP 索引筛选同质像素
                MliValue  = MliValue(SHP.PixelInd(:,num));
                % 计算均值作为滤波后的值
                Mliimg(kk,jj,ii) = mean(MliValue);
                num=num+1;
            end
        end         
        fprintf(' ADP. DESPECKLING: %d / %d is finished...\n',ii,npages);
    end     
end

t=toc;
disp(['DeSpeckling operation completed in ',int2str(t/60),' min(s).']);
disp('Done!'); 