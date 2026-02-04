function []=SHP_SelPoint_v2()
%SHP_SelPoint_v2 用于在 SAR 强度图像堆栈上选择同质像素 (SHP)
%   用法:
%       []=SHP_SelPoint_v2()
%   输入:
%       - param_list:  由 parmslist_generator 生成的参数文件
%
%   [1] Distributed scatterer interferometry with the refinement of spatiotemporal coherence
%        Mi Jiang, Andrea Monti-Guarnieri
%        IEEE Transactions on Geoscience and Remote Sensing vol. 58, no. 6, 
%        June 2020, pp. 3977-3987
%
%   [2] Fast Statistically Homogeneous Pixel Selection for Covariance Matrix Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding, Ramon F. Hanssen, Rakesh Malhotra and Ling Chang,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 53, no. 3,
%        March 2015, pp. 1213-1224.
% 
%   [3] The potential of more accurate InSAR covariance matrix estimation for land cover mapping
%        Mi Jiang, Bin Yong, Xin Tian, Rakesh Malhotra, Rui Hu, Zhiwei Li, Zhongbo Yu and Xinxin Zhang,
%        ISPRS Journal of Photogrammetry and Remote Sensing Vol. 126,
%        April 2017, pp. 120-128.
% 
% 
%   提示: 当使用多视 (N-looks) 强度图像时，只能使用 BWS 算法！
%
%
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University,  
%
%   ============================================================
%   11/2021 MJ change input/ouptput opinion
%   ============================================================


tic;

%读取参数文件
paramslistname = 'param_list';
if ~exist(paramslistname,'file')
    error('please prepare param_list file using parmslist_generator');
end

params = readDSIpro(paramslistname);

%输入参数解析
nlines  = params.nlines;
nwidths = params.nwidths;
mlipath = params.mli_path;
CalWin  = params.shp_calwin; %[方位向 距离向]
Alpha   = params.shp_alpha;
EstAgr  = params.shp_method;
disp('------------Parameter Info.------------');
disp(['Lines (azimuth):        ', num2str(nlines)]);
disp(['Widths (range):         ', num2str(nwidths)]);
disp(['SHP algorithm:          ', EstAgr]);
disp(['Sliding window size:    ', '[',num2str(CalWin(1)),' ',num2str(CalWin(2)),']']);
disp(['Alpha:                  ', num2str(Alpha)]);		
disp('--------------------------------------- ');

%输出文件名定义
shpftname = 'shp.ft'; % SHP 足迹 (索引) 文件
broname   = 'shp.bro';% 同质像素 (兄弟像素) 数量文件

% 读取强度图像堆栈
mlistack = ImgRead(mlipath,'mli',nlines,'float32'); % 读取后缀为 mli 的强度序列
mlistack = mlistack.datastack;
mlistack = single(mlistack);
[~,~,npages] = size(mlistack);

%参数准备:
RadiusRow=(CalWin(1)-1)/2; % 窗口行半径
RadiusCol=(CalWin(2)-1)/2; % 窗口列半径
InitRow=(CalWin(1)+1)/2; % InitRow 是中心行
InitCol=(CalWin(2)+1)/2; % InitCol 是中心列

%边缘镜像填充
mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');
meanmli = mean(mlistack,3); % 计算时间平均强度
[nlines_EP,nwidths_EP]= size(meanmli);
PixelInd=false(CalWin(1)*CalWin(2),nlines*nwidths); % 初始化结果矩阵

%估计 SHP
num=1;
p=1;
all = nlines*nwidths;
all_step = floor(all/10);

if strcmpi(EstAgr,'FaSHPS') 
    
    %设置参数统计的阈值
    LRT_nl = 3;
    LRT_nw = 3;
    if RadiusRow<LRT_nl
        LRT_nl=1;
    end
    if RadiusCol<LRT_nw
        LRT_nw=1;
    end

    %计算临界区 (Critical region)
    CR_lo = finv(Alpha/2,2*npages,2*npages);
    CR_up = finv(1-Alpha/2,2*npages,2*npages);
    Galpha_L = gaminv(Alpha/2,npages,1);
    Galpha_U = gaminv(1-Alpha/2,npages,1);    
    
    % 遍历每一个像素
    for kk=InitCol:nwidths_EP-RadiusCol
        for ll=InitRow:nlines_EP-RadiusRow       
            %初始估计 (似然比检验 Likelihood-ratio test)
            temp = meanmli(ll-LRT_nl:ll+LRT_nl,kk-LRT_nw:kk+LRT_nw);
            T = meanmli(ll,kk)./temp;
            T = T>CR_lo&T<CR_up;
            SeedPoint = mean(temp(T));
            %迭代 (Gamma 置信区间)
            MeanMatrix = meanmli(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol);
            SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %检查成员资格
            SeedPoint(InitRow,InitCol)=true;
            %连通性分析
            LL = bwlabel(SeedPoint); %double
            PixelInd(:,num)=LL(:)==LL(InitRow,InitCol);  
            num=num+1;
            if num == all_step * p
                disp(['progress: ', num2str(10*p),'%']);
                p = p+1;
            end
        end
    end
else %BWS (非参数统计)
    for kk=InitCol:nwidths_EP-RadiusCol
        for ll=InitRow:nlines_EP-RadiusRow
            Matrix = mlistack(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol,:);
            Ref = Matrix(InitRow,InitCol,:);
            T = BWStest(repmat(Ref(:),[1,CalWin(1)*CalWin(2)])...
                ,reshape(Matrix,[CalWin(1)*CalWin(2),npages])',Alpha);   
            temp=reshape(~T,[CalWin(1),CalWin(2)]);
            %连通性分析
            LL = bwlabel(temp);
            PixelInd(:,num)=LL(:)==LL(InitRow,InitCol);       
            num=num+1;
            if num == all_step * p
                disp(['progress: ', num2str(10*p),'%']);
                p = p+1;
            end              
        end
    end   
end

%生成 SHP 数量图            
BroNum = sum(PixelInd,1);
BroNum = uint16(reshape(BroNum(:),[nlines,nwidths]));    
          
%输出结果到文件
if exist(shpftname,'file')
    delete(shpftname);
end
fwritebkj(PixelInd',shpftname,'uchar','b'); % 写入 SHP 索引 (转置后写入)
if exist(broname,'file')
    delete(broname);
end
fwritebkj(BroNum,broname,'uint16','b'); % 写入 SHP 数量
t=toc;

figure;imagesc(BroNum);axis image off;colormap jet;
ti=title ('Homogeneous Pixel Number');colorbar;
set(ti,'fontweight','bold');
disp(['SHP_SelPoint operation completed in ',int2str(t/60),' min(s).']);
disp('Done!');

% =========================================================================
% SHP_SelPoint_v2 与 SHP_SelPoint (v1) 的主要区别说明
% =========================================================================
%
% 1. 输入方式 (Input):
%    - SHP_SelPoint (v1): 标准 MATLAB 函数，通过函数参数直接接收数据。
%      调用示例: [SHP] = SHP_SelPoint(mlistack, [15 15], 0.05, 'FaSHPS');
%      适用场景: 适合在代码中直接调用，进行交互式分析或调试。
%
%    - SHP_SelPoint_v2 (v2): 不接收输入参数，依赖外部参数文件 'param_list'。
%      调用示例: SHP_SelPoint_v2();
%      适用场景: 适合集成到自动化批处理流程中，通过文件系统读取配置。
%
% 2. 输出方式 (Output):
%    - SHP_SelPoint (v1): 将结果作为 MATLAB 结构体变量 (SHP) 返回到工作区。
%      包含字段: .PixelInd (索引), .BroNum (数量), .CalWin (窗口)。
%
%    - SHP_SelPoint_v2 (v2): 将结果直接写入磁盘文件，不返回变量。
%      输出文件: 'shp.ft' (SHP索引), 'shp.bro' (SHP数量)。
%      适用场景: 适合处理大数据量，避免占用过多内存，或供后续流程读取文件使用。
%
% 3. 总结:
%    - v1 更灵活，适合开发和研究。
%    - v2 更像是一个独立的“处理模块”，适合工程化和批处理。
% =========================================================================

