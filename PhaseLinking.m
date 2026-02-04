function [] = PhaseLinking(datalistname)
%PhaseLinking 单视条件下的干涉相位优化 (Phase Linking)
% 
%   用法:
%       [] = PhaseLinking('C:\Users\DELL\Desktop\DSIpro\data_list');
%
%   输入:
%   - datalistname: 包含 2 列 (主影像 辅影像) 的文件，例如:
%                   20170526 20170213
%                   20170526 20170219
%                   20170526 20170303
%                   20170526 20170315
%                   ...
%   - param_list:  由 parmslist_generator 生成的参数文件 (在代码内部读取)
% 
%   参考文献:
%   [1] Distributed Scatterer Interferometry with The Refinement of Spatiotemporal Coherence, 
%       Jiang Mi, Monti-Guarnieri Andrea. IEEE Transactions on Geoscience and Remote Sensing, 58(6): 3977-3987, 2020.
% 
%
% 
% 
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University  
%
%   ======================================================================
%   10/2018 MJ: 添加 iw (强度加权相干性估计) 选项
%   01/2019 MJ: 添加全网络平均相干性
%   04/2019 MJ: 移除 coh_opt 和 iw 选项
%   04/2019 MJ: 修复整个堆栈中像素为 0 的 bug
%   12/2020 MJ: 移除 SHP 文件，将 SHP 选择和 PL 结合 (注：代码中似乎仍读取 shp.ft)
%   01/2021 MJ: 将单主影像相干性序列添加到 *.pl 文件中
%   ======================================================================


if nargin < 1
    help PhaseLinking
    return;
end

tic;

%读取参数文件
paramslistname = 'param_list';
if ~exist(paramslistname,'file')
    error('please prepare param_list file using parmslist_generator');
end
params = readDSIpro(paramslistname);

%输入参数
fpath   = params.diff_path; % 差分干涉图路径
if ~isempty(strmatch(fpath(end),filesep))
    fpath=fpath(1:end-1);
end
plpath  = params.pl_path; % 输出路径
if ~isempty(strmatch(plpath(end),filesep))
    plpath=plpath(1:end-1);
end
nlines  = params.nlines;
nwidths = params.nwidths;
CalWin  = params.shp_calwin; %[方位向 距离向]
blksize = params.pl_blksize; % 分块处理的块数
shpname = 'shp.ft'; % SHP 索引文件 (由 SHP_SelPoint 生成)
datalist= load(datalistname); 

%输出文件名
networkcohname = 'fullnetworkcoh.mat'; % 全网络平均相干性
pcohname = 'pcoh.mat'; % 后验相干性 (拟合优度)
maskname = 'MASK.mat'; % 掩膜文件

%创建输出文件夹
if exist(plpath,'dir')
    rmdir(plpath,'s');
end
mkdir(plpath);

%分块处理设置 (Block processing)
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;  
% 计算每个块的列索引范围
widthidx = floor(linspace(1,nwidths,blksize+1));
% 如果块太小 (小于窗口半径)，则强制不分块
if min(diff(widthidx)<=RadiusCol)
    blksize =1;
end
fprintf('Estimate Covariance Matrix (%d block)...\n',blksize);    
block_size=zeros(blksize,2);
for i=1:blksize
    if i~=blksize
        block_size(i,:)=[widthidx(i),widthidx(i+1)-1];
    else
        block_size(i,:)=[widthidx(i),nwidths];
    end
end
% 计算读取 SHP 文件时的线性索引范围
shpind(:,1) = sub2ind([nlines,nwidths], ones(blksize,1), block_size(:,1));
shpind(:,2) = sub2ind([nlines,nwidths], ones(blksize,1)*nlines, block_size(:,2)); 

npages = size(datalist,1);
ref_idx = datalist(:,2)==datalist(1,1);   % 找到时间序列中主影像的索引 (通常是第一个)
pcoh = zeros(nlines,nwidths,'single');% 后验相干性矩阵
temp_mask = true(npages,npages);
temp_mask = triu(temp_mask,1); % 上三角掩膜，用于计算后验相干性
fullnetworkcoh = zeros(npages,npages); % 全网络相干性矩阵
% 掩膜文件初始化
MASK = true(nlines,nwidths);

num=1;
p=1;
all = nlines*nwidths;
all_step = floor(all/10);

% 开始分块循环
for ii=1:blksize
    
    % 确定当前块读取数据的范围 (包含边缘重叠区)
    temp_start = block_size(ii,1)-RadiusCol;
    temp_end   = block_size(ii,2)+RadiusCol;
    if temp_start < 1
        temp_start = 1;
    end
    if temp_end > nwidths
        temp_end = nwidths;
    end
    
    % 读取差分干涉图堆栈
    imgstack = ones(nlines,temp_end-temp_start+1,npages,'single');
    for jj=1:npages 
        if datalist(jj,2)~=datalist(1,1)
            diffname = [fpath,filesep,num2str(datalist(jj,1)),'_',num2str(datalist(jj,2)),'.diff'];
            % 读取指定列范围的数据
            imgstack(:,:,jj)   = freadbkj(diffname,nlines,'cpxfloat32','b',1,nlines,temp_start,temp_end);  %master - slave m-1, m-2, m-n_ifg 
        end       
    end  
    imgstack = conj(imgstack); % 取共轭，变为 slave-master: 1-m, 2-m,...,n_ifg-m
    ix=imgstack~=0;
    imgstack(ix)=imgstack(ix)./abs(imgstack(ix)); % 归一化幅度，只保留相位信息
    clear ix
    
    % 当前块的实际宽度 (不含边缘)
    local_width = block_size(ii,2) - block_size(ii,1)+1;

    % 掩膜处理：如果在整个时间序列中像素值都为 0，则标记为无效
    tmp = sum(imgstack,3)~=0;
    MASK(:,temp_start:temp_end)=logical(MASK(:,temp_start:temp_end).*tmp);

    % 读取 SHP 足迹 (SHP 索引)
    % 注意：这里假设 shp.ft 文件存储的是每个像素在窗口内的 SHP 掩码
    SHP = logical(freadbkj(shpname,nlines*nwidths,'uchar','b',shpind(ii,1),shpind(ii,2),1,CalWin(1)*CalWin(2)));
    
    % 掩膜处理：如果像素没有兄弟节点 (SHP 数量为 0 或 1)，则标记为无效
    tmp = reshape(sum(SHP,2),nlines,local_width);
    tmp = tmp ~= 1;
    MASK(:,block_size(ii,1):block_size(ii,2))=logical(MASK(:,block_size(ii,1):block_size(ii,2)).*tmp);
    clear tmp
    
    % 边缘填充 (处理块边缘和图像边缘)
    imgstack = padarray(imgstack,[RadiusRow 0],'symmetric');
    if temp_start==1
        imgstack = padarray(imgstack,[0 RadiusCol],'symmetric','pre');
    end
    if temp_end == nwidths
        imgstack = padarray(imgstack,[0 RadiusCol],'symmetric','post');
    end
    
    % 临时变量
    num_shp = 1;
    % 优化后的相位
    optphase = zeros(nlines,local_width,npages,'single');
    % 相干性
    coh_pl = optphase;
    
    % Phase Linking 核心循环
    for jj=1:local_width
        for kk=1:nlines
            if MASK(kk,jj+block_size(ii,1)-1)
                
                x_patch = jj+RadiusCol;
                y_patch = kk+RadiusRow;   
                
                % 1. 协方差矩阵估计 (Covariance Matrix, SCM)
                % 提取窗口内的数据
                Z = imgstack(y_patch-RadiusRow:y_patch+RadiusRow,x_patch-RadiusCol:x_patch+RadiusCol,:);
                Z = reshape(Z,[CalWin(1)*CalWin(2),npages]).';
                % 使用 SHP 筛选同质像素
                Z = Z(:,SHP(num_shp,:));
                % 计算样本协方差矩阵 Z*Z' / L
                CpxCoh = double(Z*Z'./size(Z,2));
                
                % 相干性矩阵 (幅度)
                Coh = abs(CpxCoh);
                
                % 累加全网络相干性 (用于统计)
                fullnetworkcoh = fullnetworkcoh + Coh;
                
                % 保存单主影像时间序列相干性图 (相对于参考影像)
                coh_pl(kk,jj,:) = Coh(:,ref_idx);
                
                % 2. 正定性检测与修正
                % 确保协方差矩阵是正定的，以便求逆
                [~,r] = chol(Coh); 
                e=1e-6;
                while ~(r == 0 && rank(Coh) == npages)
                    Coh = Coh+eye(npages)*e; % 对角线加载 (Diagonal Loading)
                    [~,r] = chol(Coh);
                    e=2*e;
                end   
                
                % 3. 相位优化 (Phase Optimization)
                % 求解特征值问题: inv(|C|) * C * v = lambda * v
                % 对应于最大似然估计 (MLE) 的近似解 (EMI 算法)
                [V,~]=eig(inv(Coh).*CpxCoh);
                V=V(:,1);% 取对应于最小特征值的特征向量 (注意：MATLAB eig 返回的特征值通常是排序的，这里取第一个可能是最小的，需确认 EMI 算法细节，通常是最大特征值对应信号子空间，但这里涉及求逆，可能是最小特征值)
                % 实际上，EMI 算法通常求解 coherence matrix 的最大特征向量。
                % 但这里是 inv(Coh).*CpxCoh，这类似于白化处理。
                
                % 设置参考相位为 0 (归一化)
                optphase(kk,jj,:)=V/V(ref_idx);
                
                % 4. 后验相干性 (Posteriori Coherence) / 拟合优度
                % 计算优化后相位与原始观测相位的拟合程度
                rho = exp(1j*(angle(CpxCoh)-angle(V*V'))); 
                rho = rho(temp_mask); % 取上三角部分
                rho = sum(rho); 
                pcoh(kk,jj+block_size(ii,1)-1) = rho;
                          
            end
            num=num+1;num_shp=num_shp+1;
            if num == all_step * p
                disp(['phase linking process: ', num2str(10*p),'%']);
                p = p+1;
            end                  
        end
    end
    % 输出结果
    % PL 部分
    fprintf('Output estimated %d PL phases for Block: %d \n', npages,ii); 
    % 将优化后的相位与相干性幅度结合，保存为复数
    optphase = coh_pl.*exp(-1j*angle(optphase)); 
    for ll = 1:npages
        if blksize==1
            filename=[plpath,filesep,num2str(datalist(ll,1)),'_',num2str(datalist(ll,2)),'.pl'];
        else
            % 如果分块，先保存为临时文件
            filename=[plpath,filesep,num2str(datalist(ll,1)),'_',num2str(datalist(ll,2)),'_',num2str(ii),'.pl'];
        end
        fwritebkj(optphase(:,:,ll), filename,'cpxfloat32','b');
    end

end

clear optphase imgstack

% 输出统计结果
% 后验相干性归一化 (0-1)
pcoh = real(pcoh)*2/npages/(npages-1); 
save(pcohname,'pcoh');
% 全网络平均相干性归一化
fullnetworkcoh = fullnetworkcoh/sum(MASK(:));
%N(N-1)/2 mean coherence in spatial
save(networkcohname,'fullnetworkcoh');
save(maskname,'MASK');

% 合并分块的 PL 相位文件
if blksize~=1
    tmp = zeros(nlines,nwidths,'single');
    for ll=1:npages
        filename = [plpath,filesep,num2str(datalist(ll,1)),'_',num2str(datalist(ll,2))];
        for ii=1:blksize
            % 读取各个块
            tmp(:,block_size(ii,1):block_size(ii,2))=freadbkj([filename,'_',num2str(ii),'.pl'],nlines,'cpxfloat32','b');
        end
        % 写入合并后的文件
        fwritebkj(tmp,[filename,'.pl'],'cpxfloat32','b');
        % 删除临时文件
        delete([filename,'*_*.pl']);
    end
end

figure;imagesc(pcoh,[0,1]);
ti=title ('Posteriori coherence');
set(ti,'fontweight','bold');colormap jet;colorbar


t=toc;    
disp(['PhaseLinking operation completed in ',int2str(t/60),' min(s).']);
disp('Done!');
      