function []=DsC2Spat(datalistname)
%DsC2Spat 选择分布式散射体 (DS) 候选点，构建空间网络，并在雷达-多普勒坐标系下使用最短路径 (SP) 算法优化边缘
%   用法:
%       [] = DsC2Spat('C:\Users\DELL\Desktop\DSIpro\data_list');
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
%   [1] Distributed scatterer interferometry with the refinement of spatiotemporal coherence
%        Mi Jiang, Andrea Monti-Guarnieri
%        IEEE Transactions on Geoscience and Remote Sensing vol. 58, no. 6, 
%        June 2020, pp. 3977-3987
% 
%   [2] Sentinel-1 TOPS co-registration over low-coherence areas and its 
%       application to velocity estimation using the all pairs shortest path algorithm
%        Mi Jiang
%        Journal of Geodesy vol. 94 no. 10,
%        Oct.  2020, pp. 1-15
% 
% 
% 
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University, 


if nargin < 1
    help DsC2Spat
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
fpath   = params.pl_path; % Phase Linking (PL) 结果路径
if ~isempty(strmatch(fpath(end),filesep))
    fpath=fpath(1:end-1);
end
pcohname = 'pcoh.mat';
maskname = 'MASK.mat';
bronumname ='shp.bro';
datalist =load(datalistname);
bro_thre = params.pl_bro; %pixels
pcoh_thre= params.pl_pcoh;
dist_thre= params.apsp_dist; %m
pix_size = params.apsp_pixsize;
nlines  = params.nlines;
nwidths = params.nwidths;
az_spacing = params.az_spacing; %m
gr_spacing = params.rg_spacing/sind(params.inc); %ground spacing %m
disp('------------Parameter Info.------------');
disp(['Lines (azimuth):        ', num2str(nlines)]);
disp(['Widths (range):         ', num2str(nwidths)]);
disp(['Azimuth pixel spacing:  ', num2str(az_spacing)]);
disp(['Ground pixel spacing:   ', num2str(gr_spacing)]);
disp(['SHP number:             ', num2str(bro_thre)]);
disp(['Posteriori coherence:   ', num2str(pcoh_thre)]);
disp(['Resample size:          ', num2str(pix_size)]);
disp('--------------------------------------- ');

%输出文件名
dscname = 'ph.mat'; % 分布式散射体 (DS) 的相位 [n_ps, n_image]
ijname  = 'ij.mat'; % 分布式散射体的位置 [az rg]

%加载数据
%同质像素数量
bronum = freadbkj(bronumname,nlines,'uint16','b');

%后验相干性
pcoh=load(pcohname);
pcoh=pcoh.pcoh;

%掩膜 (MASK)
MASK = load(maskname);
MASK = MASK.MASK;

%阈值筛选 (选择 DS 候选点)
% 条件: 后验相干性 >= 阈值 AND 同质像素数量 >= 阈值 AND 在掩膜内
idx = pcoh>=pcoh_thre&bronum>=bro_thre&MASK;

%获取候选点的 ij 坐标
[X,Y] = meshgrid(1:nwidths,1:nlines);
ij=[Y(idx),X(idx)]; %[方位向 距离向] 
ij = single(ij);
clear Y X
pcoh = pcoh(idx); % 保留候选点的相干性

%读取 PL 相位 (Phase Linking 优化后的相位)
npages = size(datalist,1);
n_ps   = size(ij,1); 
ph     = ones(n_ps,npages,'single');
for ii=1:npages
    % 只有当主辅影像不同时才读取 (避免自干涉)
    if datalist(ii,2)~=datalist(1,1)
        plname = [fpath,filesep,num2str(datalist(ii,1)),'_',num2str(datalist(ii,2)),'.pl'];
        tmp    = freadbkj(plname,nlines,'cpxfloat32','b');
        ph(:,ii) = tmp(idx);
    end
end
clear tmp
save(dscname,'ph'); % 保存相位 (包含相干性信息，因为是复数)
save(ijname,'ij');  % 保存坐标

%通过最大化后验相干性重采样到更大的网格 (降采样以减少计算量)
if ~isempty(pix_size)

    rY = pix_size/az_spacing; % 方位向重采样比率
    rX = pix_size/gr_spacing; % 距离向重采样比率

    f2Y    = floor(nlines/rY);
    f2X    = floor(nwidths/rX);
    idx    = zeros(f2Y,f2X);

    startY = 1;
    for ii=1:f2Y
        startX = 1;
      for jj=1:f2X
          % 查找当前大网格内的所有 DS 候选点
          ix = find(ij(:,1)>=startY&ij(:,1)<ii*rY&ij(:,2)>=startX&ij(:,2)<jj*rX);
          % 在网格内选择后验相干性最大的点作为代表
          if ~isempty(ix)
              [~, local_idx_best_one] = max(pcoh(ix));
              idx(ii,jj) = ix(local_idx_best_one);
          end
          startX = startX+rX;

      end
      startY = startY+rY;
    end
    idx=idx(:);
    idx(idx==0)=[]; % 移除空网格
    ij = ij(idx,:); % 更新坐标
    ph = ph(idx,:); % 更新相位

end

fprintf('Refine spatial netwrok using SP algorithm.\n');

%构建边缘集 (Edge set) - 全连接网络 (但在距离阈值内)
[xx,yy] = meshgrid(1:size(ij,1),1:size(ij,1));
xx = triu(xx,1); % 上三角矩阵，避免重复边和自环
yy = triu(yy,1);
edgs = [yy(yy~=0),xx(xx~=0)];
edgs = sortrows(edgs,1);
% 计算所有边的物理距离
all_edgs = sqrt((ij(edgs(:,1),1)-ij(edgs(:,2),1)).^2*az_spacing^2 + ...
    (ij(edgs(:,1),2)-ij(edgs(:,2),2)).^2*gr_spacing^2);
% 仅保留距离小于阈值的边
edgs = edgs(all_edgs < dist_thre,:);

% 归一化相位 (只保留相位信息，幅度设为 1)
ix = ph~=0;
ph(ix) = ph(ix)./abs(ph(ix));  

%构建三角网 (TIN) 作为参考网络
tr = delaunayTriangulation(double(ij(:,1)),double(ij(:,2)));      
tinedgs = edges(tr); %[起点 终点]
% 计算 TIN 边的距离并筛选
all_tinedgs = sqrt((ij(tinedgs(:,1),1)-ij(tinedgs(:,2),1)).^2*az_spacing^2 + ...
                          (ij(tinedgs(:,1),2)-ij(tinedgs(:,2),2)).^2*gr_spacing^2);
tinedgs = tinedgs(all_tinedgs < dist_thre,:);

%计算 TIN 边的时间相干性 (仅作为演示，未进行搜索解算)
tempcoh_tin = ph(tinedgs(:,2),:).*conj(ph(tinedgs(:,1),:)); % 差分相位
tempcoh_tin = tempcoh_tin.';
tempcoh_tin = abs(sum(tempcoh_tin))./sum(abs(tempcoh_tin)); % 时间相干性
% 构建 TIN 图，权重为 -10*log10(相干性)，相干性越高权重越小 (距离越短)
GTIN = graph(tinedgs(:,1),tinedgs(:,2),-10*log10(tempcoh_tin'),size(ij,1));    

%绘制 TIN 网络
figure;
p = plot(GTIN,'XData',ij(:,2),'YData',ij(:,1),'EdgeCData',tempcoh_tin','EdgeColor','flat');colormap jet
p.MarkerSize = 1;grid on;
caxis([0,1]);    
xlim([min(ij(:,2))-20,max(ij(:,2))+20]);ylim([min(ij(:,1))-20,max(ij(:,1))+20]);
title('(a) Reference network before SP');caxis([0,1]);colorbar;axis ij

%使用 SP (最短路径) 算法优化 TIN 边的路径
% 在全连接图 (距离限制内) 中寻找连接 TIN 边两个端点的最短路径
num=1;
startpath=[];
endpath=[];
% 计算全连接图 (edgs) 的时间相干性
tempcoh = ph(edgs(:,2),:).*conj(ph(edgs(:,1),:));
tempcoh = tempcoh.';
tempcoh = abs(sum(tempcoh))./sum(abs(tempcoh));
% 构建全连接图 G
G=graph(edgs(:,1),edgs(:,2),-10*log10(tempcoh')); 

% 对每一条 TIN 边，在 G 中寻找最短路径替代它
for jj=1:size(tinedgs,1)
    [tmp]  = shortestpath(G,tinedgs(jj,1),tinedgs(jj,2));
    tmpend = num + length(tmp) - 2;
    startpath(num:tmpend) = tmp(1:end-1);
    endpath(num:tmpend)   = tmp(2:end);
    num = tmpend + 1;
end    

% 构建优化后的网络 (SP Network)
qfrom = min(startpath,endpath);
qto   = max(startpath,endpath);
ft    = [qfrom.', qto.'];
ft    = unique(ft,'rows');% 移除重复边
[~,Lia] = ismember(ft,edgs,'rows'); % 找到优化后的边在原边集中的索引
GSP=graph(ft(:,1),ft(:,2),-10*log10(tempcoh(Lia)),size(ij,1));

%绘制 SP 网络
figure;
p = plot(GSP,'XData',ij(:,2),'YData',ij(:,1),'EdgeCData',tempcoh(Lia)','EdgeColor','flat');
colormap(jet(65));
p.MarkerSize = 1;
caxis([0,1]);
grid on;
xlim([min(ij(:,2))-20,max(ij(:,2))+20]);ylim([min(ij(:,1))-20,max(ij(:,1))+20]);
title('(b) Reference network after SP');axis ij;caxis([0,1]);colorbar;

%对比 TIN 和 SP 网络的相干性分布
figure;
htin = histogram(tempcoh_tin(:));
htin.Normalization = 'probability';
htin.BinWidth = 0.01;
htin.DisplayStyle='stairs';
hold on;
hapsp = histogram(tempcoh(Lia));
hapsp.Normalization = 'probability';
hapsp.BinWidth = 0.01;
hapsp.DisplayStyle='stairs';
xlim([0,1]);xlabel('Coherence');ylabel('Probability');legend('TIN','SP')
grid on

t=toc;    
disp(['DsC2Spat operation completed in ',int2str(t/60),' min(s).']);
disp('Done!');