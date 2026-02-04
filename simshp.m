function simshp(spnum,CalWin,sigma)
%simshp 用于验证同质像素选择算法的准确性
%   比较了四种方法: GLRT, KS, BWS 和 SHP (FaSHPS)。
%   结果与 JIANG et al., ISPRS, 2017 发表的结果相似。
%
%   用法:
%       simshp(spnum,CalWin,sigma);
%   
%
%   输入:
%   - spnum:    堆栈大小 (页数)
%   - CalWin:   图像窗口大小 [高度, 宽度]
%   - sigma:    无噪声幅值 (Rayleigh 分布参数的基准)
%
%
%
% 
%   参考文献:
%   [1] The potential of more accurate InSAR covariance matrix estimation for land cover mapping
%        Mi Jiang, Bin Yong, Xin Tian, Rakesh Malhotra, Rui Hu, Zhiwei Li, Zhongbo Yu and Xinxin Zhang,
%        ISPRS Journal of Photogrammetry and Remote Sensing Vol. 126,
%        April 2017, pp. 120-128.
%
%   
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University,

if nargin < 3 || isempty(sigma)
    sigma = 200;
end

if nargin < 2 || isempty(CalWin)
    CalWin = [15 15]; %img size
end

if nargin < 1 || isempty(spnum)
    spnum = 30;
end

ratio=1:.2:3; % 两个区域的参数比率，用于模拟异质性 (1 表示同质，>1 表示异质)
rep=10000; % 蒙特卡洛模拟重复次数
alpha =0.05; % 显著性水平

% Rayleigh 分布的真实参数
sigma=sqrt(sigma/2);
Matrix=sigma*ones(CalWin(1),CalWin(2)); 

InitRow=(CalWin(1)+1)/2; % InitRow 是中心行
InitCol=(CalWin(2)+1)/2; % InitCol 是中心列

%------------------------%
h1=zeros(length(ratio),rep); %GLRT 结果存储
h2=h1; %KS 结果存储
h3=h1; %BWS 结果存储
h4=h1; %SHP 结果存储
kstemp=zeros(CalWin(1),CalWin(2));

% 计算临界值 (Critical Values)
CR = GLRT_CR(spnum,alpha); % GLRT 的临界值
CRSHP=finv([alpha/2,1-alpha/2],2*spnum,2*spnum); % SHP (F-test) 的临界值

num=1;
p=1;
all = length(ratio)*rep;
all_step = floor(all/10);

% 遍历不同的比率 (模拟不同程度的异质性)
for ii=1:length(ratio)
    % 构造异质窗口：右半部分像素的参数乘以 ratio
    Matrix(:,InitCol+1:end)=ratio(ii)*sigma;
    NEWMatrix = repmat(Matrix,[1,1,spnum]);
    
    % 蒙特卡洛模拟循环
    for jj=1:rep
        % 生成 Rayleigh 分布噪声并调制信号
        Noise = random('rayl', 1,[CalWin(1),CalWin(2),spnum]);    
        NoiseAdd =Noise.*NEWMatrix;
        
        % 获取中心参考像素的时间序列
        Ref = NoiseAdd(InitRow,InitCol,:);
        
        % 1. GLRT (广义似然比检验)
        htemp = glrt(NoiseAdd,Ref,CR);
        h1(ii,jj)=sum(htemp(:)); % 记录拒绝原假设(认为是异质)的像素数量
        
        % 2. KS (Kolmogorov-Smirnov 检验)
        for ll=1:CalWin(1)
            for kk=1:CalWin(2)
                temp = NoiseAdd(ll,kk,:);
                % kstest2 比较两个样本是否来自同一分布
                kstemp(ll,kk)=kstest2(Ref(:),temp(:));
            end
        end
        h2(ii,jj)=sum(kstemp(:));
        
        % 3. BWS (Baumgartner-WeiB-Schindler 检验)
        Refarray = repmat(Ref(:),[1,CalWin(1)*CalWin(2)]);
        temparray= reshape(NoiseAdd,[CalWin(1)*CalWin(2),spnum])';
        bwstemp  = BWStest(Refarray,temparray,alpha);
        h3(ii,jj)=sum(bwstemp(:));       
        
        % 4. SHP (基于 FaSHPS 的方法)
        % 注意：这里传入的是强度 (幅度的平方)
        htemp = shp(NoiseAdd.^2,Ref.^2,CRSHP,alpha);
        h4(ii,jj)=sum(htemp(:));         
        
        num=num+1;
        if num == all_step * p;
            disp(['progress: ', num2str(10*p),'%']);
            p = p+1;
        end        
        
    end
end

% 计算拒绝率 (Rejection Rate) = 拒绝数量 / 总像素数
h1=h1/CalWin(1)/CalWin(2);
h2=h2/CalWin(1)/CalWin(2);
h3=h3/CalWin(1)/CalWin(2);
h4=h4/CalWin(1)/CalWin(2);

% 统计均值和标准差
meanh1=mean(h1,2);
stdh1 =std(h1,0,2);
meanh2=mean(h2,2);
stdh2 =std(h2,0,2);
meanh3=mean(h3,2);
stdh3 =std(h3,0,2);
meanh4=mean(h4,2);
stdh4 =std(h4,0,2);
figure;plot(ratio,meanh1,'diamond-',ratio,meanh2,'x--',ratio,meanh3,'*--',ratio,meanh4,'o-.');grid on;legend('GLRT','KS','BWS','SHP');ylabel('Mean rejection');xlabel('\sigma_1/\sigma_2');
figure;plot(ratio,stdh1,'diamond-',ratio,stdh2,'x--',ratio,stdh3,'*--',ratio,stdh4,'o-.');grid on;legend('GLRT','KS','BWS','SHP');ylabel('Std. rejection');xlabel('\sigma_1/\sigma_2');



%%
function y1 = glrt(NoiseAdd,ref,CR)
% GLRT 算法实现
[L,P,spnum] = size(NoiseAdd);

ref=repmat(ref,[L,P,1]);
% 计算参考像素和目标像素的二阶矩 (能量)
spest1=1/spnum*sum(ref.^2,3);
spest2=1/spnum*sum(NoiseAdd.^2,3);                
sigmamean=(spest1+spest2)/2;
% 计算 GLRT 统计量 T
T= 2*spnum*log(sigmamean)-spnum*log(spest1)-spnum*log(spest2);
% 判断是否大于临界值
y1 = T>CR;


%%
function y1 = shp(data,ref,CRSHP,Alpha)
% SHP (FaSHPS) 算法实现
[L,P,npages] = size(data);
LRT_nl = 3; %7*7 window size
LRT_nw = 3; 

InitRow=(L+1)/2; % InitRow is CenterRow
InitCol=(P+1)/2;

% Gamma 分布置信区间临界值
Galpha_L = gaminv(Alpha/2,npages,1);
Galpha_U = gaminv(1-Alpha/2,npages,1);

% 1. 初始估计 (LRT)
% 在小窗口内进行 F 检验
Matrix = data(InitRow-LRT_nl:InitRow+LRT_nl,InitCol-LRT_nw:InitCol+LRT_nw,:);
ref = mean(ref,3); % 参考像素的时间平均
temp = mean(Matrix,3); % 邻域像素的时间平均
T = ref./temp; % 均值比
T = T>CRSHP(1)&T<CRSHP(2); % F 检验
SeedPoint = mean(temp(T)); % 初始种子点均值

% 2. 迭代 (Gamma 置信区间)
MeanMatrix = mean(data,3);
% 检查所有像素是否在置信区间内
SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
y1 = ~SeedPoint; % 返回拒绝掩码 (1 表示非同质)

function C = GLRT_CR(imgnum,alpha)
%GLRT_CR 用于查找 Rayleigh 分布下双样本广义似然比检验 (GLRT) 的临界区域
%
%   用法:
%       C = GLRT_CR(imgnum,alpha);
%   
%
%   输入:
%   - imgnum:   样本数量 (时间点数)
%   - alpha:    显著性水平 (0~1)，默认为 0.05
%
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University,  

if nargin < 2
    alpha =.05;
end

if nargin < 1
    help GLRT_CR 
    return
end

sigma = 200/sqrt(2);
rep=10000;

% 模拟零假设下的分布 (两个样本来自同一分布)
a1=random('rayl',sigma,[imgnum,rep]);
a2=random('rayl',sigma,[imgnum,rep]);
spest1=1/imgnum*sum(a1.^2);
spest2=1/imgnum*sum(a2.^2);
sigmamean=(spest1+spest2)/2;
% 计算统计量分布
T= 2*imgnum*log(sigmamean)-imgnum*log(spest1)-imgnum*log(spest2);
% 获取分位数作为临界值
C=quantile(T,1-alpha);