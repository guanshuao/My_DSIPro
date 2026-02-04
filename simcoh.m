function Coh=simcoh(B,Lshp)
%simcoh 用于验证非参数 Bootstrapping 方法在减少相干性估计偏差方面的准确性
%
%
%   用法:
%       Coh=simcoh(B,Lshp)
%   
%
%   输入:
%   - B:    Bootstrap 重采样次数 (Replication)
%   - Lshp: SLC 样本大小
%   输出:
%   - Coh:  一个结构体，包含 Bootstrapping 校正前后的相干性幅值统计结果
%
%   参考文献:
%   [1] Fast Statistically Homogeneous Pixel Selection for Covariance Matrix Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding, Ramon F. Hanssen, Rakesh Malhotra and Ling Chang,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 53, no. 3,
%        March 2015, pp. 1213-1224.
% 
%   [2] Hybrid Approach for Unbiased Coherence Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding and Zhiwei Li,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 52, no. 5,
%        May  2014, pp. 2459-2473.
%
%   
%   本工具箱仅供研究使用，请在任何产生的出版物中引用上述论文。
%
%   Mi JIANG, Sun Yat-sen University,

replication = 10000; % 蒙特卡洛模拟的重复次数，用于统计平均性能
coh_true=0:0.1:1;     % 真实的相干性值，从 0 到 1，步长 0.1
UBCoh=zeros(replication,length(Lshp)); % 初始化去偏后的相干性矩阵 (注意：这里 length(Lshp) 似乎应该是 length(coh_true)，除非 Lshp 是个向量？看代码逻辑 Lshp 应该是标量样本数，这里初始化大小可能有点问题，或者意图是针对不同的 Lshp? 但下面循环是用 length(coh_true))
% 修正理解：如果 Lshp 是标量，length(Lshp) 为 1。如果代码原本是想存不同真值的模拟结果，应该是 zeros(replication, length(coh_true))。
% 让我们看循环：BCoh(ii,jj) = ... 其中 jj 是 coh_true 的索引。所以列数应该是 length(coh_true)。
% 原代码: UBCoh=zeros(replication,length(Lshp)); 
% 如果 Lshp 是标量，这只有 1 列。但下面循环 jj 跑到 11 (0:0.1:1)。这会触发自动扩容，或者报错。
% 假设 Lshp 是标量输入。这里可能是个笔误，或者 MATLAB 会自动扩展数组。
% 既然我是加注释，我保留原代码，但加上注释说明。

BCoh=UBCoh; % 初始化有偏相干性矩阵

% 遍历每一个设定的真实相干性值
for jj=1:length(coh_true)  
    tic;
    % 进行多次蒙特卡洛模拟
    for ii=1:replication
    % 1. SLC 数据模拟 (Simulate SLC data)
    % 生成两个复高斯随机信号 s1 和 x
    s1 = (randn(Lshp,1) + 1j * randn(Lshp,1))/sqrt(2);
    x  = (randn(Lshp,1) + 1j * randn(Lshp,1))/sqrt(2);
    % 根据真实相干性 coh_true(jj) 生成相关的信号 s2
    % s2 由 s1 的相关部分和 x 的不相关部分组成
    s2 = coh_true(jj).*s1 + sqrt(1 - coh_true(jj).^2).*x;
    
    InterfValue =s1.*conj(s2); % 干涉图 (共轭相乘)
    
    % 2. 计算有偏相干性 (Biased coherence estimation)
    % 标准的样本相干性计算公式
    BCoh(ii,jj)=abs(sum(InterfValue))./sqrt(sum(abs(s1).^2).*sum(abs(s2).^2));
    
    % 3. 非参数 Bootstrapping (Non-parametric Bootstrapping)
    % 生成 B 组重采样索引，范围 1~Lshp，有放回抽样
    Idx   = randi(Lshp,Lshp,B); 
    Bts1 = s1(Idx); % 重采样 s1
    Bts2 = s2(Idx); % 重采样 s2
    
    % 计算 B 次 Bootstrap 样本的相干性
    % sum(..., 1) 对每一列(每一次重采样)求和
    Btcoh=abs(sum(Bts1.*conj(Bts2),1)./sqrt(sum(abs(Bts1).^2,1).*sum(abs(Bts2).^2,1)));
    
    % 计算去偏相干性 (Unbiased coherence)
    % 公式: 2 * 原始估计 - Bootstrap均值
    UBCoh(ii,jj) = 2*BCoh(ii,jj)-mean(Btcoh);     
    end   
    time=toc;
    fprintf('Processing coherence: %d / %d, time = %.0f sec\n',jj,length(coh_true),time);
end

% 统计结果
Coh.UBCoh_mean = mean(UBCoh); % 去偏相干性均值
Coh.BCoh_mean = mean(BCoh);   % 有偏相干性均值

% 绘图 1: 均值对比
figure;
plot(coh_true,coh_true,'LineWidth',2); % 理想的 1:1 线
hold on;
plot(coh_true,Coh.BCoh_mean,'--*');    % 有偏估计均值
plot(coh_true,Coh.UBCoh_mean,'-.o');   % 去偏估计均值
grid on
xlabel('True coherence (真实相干性)');
ylabel('Mean coherence magnitude (平均相干性幅值)');
xlim([0,1]);ylim([0,1]);
legend('truth','bias','bootstrap')

% 绘图 2: 标准差对比
figure;%std
hold on;
plot(coh_true,std(BCoh),'--*');  % 有偏估计标准差
plot(coh_true,std(UBCoh),'-.o'); % 去偏估计标准差
grid on
% title('Standard Deviation (标准差)');
% legend('bias','bootstrap');

% 计算均方根误差 (RMSE)
% RMSE = sqrt(Bias^2 + Variance) = sqrt((Mean - True)^2 + Std^2)
RMSE_b = std(BCoh).^2+(Coh.BCoh_mean-coh_true).^2;
RMSE_b = sqrt(RMSE_b); % 有偏估计 RMSE

RMSE = std(UBCoh).^2+(Coh.UBCoh_mean-coh_true).^2;
RMSE = sqrt(RMSE);     % 去偏估计 RMSE

% 绘图 3: RMSE 对比
figure;
plot(coh_true,RMSE_b,'-.o');
hold on;
plot(coh_true,RMSE,'--*');
% title('RMSE (均方根误差)');
% legend('bias RMSE','bootstrap RMSE');
% xlabel('True coherence');