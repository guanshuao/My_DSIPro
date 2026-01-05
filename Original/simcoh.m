function Coh=simcoh(B,Lshp)
%This function is used to validate the accuracy of bias mitigation by means
%of non-parametric bootstrapping method. 
%
%
%   Usage:
%       Coh=simcoh(B,Lshp)
%   
%
%   Inputs:
%   - B:    The bootstrap replication
%   - Lshp: The SLC sample size 
%   Outputs:
%   - Coh:  A struct including coherence magnitude before and after bootstrapping correction    
%
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
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen University,

replication = 10000;
coh_true=0:.1:1;
UBCoh=zeros(replication,length(Lshp));
BCoh=UBCoh;
for jj=1:length(coh_true)  
    tic;
    for ii=1:replication
    %slc simulate
    s1 = (randn(Lshp,1) + 1j * randn(Lshp,1))/sqrt(2);
    x  = (randn(Lshp,1) + 1j * randn(Lshp,1))/sqrt(2);
    s2 = coh_true(jj).*s1 + sqrt(1 - coh_true(jj).^2).*x;
    InterfValue =s1.*conj(s2);
    %baised coherence
    BCoh(ii,jj)=abs(sum(InterfValue))./sqrt(sum(abs(s1).^2).*sum(abs(s2).^2));
    %Non-parametric Bootstrapping
    Idx   = randi(Lshp,Lshp,B);
    Bts1 = s1(Idx);
    Bts2 = s2(Idx);
    Btcoh=abs(sum(Bts1.*conj(Bts2),1)./sqrt(sum(abs(Bts1).^2,1).*sum(abs(Bts2).^2,1)));
    UBCoh(ii,jj) = 2*BCoh(ii,jj)-mean(Btcoh);     
    end   
    time=toc;
    fprintf('Processing coherence: %d / %d, time = %.0f sec\n',jj,length(coh_true),time);
end

Coh.UBCoh_mean = mean(UBCoh);
Coh.BCoh_mean = mean(BCoh);
figure;plot(coh_true,coh_true,'LineWidth',2);hold on;plot(coh_true,Coh.BCoh_mean,'--*');plot(coh_true,Coh.UBCoh_mean,'-.o');grid on
xlabel('True coherence');ylabel('Mean coherence magnitude');xlim([0,1]);ylim([0,1]);
legend('truth','bias','bootstrap')

figure;%std
hold on;plot(coh_true,std(BCoh),'--*');plot(coh_true,std(UBCoh),'-.o');grid on
RMSE_b = std(BCoh).^2+(Coh.BCoh_mean-coh_true).^2;
RMSE_b = sqrt(RMSE_b);
RMSE = std(UBCoh).^2+(Coh.UBCoh_mean-coh_true).^2;
RMSE = sqrt(RMSE);
figure;plot(coh_true,RMSE_b,'-.o');hold on;plot(coh_true,RMSE,'--*');


