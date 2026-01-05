function simshp(spnum,CalWin,sigma)
%This function is used to validate the accuracy of homogeneous pixel selection algorithm
%Four methods: GLRT,KS,BWS and SHP are given, the results are similar with
%that published in JIANG et al., ISPRS, 2017.
%
%   Usage:
%       simulation(spnum,CalWin,sigma);
%   
%
%   Inputs:
%   - spnum:    The stack size (page)
%   - CalWin:   The image size (height,width)
%   - sigma:    Noise-free amplitude
%
%
%
% 
%   [1] The potential of more accurate InSAR covariance matrix estimation for land cover mapping
%        Mi Jiang, Bin Yong, Xin Tian, Rakesh Malhotra, Rui Hu, Zhiwei Li, Zhongbo Yu and Xinxin Zhang,
%        ISPRS Journal of Photogrammetry and Remote Sensing Vol. 126,
%        April 2017, pp. 120-128.
%
%   
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen University,

if nargin < 3
    sigma =200;
end

if nargin < 2
    CalWin =[15 15]; %img size
end

if nargin < 1
    spnum =30;
end

ratio=1:.2:3;
rep=10000; %replication
alpha =0.05;

%the true parameter of the Rayleigh distribution
sigma=sqrt(sigma/2);
Matrix=sigma*ones(CalWin(1),CalWin(2)); 

InitRow=(CalWin(1)+1)/2; % InitRow is CenterRow
InitCol=(CalWin(2)+1)/2;

%------------------------%
h1=zeros(length(ratio),rep); %GLRT
h2=h1; %KS
h3=h1; %BWS
h4=h1; %SHP
kstemp=zeros(CalWin(1),CalWin(2));
CR = GLRT_CR(spnum,alpha);
CRSHP=finv([alpha/2,1-alpha/2],2*spnum,2*spnum);
num=1;
p=1;
all = length(ratio)*rep;
all_step = floor(all/10);
for ii=1:length(ratio)
    Matrix(:,InitCol+1:end)=ratio(ii)*sigma;
    NEWMatrix = repmat(Matrix,[1,1,spnum]);
    for jj=1:rep
        Noise = random('rayl', 1,[CalWin(1),CalWin(2),spnum]);    
        NoiseAdd =Noise.*NEWMatrix;
        Ref = NoiseAdd(InitRow,InitCol,:);
        %GLRT
        htemp = glrt(NoiseAdd,Ref,CR);
        h1(ii,jj)=sum(htemp(:));
        %KS
        for ll=1:CalWin(1)
            for kk=1:CalWin(2)
                temp = NoiseAdd(ll,kk,:);
                kstemp(ll,kk)=kstest2(Ref(:),temp(:));
            end
        end
        h2(ii,jj)=sum(kstemp(:));
        %BWS
        Refarray = repmat(Ref(:),[1,CalWin(1)*CalWin(2)]);
        temparray= reshape(NoiseAdd,[CalWin(1)*CalWin(2),spnum])';
        bwstemp  = BWStest(Refarray,temparray,alpha);
        h3(ii,jj)=sum(bwstemp(:));       
        %SHP
        htemp = shp(NoiseAdd.^2,Ref.^2,CRSHP,alpha);
        h4(ii,jj)=sum(htemp(:));         
        
        num=num+1;
        if num == all_step * p;
            disp(['progress: ', num2str(10*p),'%']);
            p = p+1;
        end        
        
    end
end
h1=h1/CalWin(1)/CalWin(2);
h2=h2/CalWin(1)/CalWin(2);
h3=h3/CalWin(1)/CalWin(2);
h4=h4/CalWin(1)/CalWin(2);


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

[L,P,spnum] = size(NoiseAdd);

ref=repmat(ref,[L,P,1]);
spest1=1/spnum*sum(ref.^2,3);
spest2=1/spnum*sum(NoiseAdd.^2,3);                
sigmamean=(spest1+spest2)/2;
T= 2*spnum*log(sigmamean)-spnum*log(spest1)-spnum*log(spest2);
y1 = T>CR;


%%
function y1 = shp(data,ref,CRSHP,Alpha)
[L,P,npages] = size(data);
LRT_nl = 3; %7*7 window size
LRT_nw = 3; 

InitRow=(L+1)/2; % InitRow is CenterRow
InitCol=(P+1)/2;

Galpha_L = gaminv(Alpha/2,npages,1);
Galpha_U = gaminv(1-Alpha/2,npages,1);

%Initial estimation (LRT)
Matrix = data(InitRow-LRT_nl:InitRow+LRT_nl,InitCol-LRT_nw:InitCol+LRT_nw,:);
ref = mean(ref,3);
temp = mean(Matrix,3);
T = ref./temp;
T = T>CRSHP(1)&T<CRSHP(2);
SeedPoint = mean(temp(T));
%iteration (Gamma Confidence interval)
MeanMatrix = mean(data,3);
SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
y1 = ~SeedPoint;

function C = GLRT_CR(imgnum,alpha)
%This function is used to find the critical region of two-sample Generalized Likelihood
%Ratio Test (GLRT) under Rayleigh distribution
%
%   Usage:
%       C = GLRT_CR(imgnum,alpha);
%   
%
%   Inputs:
%   - imgnum:   the sample number 
%   - alpha:    A value between 0 and 1 specifying the
%               significance level. Default is 0.05 for 5% significance.
%
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
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

a1=random('rayl',sigma,[imgnum,rep]);
a2=random('rayl',sigma,[imgnum,rep]);
spest1=1/imgnum*sum(a1.^2);
spest2=1/imgnum*sum(a2.^2);
sigmamean=(spest1+spest2)/2;
T= 2*imgnum*log(sigmamean)-imgnum*log(spest1)-imgnum*log(spest2);
C=quantile(T,1-alpha);

