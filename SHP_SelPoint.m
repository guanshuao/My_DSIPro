function [SHP]=SHP_SelPoint(mlistack,CalWin,Alpha,EstAgr)
%This function is used to select homogeneous pixels (SHP) on SAR intensity stack
%   Usage:
%       [SHP]=SHP_SelPoint(mlistack,CalWin,Alpha,EstAgr)
%   
%
%   Inputs:
%   - mlistack: A height by width by page matrix
%   - CalWin:   Fixed boxcar window size
%   - Alpha:    A value between 0 and 1 specifying the
%               significance level. Default is 0.05 for 5% significance.
%   - EstAgr:   Pixel selection algorithm: 1) FaSHPS; 2) BWS (support Alpha
%   =.05 and .01 only)
%   Outputs:
%   - SHP.PixelInd: A CalWin(1)*CalWin(2) by size(mlistack,1)*size(mlistack,2) array with elements of type logical, containing a SHPs set per pixel 
%   - SHP.BroNum:   The SHP number per pixel (reference pixel is not included) 
%   - SHP.CalWin:   Fixed boxcar window size
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
%   Tip: when N-looks intensity images are used, only BWS can be used@!  
%
%
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen University,  
%
%   ======================================================================
%   11/2021 MJ remove SHP.BroNum-1 (the brother pixel becomes itself when BroNum=1)
%   ======================================================================


if nargin < 4
    EstAgr='FaSHPS';
end

if nargin < 3
    Alpha = 0.05;
end

if nargin < 2
    CalWin = [15 15];
end

if nargin < 1
    help SHP_SelPoint
    return;
end

tic;

if length(size(mlistack))~=3
    error('Please input 3D matrix...');
end

[nlines,nwidths,npages] = size(mlistack);
mlistack=single(mlistack);

%Parameter prepare:
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;
InitRow=(CalWin(1)+1)/2; % InitRow is CenterRow
InitCol=(CalWin(2)+1)/2; % InitCol is CenterCol

%Edge mirror-image
mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');
meanmli = mean(mlistack,3);
[nlines_EP,nwidths_EP]= size(meanmli);
SHP.PixelInd=false(CalWin(1)*CalWin(2),nlines*nwidths);

%estimate SHPs
num=1;
p=1;
all = nlines*nwidths;
all_step = floor(all/10);

if strcmpi(EstAgr,'FaSHPS') 
    
    %Set thresholds for parameter statistics
    LRT_nl = 3;
    LRT_nw = 3;
    if RadiusRow<LRT_nl
        LRT_nl=1;
    end
    if RadiusCol<LRT_nw
        LRT_nw=1;
    end

    %Critical region
    CR_lo = finv(Alpha/2,2*npages,2*npages);
    CR_up = finv(1-Alpha/2,2*npages,2*npages);
    Galpha_L = gaminv(Alpha/2,npages,1);
    Galpha_U = gaminv(1-Alpha/2,npages,1);    
    
    for kk=InitCol:nwidths_EP-RadiusCol
        for ll=InitRow:nlines_EP-RadiusRow       
            %Initial estimation (Likelihood-ratio test)
            temp = meanmli(ll-LRT_nl:ll+LRT_nl,kk-LRT_nw:kk+LRT_nw);
            T = meanmli(ll,kk)./temp;
            T = T>CR_lo&T<CR_up;
            SeedPoint = mean(temp(T));
            %iteration (Gamma Confidence interval)
            MeanMatrix = meanmli(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol);
            SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
            SeedPoint(InitRow,InitCol)=true;
            %connection
            LL = bwlabel(SeedPoint); %double
            SHP.PixelInd(:,num)=LL(:)==LL(InitRow,InitCol);  
            num=num+1;
            if num == all_step * p
                disp(['progress: ', num2str(10*p),'%']);
                p = p+1;
            end
        end
    end
else %BWS (non-parameter statistics)
    for kk=InitCol:nwidths_EP-RadiusCol
        % 展示进度和当前时间
        disp(['Processing column ', num2str(kk-InitCol+1), ' of ', num2str(nwidths), ' at ', datestr(now,'HH:MM:SS')]);
        for ll=InitRow:nlines_EP-RadiusRow
            Matrix = mlistack(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol,:);
            Ref = Matrix(InitRow,InitCol,:);
            T = BWStest(repmat(Ref(:),[1,CalWin(1)*CalWin(2)])...
                ,reshape(Matrix,[CalWin(1)*CalWin(2),npages])',Alpha);   
            temp=reshape(~T,[CalWin(1),CalWin(2)]);
            %connection
            LL = bwlabel(temp);
            SHP.PixelInd(:,num)=LL(:)==LL(InitRow,InitCol);       
            num=num+1;
            if num == all_step * p
                disp(['progress: ', num2str(10*p),'%']);
                p = p+1;
            end              
        end
    end   
end
%SHPs map            
SHP.BroNum = sum(SHP.PixelInd,1);
SHP.BroNum = uint16(reshape(SHP.BroNum(:),[nlines,nwidths]));    
SHP.CalWin = CalWin;            
t=toc;

figure;imagesc(SHP.BroNum);axis image off;colormap jet;
ti=title ('Homogeneous Pixel Number');colorbar;
set(ti,'fontweight','bold');
disp(['SHP_SelPoint operation completed in ',int2str(t/60),' min(s).']);
disp('Done!');

    