function Mliimg = DeSpeckling(mlistack,SHP)
%This function filters the intensity stack image, returns boxcar filtered intensity stack image without SHP file
%supports N-looks inensity images.
%
%   Usage:
%      Mliimg = DeSpeckling(mlistack,SHP)
%   
%
%   Inputs:
%   - mlistack: A height by width by page (real) matrix,e.g., SAR single-look
%               intensity series
%   - SHP:      See script "SHP_SelPoint.m" for details
%
%
%   Outputs:
%   - Mliimg:   filtered intensity images 
%
%
%   [1] Fast Statistically Homogeneous Pixel Selection for Covariance Matrix Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding, Ramon F. Hanssen, Rakesh Malhotra and Ling Chang,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 53, no. 3,
%        March 2015, pp. 1213-1224.
% 
%
% 
%
% 
% 
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen University,  

%   ===========================================
%   MJ 11/2021 remove Acc and Gaussian kernel
%   ===========================================

BoxFilt=false;


if nargin < 2
    BoxFilt=true;
end

if nargin < 1
    help DeSpeckling
    return
end

tic;
[nlines,nwidths,npages]=size(mlistack);
Mliimg = mlistack;
if BoxFilt
    h = fspecial('average',[7 7]); 
    for ii=1:npages
        Mliimg(:,:,ii) = filter2(h,mlistack(:,:,ii));
        fprintf('BOXCAR DESPECKLING: %d / %d is finished...\n',ii,npages);
    end     
else
    CalWin =SHP.CalWin;
    RadiusRow=(CalWin(1)-1)/2;
    RadiusCol=(CalWin(2)-1)/2;  
    mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');
    %Despeckling    
    for ii=1:npages
        temp = mlistack(:,:,ii);
        num=1;
        for jj = 1:nwidths
            for kk= 1:nlines
                x_global  = jj+RadiusCol;
                y_global  = kk+RadiusRow;
                MliValue  = temp(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                MliValue  = MliValue(SHP.PixelInd(:,num));
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
        
   
