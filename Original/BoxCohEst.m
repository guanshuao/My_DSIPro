function [Coh,Ph] = BoxCohEst(mlistack,mlilist,infstack,inflist,CalWin)
%This function is used to estimate coherence and phase by Boxcar filtering
%supports N-looks images.
%   Usage:
%       [Coh,Ph] = BoxCohEst(mlistack,mlilist,infstack,inflist,CalWin)
%
%
%   Inputs:
%   - mlistack: A height by width by page (real) matrix,e.g., SAR single-look
%               intensity series
%   - mlilist:  A n*1 file list including the intensity images named by <yyyymmdd>
%   - infstack: A height by width by page (complex) matrix,e.g., InSAR
%               single-look interferogram after removing phase gradient (topo., disp.)
%   - inflist:  A n*2 file list including complex interferogram named by <yyyymmdd yyyymmdd>
%   - CalWin:   The window used to collect the samples for coherence
%               estimation [az rg] 
%
%   Outputs:
%   - Coh:      Coherence magnitude stack (real)
%   - Ph:       Interferometric phase (real)           
%
% 
%
% 
% 
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen University,  

%   ======================================================================
%   11/2021 MJ replace Gaussian lowpass filter with average filter
%   ======================================================================

if nargin < 5
   CalWin= [7 7]; %[az, rg]
end

if nargin < 4
    help BoxCohEst
    return
end

tic;
[nlines,nwidths,npages]=size(infstack);
Coh=zeros(nlines,nwidths,npages,'single');
[~,idx]=ismember(inflist,mlilist);
h = fspecial('average',CalWin); % N-looks data can only use 'average' 
CohEstAgr='fast';
if strcmpi(CohEstAgr,'fast')
    for ii=1:npages
        % intensity image of m and s
        m1 = mlistack(:,:,idx(ii,1));
        m2 = mlistack(:,:,idx(ii,2));
        % numerator of coherence estimator
        nu = filter2(h,sqrt(m1.*m2).*infstack(:,:,ii));
        de1 = filter2(h,m1);
        de2 = filter2(h,m2);        
        Coh(:,:,ii) = nu./sqrt(de1.*de2);
    end   
else %Note: conventional programming 
    
    RadiusRow=(CalWin(1)-1)/2;
    RadiusCol=(CalWin(2)-1)/2;

    mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');    
    infstack = padarray(infstack,[RadiusRow RadiusCol],'symmetric');   
    
    for ii=1:npages
        % intensity image of m and s
        m1 = mlistack(:,:,idx(ii,1));
        m2 = mlistack(:,:,idx(ii,2));
        % numerator of coherence estimator
        nu = sqrt(m1.*m2).*infstack(:,:,ii);
        for jj=1:nwidths
            for kk=1:nlines                
                x_patch = jj+RadiusCol;
                y_patch = kk+RadiusRow;        
                m1_tmp  = m1(y_patch-RadiusRow:y_patch+RadiusRow,x_patch-RadiusCol:x_patch+RadiusCol); 
                m2_tmp  = m2(y_patch-RadiusRow:y_patch+RadiusRow,x_patch-RadiusCol:x_patch+RadiusCol);
                nu_tmp  = nu(y_patch-RadiusRow:y_patch+RadiusRow,x_patch-RadiusCol:x_patch+RadiusCol);
                m1_tmp  = mean(m1_tmp(:));
                m2_tmp  = mean(m2_tmp(:));
                nu_tmp  = mean(nu_tmp(:));
                Coh(kk,jj,ii) = nu_tmp./sqrt(m1_tmp.*m2_tmp);
            end
        end
    end   
end
    

if nargout>1 %filtering mli    
    Ph=angle(Coh);
end

Coh = abs(Coh);

t=toc;
disp(['BoxCohEst operation completed in ',num2str(t/60),' min(s).']);
disp('Done!');    
