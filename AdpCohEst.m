function [Coh,Ph] = AdpCohEst(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode)
%This function is used to estimate coherence and interferometric phase using SHP,
%supports N-looks images.
%
%   Usage:
%       [Coh,Ph] = AdpCohEst(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode)
%   or:
%       [Coh] = AdpCohEst(mlistack,mlilist,infstack,inflist,SHP,BiasCorr,OutputMode)
%
%
%   Inputs:
%   - mlistack: A height by width by page (real) matrix,e.g., SAR single-look
%               intensity series
%   - mlilist:  A n*1 file list including the intensity images named by <yyyymmdd>
%   - infstack: A height by width by page (complex) matrix,e.g., InSAR
%               single-look interferogram after removing phase gradient (topo., disp.)
%   - inflist:  A n*2 file list including complex interferogram named by <yyyymmdd yyyymmdd>
%   - SHP:      See script "SHP_SelPoint.m" for details
%   - BiasCorr: Open log-moment estimate to mitigate both bias and variance of
%               coherence observations 
%   - OutputMode: Output mode, 'stack' or 'average'
%
%   Outputs:
%   - Coh:      Coherence magnitude stack (real)
%   - Ph:       Interferometric phase (real)         
%
%   [1] Delineation of built-up land change from SAR stack by analysing the coefficient of variation
%        M Jiang, A Hooper, X Tian, J Xu, SN Chen, ZF Ma, X Cheng,
%        ISPRS Journal of Photogrammetry and Remote Sensing 169, Nov. 2020,
%        pp. 93-108.
%
%   [2] Hybrid Approach for Unbiased Coherence Estimation for Multitemporal InSAR
%        Mi Jiang, Xiaoli Ding and Zhiwei Li,
%        IEEE Transactions on Geoscience and Remote Sensing vol. 52, no. 5,
%        May  2014, pp. 2459-2473.
% 
%
% 
% 
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen University,  

%   =====================================================
%   MJ 11/2021 remove boxcar, bootstrapping and Acc
%   MJ 11/2021 add bias correction opinion and phase output
%   =====================================================
if nargin < 7 || isempty(OutputMode) % 如果未提供 OutputMode 参数，则使用默认值 'stack'
   OutputMode = 'stack'; % 'stack' or 'average'
end

if nargin < 6 || isempty(BiasCorr) % 如果未提供 BiasCorr 参数，则使用默认值 'n'
   BiasCorr = 'n'; % 'y' or 'n'
end

if nargin < 5
    help AdpCohEst
    return
end

tic;
[nlines,nwidths,npages]=size(infstack);
[~,idx]=ismember(inflist,mlilist);

CalWin =SHP.CalWin;
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;   

if strcmpi(OutputMode,'stack')
    % stack 模式：保存完整堆栈
    Coh=zeros(nlines,nwidths,npages,'single');
    
    for ii=1:npages
        m1 = mlistack(:,:,idx(ii,1));
        m2 = mlistack(:,:,idx(ii,2));
        Intf= sqrt(m1.*m2).*infstack(:,:,ii);

        m1 = padarray(m1,[RadiusRow RadiusCol],'symmetric');
        m2 = padarray(m2,[RadiusRow RadiusCol],'symmetric');
        Intf= padarray(Intf,[RadiusRow RadiusCol],'symmetric');           
        nu = zeros(nlines,nwidths,'single');
        de1=nu;
        de2=nu;
        num=1;
        for jj = 1:nwidths
            for kk= 1:nlines
                x_global  = jj+RadiusCol;
                y_global  = kk+RadiusRow;
                MasterValue= m1(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                SlaveValue = m2(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                InterfValue= Intf(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                MasterValue= MasterValue(SHP.PixelInd(:,num));
                SlaveValue = SlaveValue(SHP.PixelInd(:,num));
                InterfValue= InterfValue(SHP.PixelInd(:,num));
                nu(kk,jj)  = mean(InterfValue, 'omitnan');
                de1(kk,jj) = mean(MasterValue, 'omitnan');
                de2(kk,jj) = mean(SlaveValue, 'omitnan');       
                num=num+1;
            end
        end
        denominator = sqrt(de1.*de2);
        denominator(denominator == 0) = NaN;
        Coh(:,:,ii) = nu./denominator;   
        fprintf('adp. coherence: %3d / %d is finished...\n',ii,npages);
    end

    % 处理NaN和Inf值
    Coh(isnan(Coh) | isinf(Coh)) = 0; 
    
    % 先保存复数相干性用于后续相位计算
    CohComplex = Coh;
    
    % 取相干性幅度
    Coh = abs(Coh);

    % Bias mitigation for coherence magnitude (stack mode)
    if strcmpi(BiasCorr,'y')
        tmp = padarray(Coh,[RadiusRow RadiusCol],'symmetric');
        for ii=1:npages 
            num  = 1;        
            for jj = 1:nwidths
                for kk= 1:nlines
                    x_global  = jj+RadiusCol;
                    y_global  = kk+RadiusRow;
                    CohValue  = tmp(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol,ii);
                    CohValue  = CohValue(SHP.PixelInd(:,num));
                    CohValue(CohValue <= 0) = NaN;
                    Coh(kk,jj,ii) = mean(log(CohValue), 'omitnan');
                    num=num+1;
                end
            end 
            fprintf('bias mitigation: %3d / %d is finished...\n',ii,npages);  
        end
        Coh = exp(Coh);
        Coh(isnan(Coh) | isinf(Coh)) = 0;
    end
    
    if nargout > 1
        Ph = angle(CohComplex);
    end
end

if strcmpi(OutputMode,'average')
    % average 模式：使用中间变量累加，节省内存
    CohSum = zeros(nlines,nwidths,'single');
    
    for ii=1:npages
        m1 = mlistack(:,:,idx(ii,1));
        m2 = mlistack(:,:,idx(ii,2));
        Intf= sqrt(m1.*m2).*infstack(:,:,ii);

        m1 = padarray(m1,[RadiusRow RadiusCol],'symmetric');
        m2 = padarray(m2,[RadiusRow RadiusCol],'symmetric');
        Intf= padarray(Intf,[RadiusRow RadiusCol],'symmetric');           
        nu = zeros(nlines,nwidths,'single');
        de1=nu;
        de2=nu;
        num=1;
        for jj = 1:nwidths
            for kk= 1:nlines
                x_global  = jj+RadiusCol;
                y_global  = kk+RadiusRow;
                MasterValue= m1(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                SlaveValue = m2(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                InterfValue= Intf(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                MasterValue= MasterValue(SHP.PixelInd(:,num));
                SlaveValue = SlaveValue(SHP.PixelInd(:,num));
                InterfValue= InterfValue(SHP.PixelInd(:,num));
                nu(kk,jj)  = mean(InterfValue, 'omitnan');
                de1(kk,jj) = mean(MasterValue, 'omitnan');
                de2(kk,jj) = mean(SlaveValue, 'omitnan');       
                num=num+1;
            end
        end
        denominator = sqrt(de1.*de2);
        denominator(denominator == 0) = NaN;
        CohSingle = nu./denominator;
        
        % 处理NaN和Inf值
        CohSingle(isnan(CohSingle) | isinf(CohSingle)) = 0;
        
        % 取相干性幅度
        CohSingle = abs(CohSingle);
        
        % Bias mitigation for current coherence image (average mode)
        if strcmpi(BiasCorr,'y')
            tmp = padarray(CohSingle,[RadiusRow RadiusCol],'symmetric');
            CohCorrected = zeros(nlines,nwidths,'single');
            num = 1;
            for jj = 1:nwidths
                for kk= 1:nlines
                    x_global  = jj+RadiusCol;
                    y_global  = kk+RadiusRow;
                    CohValue  = tmp(y_global-RadiusRow:y_global+RadiusRow,x_global-RadiusCol:x_global+RadiusCol);
                    CohValue  = CohValue(SHP.PixelInd(:,num));
                    CohValue(CohValue <= 0) = NaN;
                    CohCorrected(kk,jj) = mean(log(CohValue), 'omitnan');
                    num=num+1;
                end
            end
            CohSingle = exp(CohCorrected);
            CohSingle(isnan(CohSingle) | isinf(CohSingle)) = 0;
        end
        
        % 累加到中间变量
        CohSum = CohSum + CohSingle;
        fprintf('adp. coherence: %3d / %d is finished...\n',ii,npages);
    end
    
    % 计算平均值
    Coh = CohSum / npages;
    
    if nargout > 1
        warning('average 模式下只输出 Coh，Ph 输出将被忽略。');
    end
end

t=toc;
disp(['AdpCohEst operation completed in ',int2str(t/60),' min(s).']);
disp('Done!');
