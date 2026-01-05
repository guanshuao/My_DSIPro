function ImgWrite(matstack,filename,imgpath,suffixname,bkformat,machinefmt)
%Write data stack to a specified directory
%
%   Usage:
%       ImgWrite(matstack,filename,imgpath,suffixname,bkformat,machinefmt)
%   
%
%   Inputs:
%   - matstack:     A height by width by page matrix
%   - filename:     The file name in which each row corresponds to each page of matstack, e.g., N*1 array with <yyyymmdd> or N*2 array with [<yyyymmdd>,<yyyymmdd>]  
%   - imgpath:      The path of output (including directory name you want)
%   - suffixname:   The suffix of all files for output, e.g., cpxcc 
%   - bkformat:     The format of data for output, e.g., 'float32'
%   - machinefmt:   See fwritebkj.m for details
%
%
%   Examples:
%   To write a batch of intensity series in float32 format with Big-endian 
%   ordering machinefmt to directory named by 'Despeckle', use:
%   ImgWrite(matstack,filename,'/home/user/INSAR/COHEST/Despeckle','mli.sm','float32')
%   where filename =[20100311   
%                    20100424    
%                    20100607 
%                    ...]
%
%   For complex interferogram with GAMMA format, use: 
%   ImgWrite(matstack,filename,'/home/user/INSAR/COHEST/DIFF','diff','cpxfloat32','b')
%   where filename = [20100311    20100424
%                     20100424    20100607
%                     20100607    20100721
%                     20100721    20100903
%                     ...]
%
%
%   Mi JIANG, Sun Yat-sen University, 

if nargin < 6
    machinefmt='b'; %GAMMA software, for example
end

if nargin < 5
    bkformat='float32'; %for *mli,*cc file
end

if nargin < 4
    help ImgWrite
    return;
end

[~,~,npages]=size(matstack);
if npages~=size(filename,1)
    error('The length filename is not equivalent to the number of image.')
end

if ~isempty(strmatch(imgpath(end),filesep))
    imgpath=imgpath(1:end-1);
end

if ~exist(imgpath,'dir')
    mkdir(imgpath);
else
    k = strfind(imgpath,filesep);
    disp([imgpath(k(end)+1:end), ' directory already exists...']);
end    

disp(['The number of the output ', suffixname,' images:' num2str(npages)]);
for ii=1:npages
    tic;
    if size(filename,2)~=1
        file_name = [num2str(filename(ii,1)),'_',num2str(filename(ii,2)),'.',suffixname];
    else
        file_name = [num2str(filename(ii)),'.',suffixname];
    end
    fwritebkj(matstack(:,:,ii), [imgpath,filesep,file_name], bkformat,machinefmt);
    time=toc;
    fprintf('Writing Img %3d / %d, time = %.0f sec\n',ii,npages,time);      
end