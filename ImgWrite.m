function ImgWrite(matstack,filename,imgpath,suffixname,bkformat,machinefmt)
%ImgWrite 将数据堆栈写入指定目录
%
%   用法:
%       ImgWrite(matstack,filename,imgpath,suffixname,bkformat,machinefmt)
%   
%
%   输入:
%   - matstack:     高度 x 宽度 x 页数 的矩阵 (要写入的数据)
%   - filename:     文件名列表，每一行对应 matstack 的一页。
%                   例如: N*1 数组 <yyyymmdd> (强度图) 或 N*2 数组 [<yyyymmdd>,<yyyymmdd>] (干涉图)
%   - imgpath:      输出路径 (包含你想要的目录名)
%   - suffixname:   输出文件的后缀名，例如 'cpxcc'
%   - bkformat:     输出数据的格式，例如 'float32' (详见 fwritebkj.m)
%   - machinefmt:   字节序，例如 'b' (大端) (详见 fwritebkj.m)
%
%
%   示例:
%   将一批 float32 格式、大端字节序的强度序列写入 'Despeckle' 目录:
%   ImgWrite(matstack,filename,'/home/user/INSAR/COHEST/Despeckle','mli.sm','float32')
%   其中 filename =[20100311   
%                    20100424    
%                    20100607 
%                    ...]
%
%   对于 GAMMA 格式的复数干涉图: 
%   ImgWrite(matstack,filename,'/home/user/INSAR/COHEST/DIFF','diff','cpxfloat32','b')
%   其中 filename = [20100311    20100424
%                     20100424    20100607
%                     20100607    20100721
%                     20100721    20100903
%                     ...]
%
%
%   Mi JIANG, Sun Yat-sen University, 

if nargin < 6 || isempty(machinefmt)
    machinefmt='b'; % 默认为大端字节序 (例如 GAMMA 软件)
end

if nargin < 5 || isempty(bkformat)
    bkformat='float32'; % 默认格式为 float32 (适用于 *mli, *cc 文件)
end

if nargin < 4
    help ImgWrite
    return;
end

[~,~,npages]=size(matstack);

% 检查文件名列表长度是否与图像页数一致
if npages~=size(filename,1)
    error('The length filename is not equivalent to the number of image.')
end

% 移除路径末尾的分隔符 (如果存在)
if ~isempty(strmatch(imgpath(end),filesep))
    imgpath=imgpath(1:end-1);
end

% 如果输出目录不存在，则创建它
if ~exist(imgpath,'dir')
    mkdir(imgpath);
else
    k = strfind(imgpath,filesep);
    disp([imgpath(k(end)+1:end), ' directory already exists...']);
end    

disp(['The number of the output ', suffixname,' images:' num2str(npages)]);
for ii=1:npages
    tic;
    % 根据 filename 的列数构造输出文件名
    if size(filename,2)~=1
        % 双日期格式 (例如干涉图: 20100311_20100424.diff)
        file_name = [num2str(filename(ii,1)),'_',num2str(filename(ii,2)),'.',suffixname];
    else
        % 单日期格式 (例如强度图: 20100311.mli)
        file_name = [num2str(filename(ii)),'.',suffixname];
    end
    % 写入文件
    fwritebkj(matstack(:,:,ii), [imgpath,filesep,file_name], bkformat,machinefmt);
    time=toc;
    fprintf('Writing Img %3d / %d, time = %.0f sec\n',ii,npages,time);      
end