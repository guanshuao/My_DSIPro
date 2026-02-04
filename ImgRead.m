function Data=ImgRead(imgpath,suffixname,nline,bkformat,machinefmt,Blk)
%ImgRead 从指定目录读取数据堆栈。
%   用法:
%       Data=ImgRead(imgpath,suffixname,nline,bkformat,machinefmt);
%   
%
%   输入:
%   - imgpath:      图像集的路径 (字符串)
%   - suffixname:   图像集中所有文件的后缀，例如 'mli' (字符串)
%   - nline:        图像高度/行数 (数值)
%   - bkformat:     数据格式，详见 freadbkj.m (例如 'float32', 'cpxfloat32')
%   - machinefmt:   字节序，详见 freadbkj.m (例如 'b' 表示大端)
%   - Blk:          读取块的范围 [r0, rN, c0, cN] (起始行, 结束行, 起始列, 结束列)
%
%   输出:
%   - Data.datastack:   一个 高度 x 宽度 x 页数 的矩阵，每一页对应一个 2D 图像
%   - Data.filename:    文件名列表 (通常包含从文件名提取的日期)
%
%   示例:
%   读取一批 float32 格式、大端字节序、高度为 200 行的强度序列:
%   Data=ImgRead('/home/user/INSAR/COHEST/MLI','mli',200,'float32');
%
%   读取高度为 1800 的复数差分干涉图:
%   Data=ImgRead('/home/user/INSAR/COHEST/DIFF','diff',1800,'cpxfloat32');
%
%
%   Mi JIANG, Sun Yat-sen University, 
%
%   ======================================================================
%   012/2017 MJ add Blk
%   ======================================================================
if nargin < 6 || isempty(Blk)
    Blk =[]; % 如果未提供 Blk 参数，则初始化为空
end

if nargin < 5 || isempty(machinefmt)
    machinefmt='b'; % 默认为大端字节序 (例如 GAMMA 软件)
end

if nargin < 4 || isempty(bkformat)
    bkformat='float32'; % 默认格式为 float32 (适用于 *mli, *cc 文件)
end

if nargin < 3 || isempty(nline)
    help ImgRead % 如果参数不足，显示帮助信息
    return;
end

% 确保路径末尾有分隔符
if ~strcmp(imgpath(end),filesep)
    imgpath=[imgpath,filesep];
end

% 获取指定路径下所有匹配后缀的文件
tag_files = dir([imgpath,'*',suffixname]);
img_num = length(tag_files);
disp(['The number of the ', suffixname,' images:' num2str(img_num)]); % 显示找到的图像数量

for ii=1:img_num
    tic; % 开始计时
    if isempty(Blk)
        % 读取整个文件
        Data.datastack(:,:,ii)=single(freadbkj([imgpath,tag_files(ii).name],nline,bkformat,machinefmt));
    else
        % 读取指定块
        Data.datastack(:,:,ii)=single(freadbkj([imgpath,tag_files(ii).name],nline,bkformat,machinefmt,Blk(1),Blk(2),Blk(3),Blk(4)));
    end
    
    % 从文件名中提取数字（通常是日期）
    temp=regexp(tag_files(ii).name,'\d+','match');
    if isscalar(temp)      % 单个日期，例如 mli 文件 (yyyymmdd)
        Data.filename(ii,1)=str2double(temp{1});
    elseif length(temp)==2 % 两个日期，例如干涉图 intf (yyyymmdd_yyyymmdd)
        Data.filename(ii,1)=str2double(temp{1});
        Data.filename(ii,2)=str2double(temp{2});
    else
        error('The format of file name should be: <yyyymmdd> or <yyyymmdd_yyyymmdd>. (文件名格式错误)')
    end
    time=toc; % 结束计时
    fprintf('Reading Img %3d / %d, time = %.0f sec\n',ii,img_num,time);      
end