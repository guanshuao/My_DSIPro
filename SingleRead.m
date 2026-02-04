function img = SingleRead(filepath,nline,bkformat,machinefmt,Blk)
%SingleRead 读取单张二进制图像文件。
%   基于 ImgRead 改写，ImgRead 读取整个数据堆栈，本函数仅读取单张图像。
%
%   用法:
%       img = SingleRead(filepath,nline,bkformat,machinefmt,Blk)
%
%   输入:
%   - filepath:   图像文件的完整路径 (字符串)
%   - nline:      图像高度/行数 (数值)
%   - bkformat:   数据格式，详见 freadbkj.m (例如 'float32', 'cpxfloat32')
%   - machinefmt: 字节序，详见 freadbkj.m (例如 'b' 表示大端，与 ImgRead 默认值一致)
%   - Blk:        可选，读取块的范围 [r0, rN, c0, cN] (起始行, 结束行, 起始列, 结束列)
%
%   输出:
%   - img:        二维矩阵，存储读取的单张图像数据
%
%   示例:
%   img = SingleRead('/path/to/MLI/20100311.mli',200,'float32','b');
%   img = SingleRead('D:\Research\Coh_Est\SL\MLI\20220506.rslc',1201,'cpxfloat32');
%
%   参见: ImgRead, freadbkj

%% 参数默认值设置
if nargin < 5 || isempty(Blk)
    Blk = []; % 如果未提供 Blk 参数，则初始化为空（读取整张图像）
end

if nargin < 4 || isempty(machinefmt)
    machinefmt = 'b'; % 默认为大端字节序 (例如 GAMMA 软件)
end

if nargin < 3
    help SingleRead % 如果参数不足，显示帮助信息
    return;
end

%% 输入参数校验
% 检查 filepath 是否为有效的字符串类型
if ~ischar(filepath) && ~isstring(filepath)
    error('filepath must be a character vector or string scalar.');
end
filepath = char(filepath); % 统一转换为字符数组格式

% 检查文件是否存在
if ~exist(filepath,'file')
    error('File does not exist: %s', filepath);
end

%% 读取图像数据
tic; % 开始计时
if isempty(Blk)
    % 读取整张图像
    img = single(freadbkj(filepath,nline,bkformat,machinefmt));
else
    % 读取指定块区域
    img = single(freadbkj(filepath,nline,bkformat,machinefmt,Blk(1),Blk(2),Blk(3),Blk(4)));
end
elapsed = toc; % 结束计时

% 从文件路径中提取文件名并显示读取信息
[~,name,~] = fileparts(filepath);
fprintf('Reading Img %s, time = %.0f sec\n', name, elapsed);

end
