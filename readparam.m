function value = readparam(param,fname)
%readparam 从参数文件中读取特定参数，例如 GAMMA 软件的 par 文件
%   输入:
%       param: 要查找的参数名称 (字符串)
%       fname: 参数文件的路径 (字符串)
%   输出:
%       value: 参数对应的值 (数值)
%
%   作者: Mi JIANG, Sun Yat-sen University

if nargin < 1
    help readparam
    return;
end

% 打开文件
fileID=fopen(fname);
% 读取文件内容为字符串单元数组
params=textscan(fileID,'%s');
% 关闭文件
fclose(fileID);
% 获取单元数组中的内容
params=params{1};

% 查找参数名称在数组中的索引
% 注意: strmatch 在较新版本的 MATLAB 中已不推荐使用，建议改用 strncmp 或 validatestring，但在旧代码中常见
idx = strmatch(param,params);

% 如果未找到参数，抛出错误
if isempty(idx)
    error('empty value ...');
end

% 读取参数名称后面的一个字符串，并将其转换为数值
% 假设参数文件格式为: 参数名 参数值
value = str2double(params(idx+1));

% 如果转换结果为 NaN (非数值)，说明值不存在或格式不正确
if isnan(value)
    error('the value does not exist...');
end

