function []=parmslist_generator(diff_path,mli_path,Imgparams,shp_calwin,shp_alpha,shp_method,pl_pcoh,pl_bro,pl_blksize,apsp_dist,apsp_pixsize)
%parmslist_generator 用于生成 DSIpro 的参数列表文件
%   用法:
%       []=parmslist_generator('C:\Users\DELL\Desktop\DSIpro\SL\DIFF','C:\Users\DELL\Desktop\DSIpro\SL\MLI');
%       []=parmslist_generator;
%   运行后，你可以在当前文件夹中检查 'param_list' 文件，并手动填充参数，然后运行 readDSIpro 将参数读取到工作区中。
%
%   输入参数 (可选):
%   - diff_path:    差分干涉图路径
%   - mli_path:     强度图路径
%   - Imgparams:    图像参数结构体 (包含 lines, widths, rlks, azlks, spacing, inc 等)
%   - shp_calwin:   SHP 选择窗口大小 [方位向 距离向] (默认 [15 15])
%   - shp_alpha:    SHP 选择显著性水平 (默认 0.05)
%   - shp_method:   SHP 选择方法 (默认 'FaSHPS')
%   - pl_pcoh:      Phase Linking 后验相干性阈值 (默认 0.75)
%   - pl_bro:       SHP 数量阈值 (默认 20)
%   - pl_blksize:   Phase Linking 处理块大小 (默认 3)
%   - apsp_dist:    APSP (全对最短路径) 距离阈值 (默认 1000 m)
%   - apsp_pixsize: APSP 重采样像素大小 (默认 100 m)
%
%   Mi JIANG, Sun Yat-sen University,


if nargin < 11 || isempty(apsp_pixsize)
    apsp_pixsize = 100; %m (APSP 重采样网格大小)
end

if nargin < 10 || isempty(apsp_dist)
    apsp_dist = 1000; %m (APSP 搜索距离阈值)
end

if nargin < 9 || isempty(pl_blksize)
    pl_blksize = 3; % (Phase Linking 块大小)
end

if nargin < 8 || isempty(pl_bro)
    pl_bro = 20; %pixels (SHP 数量阈值)
end

if nargin < 7 || isempty(pl_pcoh)
    pl_pcoh = .75; % (后验相干性阈值)
end

if nargin < 6 || isempty(shp_method)
    shp_method='FaSHPS'; % (SHP 选择算法)
end

if nargin < 5 || isempty(shp_alpha)
    shp_alpha=0.05; % (显著性水平)
end

if nargin < 4 || isempty(shp_calwin)
    shp_calwin = [15 15];%[az rg] (SHP 窗口大小)
end

if nargin < 3 || isempty(Imgparams)
    Imgparams =[]; % (图像参数结构体)
end

if nargin < 2 || isempty(mli_path)
    mli_path=' '; % (强度图路径)
end

if nargin < 1 || isempty(diff_path)
    diff_path=' '; % (干涉图路径)
end

% 解析图像参数
if isempty(Imgparams)
    lines  = 0;
    widths = 0;
    rlks   = 0;
    azlks  = 0;
    rg_spacing = 0;
    az_spacing = 0;
    inc    = 0;
else
    lines  = Imgparams.nlines;
    widths = Imgparams.nwidths;
    rlks   = Imgparams.rlks;
    azlks  = Imgparams.azlks;
    rg_spacing = Imgparams.rg_spacing;
    az_spacing = Imgparams.az_spacing;
    inc    = Imgparams.inc;   
end


parmslistname = 'param_list';
fid=fopen(parmslistname,'w');

% 写入文件头
fprintf(fid,'Distributed scatterer interferometry processor (DSIpro) - Parameter File\n\n');

fprintf(fid,'date:\t%s\n',date);

% 写入图像参数
fprintf(fid,'%%image parameters\n');
fprintf(fid,'lines:        \t\t%d\n',lines);
fprintf(fid,'widths:       \t\t%d\n',widths);
fprintf(fid,'rlks:         \t\t%d\n',rlks);
fprintf(fid,'azlks:        \t\t%d\n',azlks);
fprintf(fid,'rg_spacing:   \t\t%.6f\n',rg_spacing);
fprintf(fid,'az_spacing:   \t\t%.6f\n',az_spacing);
fprintf(fid,'inc:          \t\t%.6f\n\n',inc);

% 写入路径信息
fprintf(fid,'%%source path\n');
fprintf(fid,'diff_path:    \t\t%s\n',diff_path);
fprintf(fid,'mli_path:     \t\t%s\n',mli_path);
fprintf(fid,'pl_path:      \t\t%s\n\n',[pwd,filesep,'PL']);

% 写入处理参数
fprintf(fid,'%%processing parameters\n');
fprintf(fid,'pl_bro:       \t\t%d\n',pl_bro);
fprintf(fid,'pl_pcoh:      \t\t%.2f\n',pl_pcoh);
fprintf(fid,'pl_blksize:   \t\t%d\n',pl_blksize);
fprintf(fid,'shp_calwin:   \t\t%d %d\n',shp_calwin(1),shp_calwin(2));
fprintf(fid,'shp_method:   \t\t%s\n',shp_method);
fprintf(fid,'shp_alpha:    \t\t%.2f\n',shp_alpha);
fprintf(fid,'apsp_pixsize: \t\t%d\n',apsp_pixsize);
fprintf(fid,'apsp_dist:    \t\t%d\n',apsp_dist);

fclose(fid);