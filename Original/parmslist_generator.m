function []=parmslist_generator(diff_path,mli_path,Imgparams,shp_calwin,shp_alpha,shp_method,pl_pcoh,pl_bro,pl_blksize,apsp_dist,apsp_pixsize)
%This function is used to make parameter list for DSIpro
%   Usage:
%       []=parmslist_generator('C:\Users\DELL\Desktop\DSIpro\SL\DIFF','C:\Users\DELL\Desktop\DSIpro\SL\MLI');
%       []=parmslist_generator;
%   you can check param_list in current folder, and fill the parameters manually,then
%   run readDSIpro to read parameters into workspace
%
%   Mi JIANG, Sun Yat-sen University,


if nargin < 11
    apsp_pixsize = 100; %m
end

if nargin < 10
    apsp_dist = 1000; %m
end

if nargin < 9
    pl_blksize = 3;
end

if nargin < 8
    pl_bro = 20; %pixels
end

if nargin < 7
    pl_pcoh = .75;
end

if nargin < 6
    shp_method='FaSHPS';
end

if nargin < 5
    shp_alpha=0.05;
end

if nargin < 4
    shp_calwin = [15 15];%[az rg]
end

if nargin < 3
    Imgparams =[];
end

if nargin < 2
    mli_path=' ';
end

if nargin < 1
    diff_path=' ';
end

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

fprintf(fid,'Distributed scatterer interferometry processor (DSIpro) - Parameter File\n\n');

fprintf(fid,'date:\t%s\n',date);

fprintf(fid,'%%image parameters\n');
fprintf(fid,'lines:        \t\t%d\n',lines);
fprintf(fid,'widths:       \t\t%d\n',widths);
fprintf(fid,'rlks:         \t\t%d\n',rlks);
fprintf(fid,'azlks:        \t\t%d\n',azlks);
fprintf(fid,'rg_spacing:   \t\t%.6f\n',rg_spacing);
fprintf(fid,'az_spacing:   \t\t%.6f\n',az_spacing);
fprintf(fid,'inc:          \t\t%.6f\n\n',inc);

fprintf(fid,'%%source path\n');
fprintf(fid,'diff_path:    \t\t%s\n',diff_path);
fprintf(fid,'mli_path:     \t\t%s\n',mli_path);
fprintf(fid,'pl_path:      \t\t%s\n\n',[pwd,filesep,'PL']);

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