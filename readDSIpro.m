function params = readDSIpro(fname)
%This function is used to read parameters from param_list file into
%workspace
%   Usage:
%       []=readDSIpro('C:\Users\DELL\Desktop\DSIpro\data_list');
%
%   run parmslist_generator first to generate param_list file, modify the
%   parameters desired, and then read param_list by readDSIpro.
%
%the details of param_list are given as follows,
%   lines:     		lines of reference image
%   widths:    		widths of reference image
%   rlks:      		looks in range 
%   azlks:     		looks in azimuth
%   rg_spacing:		range pixel spacing [m]
%   az_spacing:		azimuth pixel spacing [m]
%   inc:       		incidence angle [degree]
% 
%   %source path
%   proc_path:		directory under process
%   diff_path:		directory including single-master differential interferogram series (complex)
%   mli_path: 		directory including intensity SAR series 
%   pl_path:  		directory of optimal differential interferogram series (complex)
% 
%   %processing parameters
%   pl_bro:         threshold of number of brother pixels [pixels]
%   pl_pcoh:        threshold of posteriori coherence
%   pl_blksize:     block size to process big image series
%   shp_calwin:	 	sliding window size [pixels]
%   shp_method:     method for homogeneous pixel selection
%   shp_alpha:      significance level of hypothesis test
%   apsp_pixsize: 	resampling interval for spatial network construction [m]
%   apsp_dist:		maximum distance between two spatial pixels [m]
%
%   Mi JIANG, Sun Yat-sen University,
if nargin < 1
    help readDSIpro
    return;
end

if ~exist(fname, 'file')
    error('no param_list can be found!');
end

fileID=fopen(fname);
tmp=textscan(fileID,'%s');
fclose(fileID);
tmp=tmp{1};

%line
idx = strmatch('lines',tmp);
params.nlines = str2double(tmp(idx+1));

%width
idx = strmatch('widths',tmp);
params.nwidths = str2double(tmp(idx+1));

%range looks
idx = strmatch('rlks',tmp);
params.rlks = str2double(tmp(idx+1));

%azimuth looks
idx = strmatch('azlks',tmp);
params.azlks = str2double(tmp(idx+1));

%range pixel spaceing 
idx = strmatch('rg_spacing',tmp);
params.rg_spacing = str2double(tmp(idx+1));

%azimuth pixel spaceing 
idx = strmatch('az_spacing',tmp);
params.az_spacing = str2double(tmp(idx+1));

%incidence angle
idx = strmatch('inc',tmp);
params.inc = str2double(tmp(idx+1));

%diff path
idx = strmatch('diff_path',tmp);
params.diff_path = tmp{idx+1};

%mli path
idx = strmatch('mli_path',tmp);
params.mli_path = tmp{idx+1};

%pl path
idx = strmatch('pl_path',tmp);
params.pl_path = tmp{idx+1};

%threshold of number of brother pixels
idx = strmatch('pl_bro',tmp);
params.pl_bro = str2double(tmp(idx+1));

%threshold of posteriori coherence
idx = strmatch('pl_pcoh',tmp);
params.pl_pcoh = str2double(tmp(idx+1));

%block size to process big image series
idx = strmatch('pl_blksize',tmp);
params.pl_blksize = str2double(tmp(idx+1));

%sliding window size
idx = strmatch('shp_calwin',tmp);
params.shp_calwin = [str2double(tmp(idx+1)),str2double(tmp(idx+2))];

%sliding window size
idx = strmatch('shp_method',tmp);
params.shp_method = tmp{idx+1};

%sliding window size
idx = strmatch('shp_alpha',tmp);
params.shp_alpha = str2double(tmp(idx+1));

%maximum distance window
idx = strmatch('apsp_dist',tmp);
params.apsp_dist = str2double(tmp(idx+1));

%resample size
idx = strmatch('apsp_pixsize',tmp);
params.apsp_pixsize = str2double(tmp(idx+1));