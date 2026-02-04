clc;clear;close all
%This script is used to estimate the optimal phase by phase linking 
%algorithm and the spatial network using shortest path algorithm.
%
%Single-master single look time-series has a big data volume in general, which takes 
%~2GB in workspace even for 1000*2500*20 pixels in example A.To reduce 
%memory burden for larger volume, the input and output are changed to param_list
%After generate a param_list by parmslist_generator, one can set parameters
%manually, or using the script below


%--------MAKE INPUT FILES FOR 'data_list' and 'param_list'---------

%single-master time-series data under single look
curpath   = pwd;
mlipath   = [curpath,filesep,'SL',filesep,'MLI'];
diffpath  = [curpath,filesep,'SL',filesep,'DIFF'];

%read parameter file
tag_files = dir([mlipath,filesep, '*.par']);
fname     = [mlipath,filesep,tag_files(1).name]; 
ImgParms  = readGAMMAparams(fname);

%make param_list
parmslist_generator(diffpath,mlipath,ImgParms);

%make data list
datalistname = 'data_list';
if ~exist(datalistname,'file')
    tag_files = dir([diffpath,filesep,'*.diff']);
    fid       = fopen(datalistname,'w');  
    for ii =1:length(tag_files)
        temp  = regexp(tag_files(ii).name,'\d+','match');
        fprintf(fid,'%s %s\n',temp{1},temp{2}); 
    end
    fclose(fid);
end
datalist = load(datalistname);

%alternative: establish 'data_list' and 'param_list' manually
%data_list: make a file with 2 channels (master secondary), i.e,
% 20170526 20170213
% 20170526 20170219
% 20170526 20170303
% 20170526 20170315
% ...
%param_list:run parmslist_generator, and fill the parameters manually
%------------------------------------------------------------------

%SHP selection using SHP_SelPoint_v2, which is same with SHP_SelPoint but
%different O/I format

%you can adjust parameters manually in param_list for SHP selection,
%including:
%
% shp_calwin:	 	sliding window size [pixels]
% shp_method:       method for homogeneous pixel selection
% shp_alpha:        significance level of hypothesis test
%
%for example: I adjust shp_calwin from default [15 15] to [11 21], due to
%azimuth resolution is five times larger than range resolution in S1 data.
%the following code can be also used automatically:
fid = fopen('param_list','r') ;
X = fread(fid) ;
fclose(fid);
X = char(X.') ;
Y = strrep(X, '15 15', '11 21') ; %[az rg];
fid = fopen('param_list','w') ;
fwrite(fid,Y);
fclose (fid);

%then run:
SHP_SelPoint_v2

%once shp.ft is generated, one can run PhaseLinking.m to optimize the
%time-series phase, pl_blksize can be manully adjusted in param_list
%file:

% pl_blksize:       block size to process big image series
%the bklsize depends on your memory size. 
%pl_blksize = 3 means divide the data volume (including shp.ft) into 3 blocks uniformly. 
%pl_blksize = 1 implies all data volume will be inputed into memory completely.  
PhaseLinking(datalistname);
%plot results
paramslistname = 'param_list';
params = readDSIpro(paramslistname);
n_image=1;
diff = freadbkj([diffpath,filesep,num2str(datalist(n_image,1)),'_',num2str(datalist(n_image,2)),'.diff'],params.nlines,'cpxfloat32','b');
pl   = freadbkj([params.pl_path,filesep,num2str(datalist(n_image,1)),'_',num2str(datalist(n_image,2)),'.pl'],params.nlines,'cpxfloat32','b');
figure;subplot(1,2,1);imagesc(angle(diff));
t=title ('Original phase');
set(t,'fontweight','bold');
subplot(1,2,2);imagesc(angle(pl));
t=title ('Optimal phase');
set(t,'fontweight','bold');
colormap jet
%refine spatial network using shortest path algorithm
%four parameters can be manully adjusted in param_list
% pl_bro:           threshold of number of brother pixels [pixels]
% pl_pcoh:          threshold of posteriori coherence
% apsp_pixsize: 	resampling interval for spatial network construction
% apsp_dist:		maximum distance between two spatial pixels [m]
DsC2Spat(datalistname);

tip;




