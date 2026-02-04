clc;clear;close all;
%此脚本用于估计相干性、相位和去斑，支持多视（N-looks）和单视（single-look）数据

%在 MATLAB 2016a, 2018b, 2020b, 2021a (Linux 和 Windows 操作系统) 上测试通过

%----------------------------多视（N-LOOKS）-----------------------------

% 5视时间序列图像
curpath   = pwd; %获取当前路径
mlipath   = [curpath,filesep,'ML',filesep,'MLI']; % 多视强度数据路径
diffpath  = [curpath,filesep,'ML',filesep,'DIFF']; % 多视干涉数据路径（去除相位梯度，如topo和disp）

%读取参数文件
tag_files = dir([mlipath,filesep, '*.par']); % 查找多视 MLI 数据路径下的参数文件
fname     = [mlipath,filesep,tag_files(1).name]; % 获取参数文件名
nlines    = readparam('azimuth_lines',fname); %从参数文件中读取行数

%数据输入
mlistack =ImgRead(mlipath,'mli',nlines,'float32'); %读取多视 MLI 数据
intfstack =ImgRead(diffpath,'diff',nlines,'cpxfloat32'); %读取多视 DIFF 数据

%方位向 1 * 距离向 5 视图像的窗口大小
CalWin = [11 11]; %[方位向 距离向]，5视下的11*11，相当于单视下的 [11 55] 窗口大小，因为预处理时做了距离向 5 视平均
Alpha=0.05; %显著性水平
EstAgr='BWS'; %多视数据集的 SHP 选择仅支持 BWS，单视数据集支持 FaSHPS 和 BWS 两种方法
[SHP_ml]=SHP_SelPoint(mlistack.datastack,CalWin,Alpha,EstAgr);

%Boxcar 相干性估计 
[BoxCoh_ml,BoxPh_ml] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[3 3]);

%自适应相干性估计 
Biascorr='y';%进行偏差校正
[UbCoh_ml,AdpPh_ml] = AdpCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,SHP_ml,Biascorr);

%计算数据大小并输出
fprintf('data size in workspace under n-looks:\n')
tmp = whos('intfstack');
fprintf('%d\t%s %8.2f M\n',1,'intfstack',tmp.bytes/1024/1024);
tmp = whos('mlistack');
fprintf('%d\t%s %8.2f M\n',2,'mlistack ',tmp.bytes/1024/1024);
tmp = whos('SHP_ml');
fprintf('%d\t%s %8.2f M\n',3,'SHP_ml   ',tmp.bytes/1024/1024);
tmp = whos('BoxCoh_ml');
fprintf('%d\t%s %8.2f M\n',4,'BoxCoh_ml',tmp.bytes/1024/1024);
tmp = whos('BoxPh_ml');
fprintf('%d\t%s %8.2f M\n',5,'BoxPh_ml ',tmp.bytes/1024/1024);
tmp = whos('UbCoh_ml');
fprintf('%d\t%s %8.2f M\n',6,'UbCoh_ml ',tmp.bytes/1024/1024);
tmp = whos('AdpPh_ml');
fprintf('%d\t%s %8.2f M\n',7,'AdpPh_ml ',tmp.bytes/1024/1024);

%----------------------------单视（SINGLE-LOOK）-----------------------------
mlipath  = [curpath,filesep,'SL',filesep,'MLI']; 
diffpath = [curpath,filesep,'SL',filesep,'DIFF'];

tag_files = dir([mlipath,filesep, '*.par']);
fname     = [mlipath,filesep,tag_files(1).name];
nlines    = readparam('azimuth_lines',fname); %从参数文件中读取行数

mlistack =ImgRead(mlipath,'mli',nlines,'float32');
intfstack =ImgRead(diffpath,'diff',nlines,'cpxfloat32');

CalWin = [11 21]; %[方位向 距离向]
Alpha=0.05;
EstAgr='FaSHPS';
[SHP_sl]=SHP_SelPoint(mlistack.datastack,CalWin,Alpha,EstAgr);

[BoxCoh_sl,BoxPh_sl] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[7 7]);

%自适应相干性估计 
Biascorr='y';%进行偏差校正
[UbCoh_sl,AdpPh_sl] = AdpCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,SHP_sl,Biascorr);

%计算数据大小并输出
fprintf('data size in workspace under single-look:\n')
tmp = whos('intfstack');
fprintf('%d\t%s %8.2f M\n',1,'intfstack',tmp.bytes/1024/1024);
tmp = whos('mlistack');
fprintf('%d\t%s %8.2f M\n',2,'mlistack ',tmp.bytes/1024/1024);
tmp = whos('SHP_sl');
fprintf('%d\t%s %8.2f M\n',3,'SHP_sl   ',tmp.bytes/1024/1024);
tmp = whos('BoxCoh_sl');
fprintf('%d\t%s %8.2f M\n',4,'BoxCoh_sl',tmp.bytes/1024/1024);
tmp = whos('BoxPh_sl');
fprintf('%d\t%s %8.2f M\n',5,'BoxPh_sl ',tmp.bytes/1024/1024);
tmp = whos('UbCoh_sl');
fprintf('%d\t%s %8.2f M\n',6,'UbCoh_sl ',tmp.bytes/1024/1024);
tmp = whos('AdpPh_sl');
fprintf('%d\t%s %8.2f M\n',7,'AdpPh_sl ',tmp.bytes/1024/1024);

%去斑
BOXMli = DeSpeckling(mlistack.datastack(:,:,1:5));
ADPMli = DeSpeckling(mlistack.datastack(:,:,1:5),SHP_sl);

%----------------------------绘图-----------------------------
% A. 时间序列的平均相干性
% 1. 5视 
figure(1)
clf;
subplot(2,2,1);imagesc(mean(BoxCoh_ml,3),[0 1]);
t=title ('(a)Boxcar Coherence (5-looks)');
set(t,'fontweight','bold');
%可以使用 SHP 数量移除类似 PS 的点
% thre = SHP_ml.BroNum <=4;
% meancoh_ml = mean(UbCoh_ml,3);
% meancoh_ml(thre) = nan;
subplot(2,2,2);imagesc(mean(UbCoh_ml,3),[0 1]);
t=title ('(b)Adp. Coherence (5-looks)');
set(t,'fontweight','bold');
colormap jet

% 2. 单视 
subplot(2,2,3);imagesc(mean(BoxCoh_sl,3),[0 1]);
t=title ('(c)Boxcar Coherence (1-look)');
set(t,'fontweight','bold');
%可以使用 SHP 数量移除类似 PS 的点
% thre = SHP_sl.BroNum <=20;
% meancoh_sl = mean(UbCoh_sl,3);
% meancoh_sl(thre) = nan;
subplot(2,2,4);imagesc(mean(UbCoh_sl,3),[0 1]);
t=title ('(d)Adp. Coherence (1-look)');
set(t,'fontweight','bold');
colormap jet


%B. 去斑（以单视图像为例）
figure(2)
clf;
npage =1;sc=0.9;ex=.8;
subplot(3,1,1);mlishow(BOXMli(:,:,npage),sc,ex);axis image
t=title ('(a)Boxcar Despeckling');
set(t,'fontweight','bold');
subplot(3,1,2);mlishow(ADPMli(:,:,npage),sc,ex);axis image
t=title ('(b)Adp. Despeckling');
set(t,'fontweight','bold');
subplot(3,1,3);mlishow(mlistack.datastack(:,:,npage),sc,ex);axis image
t=title ('(c)Full Speckle');
set(t,'fontweight','bold');

%C. 相位 
figure(3)
clf;
npage =1;
% 1. 5视 
subplot(1,2,1);imagesc(BoxPh_ml(:,:,npage));axis image
t=title ('(a)Boxcar Phase (5-looks)');
set(t,'fontweight','bold');

subplot(1,2,2);imagesc(AdpPh_ml(:,:,npage));axis image
t=title ('(b)Adp. Phase (5-looks)');
set(t,'fontweight','bold');
colormap jet

% 2. 单视 
figure(4)
clf;
subplot(3,1,1);imagesc(BoxPh_sl(:,:,npage));axis image
t=title ('(a)Boxcar Phase (1-look)');
set(t,'fontweight','bold');

subplot(3,1,2);imagesc(AdpPh_sl(:,:,npage));axis image
t=title ('(b)Adp. Phase (1-look)');
set(t,'fontweight','bold');

subplot(3,1,3);imagesc(angle(intfstack.datastack(:,:,npage)));axis image
t=title ('(c)Original Phase (1-look)');
set(t,'fontweight','bold');
colormap jet

%D. SHP（以 5 视图像为例） 
figure(5);
clf;
imagesc(SHP_ml.BroNum);
colormap('jet');
axis image off;
colorbar
t=title ('Homogeneous Pixel Number');
set(t,'fontweight','bold');

%给定参考像素的 SHP 足迹
mlipath   = [curpath,filesep,'ML',filesep,'MLI'];
basemapname = 'mli_ave.bmp';
SHP_FootPrint([mlipath,filesep,basemapname],SHP_ml,[643,250]);

%输出到特定目录，
imgpath = [pwd,filesep,'COH'];
ImgWrite(UbCoh_sl,intfstack.filename,imgpath,'cc','float32','b');
imgpath = [pwd,filesep,'DIFFSM'];
ImgWrite(exp(1j*AdpPh_sl),intfstack.filename,imgpath,'diff.sm','cpxfloat32','b');
tip;
