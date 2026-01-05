clc;clear;close all;
%This script is used to estimate coherence, phase and de-speckling, support
%for both N-looks and single-look data

%tested on MATLAB 2016a, 2018b, 2020b, 2021a (Linux and Windows OS)

%----------------------------N-LOOKS-----------------------------

%5-looks time-series
curpath   = pwd;
mlipath   = [curpath,filesep,'ML',filesep,'MLI'];
diffpath  = [curpath,filesep,'ML',filesep,'DIFF'];

%read parameter file
tag_files = dir([mlipath,filesep, '*.par']);
fname     = [mlipath,filesep,tag_files(1).name]; 
nlines    = readparam('azimuth_lines',fname); 

%mli data input
mlistack =ImgRead(mlipath,'mli',nlines,'float32');
intfstack =ImgRead(diffpath,'diff',nlines,'cpxfloat32');

%window size for az 1*rg 5 looks images
CalWin = [11 11]; %This is equivalent to [11 55] window size under single look
Alpha=0.05;
EstAgr='BWS'; %Only BWS is supported for SHP selection for n-looks dataset
[SHP_ml]=SHP_SelPoint(mlistack.datastack,CalWin,Alpha,EstAgr);

%Boxcar coherence estimation 
[BoxCoh_ml,BoxPh_ml] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[3 3]);

%Adaptive Coherence estimation 
Biascorr='y';%make bias correction
[UbCoh_ml,AdpPh_ml] = AdpCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,SHP_ml,Biascorr);

%count data size
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

%----------------------------SINGLE-LOOK-----------------------------
%single look time-series
mlipath  = [curpath,filesep,'SL',filesep,'MLI']; 
diffpath = [curpath,filesep,'SL',filesep,'DIFF'];

%read parameter file
tag_files = dir([mlipath,filesep, '*.par']);
fname     = [mlipath,filesep,tag_files(1).name];
nlines    = readparam('azimuth_lines',fname); %read nline from parameter file

%mli data input
mlistack =ImgRead(mlipath,'mli',nlines,'float32');
intfstack =ImgRead(diffpath,'diff',nlines,'cpxfloat32');

%window size for az1*rg5 looks images
CalWin = [11 21]; %[az rg]
Alpha=0.05;
EstAgr='FaSHPS';
[SHP_sl]=SHP_SelPoint(mlistack.datastack,CalWin,Alpha,EstAgr);

%Boxcar coherence estimation 
[BoxCoh_sl,BoxPh_sl] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[7 7]);

%Adaptive Coherence estimation 
Biascorr='y';%make bias correction
[UbCoh_sl,AdpPh_sl] = AdpCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,SHP_sl,Biascorr);

%count data size
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

%Despeckling
BOXMli = DeSpeckling(mlistack.datastack(:,:,1:5));
ADPMli = DeSpeckling(mlistack.datastack(:,:,1:5),SHP_sl);

%----------------------------PLOT-----------------------------
% A. MEAN COHERENCE OF TIME-SERIES
% 1. 5-looks 
figure(1)
clf;
subplot(2,2,1);imagesc(mean(BoxCoh_ml,3),[0 1]);
t=title ('(a)Boxcar Coherence (5-looks)');
set(t,'fontweight','bold');
%May remove ps-like points using SHPs number
% thre = SHP_ml.BroNum <=4;
% meancoh_ml = mean(UbCoh_ml,3);
% meancoh_ml(thre) = nan;
subplot(2,2,2);imagesc(mean(UbCoh_ml,3),[0 1]);
t=title ('(b)Adp. Coherence (5-looks)');
set(t,'fontweight','bold');
colormap jet

% 2. single look 
subplot(2,2,3);imagesc(mean(BoxCoh_sl,3),[0 1]);
t=title ('(c)Boxcar Coherence (1-look)');
set(t,'fontweight','bold');
%May remove ps-like points using SHPs number
% thre = SHP_sl.BroNum <=20;
% meancoh_sl = mean(UbCoh_sl,3);
% meancoh_sl(thre) = nan;
subplot(2,2,4);imagesc(mean(UbCoh_sl,3),[0 1]);
t=title ('(d)Adp. Coherence (1-look)');
set(t,'fontweight','bold');
colormap jet


%B. DESPECKLING (using single look image as an example)
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

%C. PHASE 
figure(3)
clf;
npage =1;
% 1. 5-looks 
subplot(1,2,1);imagesc(BoxPh_ml(:,:,npage));axis image
t=title ('(a)Boxcar Phase (5-looks)');
set(t,'fontweight','bold');

subplot(1,2,2);imagesc(AdpPh_ml(:,:,npage));axis image
t=title ('(b)Adp. Phase (5-looks)');
set(t,'fontweight','bold');
colormap jet

% 2. single look 
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

%D. SHP (using 5-looks image as an example) 
figure(5);
clf;
imagesc(SHP_ml.BroNum);
colormap('jet');
axis image off;
colorbar
t=title ('Homogeneous Pixel Number');
set(t,'fontweight','bold');

%SHP footprint for a given reference pixel
mlipath   = [curpath,filesep,'ML',filesep,'MLI'];
basemapname = 'mli_ave.bmp';
SHP_FootPrint([mlipath,filesep,basemapname],SHP_ml,[643,250]);

%output to a specific directory,
imgpath = [pwd,filesep,'COH'];
ImgWrite(UbCoh_sl,intfstack.filename,imgpath,'cc','float32','b');
imgpath = [pwd,filesep,'DIFFSM'];
ImgWrite(exp(1j*AdpPh_sl),intfstack.filename,imgpath,'diff.sm','cpxfloat32','b');
tip;




