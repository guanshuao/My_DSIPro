clc;close all;

% curpath   = pwd;
mlipath   = ['/titan/guanshuao/Beijing_Sentinel',filesep,'SL',filesep,'MLI'];
diffpath  = ['/titan/guanshuao/Beijing_Sentinel',filesep,'SL',filesep,'DIFF'];

%read parameter file
tag_files = dir([mlipath,filesep, '*.par']);
fname     = [mlipath,filesep,tag_files(1).name]; 
nlines    = readparam('azimuth_lines',fname); 


if exist('mlistack', 'var')
    disp('mlistack变量已存在，跳过SLC堆栈读取步骤');
else
    disp('正在读取SLC堆栈...');
    slcstack = ImgRead(mlipath,'slc',nlines,'cpxfloat32');
    [~,~,num_slcstack] = size(slcstack.datastack);
    for i = 1:num_slcstack
        slcstack.datastack(:,:,i) = flip(slcstack.datastack(:,:,i));
    end
    mlistack = slcstack; 
    mlistack.datastack = abs(slcstack.datastack).^2;
    clear slcstack;
    disp('SLC堆栈读取完成');
end


if exist('intfstack', 'var')
    disp('intfstack变量已存在，跳过干涉图堆栈读取步骤');
else
    disp('正在读取干涉图堆栈...');
    intfstack = ImgRead(diffpath,'slc',nlines,'cpxfloat32');
    intfstack.datastack = intfstack.datastack ./ abs(intfstack.datastack);
    [~,~,num_intfstack] = size(intfstack.datastack);
    for i = 1:num_intfstack
        intfstack.datastack(:,:,i) = flip(intfstack.datastack(:,:,i));
    end
    disp('干涉图堆栈读取完成');
end


NumWorkers=20;
delete(gcp('nocreate'));
c = parcluster('local');
c.NumWorkers = NumWorkers;        
parpool(c, NumWorkers);
Biascorr='n';%make bias correction

load('/titan/guanshuao/Beijing_Sentinel/SL/SHP_BWS_19.mat');
[AdpCoh_BWS_19] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_19, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_19); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_19.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_19.mat', 'AdpCoh_BWS_19', '-v7.3');
clear AdpCoh_BWS_19;clear SHP_BWS_19;

load('/titan/guanshuao/Beijing_Sentinel/SL/SHP_BWS_17.mat');
[AdpCoh_BWS_17] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_17, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_17); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_17.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_17.mat', 'AdpCoh_BWS_17', '-v7.3');
clear AdpCoh_BWS_17;clear SHP_BWS_17;

load('/titan/guanshuao/Beijing_Sentinel/SL/SHP_BWS_15.mat');
[AdpCoh_BWS_15] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_15, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_15); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_15.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_15.mat', 'AdpCoh_BWS_15', '-v7.3');
clear AdpCoh_BWS_15;clear SHP_BWS_15;

load('/titan/guanshuao/Beijing_Sentinel/SL/SHP_BWS_13.mat');
[AdpCoh_BWS_13] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_13, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_13); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_13.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_13.mat', 'AdpCoh_BWS_13', '-v7.3');
clear AdpCoh_BWS_13;clear SHP_BWS_13;

load('/titan/guanshuao/Beijing_Sentinel/SL/SHP_BWS_11.mat');
[AdpCoh_BWS_11] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_11, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_11); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_11.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_11.mat', 'AdpCoh_BWS_11', '-v7.3');
clear AdpCoh_BWS_11;clear SHP_BWS_11;





[BoxCoh_9] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[9 9], 'average');
figure('Visible','off'); imagesc(BoxCoh_9); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_9.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_9.mat', 'BoxCoh_9', '-v7.3');
clear BoxCoh_9;

[BoxCoh_11] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[11 11], 'average');
figure('Visible','off'); imagesc(BoxCoh_11); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_11.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_11.mat', 'BoxCoh_11', '-v7.3');
clear BoxCoh_11;

[BoxCoh_13] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[13 13], 'average');
figure('Visible','off'); imagesc(BoxCoh_13); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_13.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_13.mat', 'BoxCoh_13', '-v7.3');
clear BoxCoh_13;

[BoxCoh_15] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[15 15], 'average');
figure('Visible','off'); imagesc(BoxCoh_15); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_15.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_15.mat', 'BoxCoh_15', '-v7.3');
clear BoxCoh_15;

[BoxCoh_17] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[17 17], 'average');
figure('Visible','off'); imagesc(BoxCoh_17); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_17.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_17.mat', 'BoxCoh_17', '-v7.3');
clear BoxCoh_17;

[BoxCoh_19] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[19 19], 'average');
figure('Visible','off'); imagesc(BoxCoh_19); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_19.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_19.mat', 'BoxCoh_19', '-v7.3');
clear BoxCoh_19;

[BoxCoh_21] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[21 21], 'average');
figure('Visible','off'); imagesc(BoxCoh_21); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_21.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_21.mat', 'BoxCoh_21', '-v7.3');
clear BoxCoh_21;

[BoxCoh_23] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[23 23], 'average');
figure('Visible','off'); imagesc(BoxCoh_23); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_23.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_23.mat', 'BoxCoh_23', '-v7.3');
clear BoxCoh_23;

[BoxCoh_25] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[25 25], 'average');
figure('Visible','off'); imagesc(BoxCoh_25); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_25.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_25.mat', 'BoxCoh_25', '-v7.3');
clear BoxCoh_25;

[BoxCoh_27] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[27 27], 'average');
figure('Visible','off'); imagesc(BoxCoh_27); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_27.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_27.mat', 'BoxCoh_27', '-v7.3');
clear BoxCoh_27;

[BoxCoh_29] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[29 29], 'average');
figure('Visible','off'); imagesc(BoxCoh_29); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_29.png');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_29.mat', 'BoxCoh_29', '-v7.3');
clear BoxCoh_29;