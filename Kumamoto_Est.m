clc;close all;

% 以下为开始时的代码，速度较慢
% datesDir = '/titan/guanshuao/Kumamoto/process/dates_resampled';
% pairsDir = '/titan/guanshuao/Kumamoto/process/pairs';
% subset = [201, 25200, 14001, 23000];


% [slcstack, ~] = readISCE2SLCStack(datesDir, 'Subset', subset);
% mlistack = slcstack; mlistack.datastack = abs(slcstack.datastack).^2;
% clear slcstack;

% [intfstack, ~] = readISCE2intStack(pairsDir, 'Subset', subset);
% ImgWrite(intfstack.datastack,intfstack.filename,'/sar/guanshuao/Kumamoto/DIFF','int','cpxfloat32', 'b');
% intfstack.datastack = intfstack.datastack ./ abs(intfstack.datastack);


if exist('powerstack', 'var')
    display('powerstack exists, skip reading SLC stack');
else
    slcstack = ImgRead('/sar/guanshuao/Kumamoto/SLC','slc',25000,'cpxfloat32','b');
    powerstack = slcstack; powerstack.datastack = abs(slcstack.datastack).^2;
    clear slcstack;
    display('SLC stack reading completed');
end


if exist('intfstack', 'var')
    display('intfstack exists, skip reading interferogram stack');
else
    intfstack = ImgRead('/sar/guanshuao/Kumamoto/DIFF','int',25000,'cpxfloat32','b');
    intfstack.datastack = intfstack.datastack ./ abs(intfstack.datastack);
    display('Interferogram stack reading completed');
end



delete(gcp('nocreate'));
NumWorkers=25;
Biascorr='n';%make bias correction





load('/sar/guanshuao/Kumamoto/Data/connection/SHP_BWS_23.mat');
[AdpCoh_BWS_23] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_23, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_23); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_23.png');
save('/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_23.mat', 'AdpCoh_BWS_23', '-v7.3');
clear AdpCoh_BWS_23;clear SHP_BWS_23;

load('/sar/guanshuao/Kumamoto/Data/connection/SHP_BWS_25.mat');
[AdpCoh_BWS_25] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_25, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_25); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_25.png');
save('/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_25.mat', 'AdpCoh_BWS_25', '-v7.3');
clear AdpCoh_BWS_25;clear SHP_BWS_25;

load('/sar/guanshuao/Kumamoto/Data/connection/SHP_BWS_27.mat');
[AdpCoh_BWS_27] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_27, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_27); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_27.png');
save('/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_27.mat', 'AdpCoh_BWS_27', '-v7.3');
clear AdpCoh_BWS_27;clear SHP_BWS_27;

load('/sar/guanshuao/Kumamoto/Data/connection/SHP_BWS_29.mat');
[AdpCoh_BWS_29] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_29, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_29); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_29.png');
save('/sar/guanshuao/Kumamoto/Data/connection/AdpCoh_BWS_29.mat', 'AdpCoh_BWS_29', '-v7.3');
clear AdpCoh_BWS_29;clear SHP_BWS_29;



load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_13.mat');
[AdpCoh_BWS_13] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_13, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_13); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_13.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_13.mat', 'AdpCoh_BWS_13', '-v7.3');
clear AdpCoh_BWS_13;clear SHP_BWS_13;



load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_17.mat');
[AdpCoh_BWS_17] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_17, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_17); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_17.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_17.mat', 'AdpCoh_BWS_17', '-v7.3');
clear AdpCoh_BWS_17;clear SHP_BWS_17;

load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_19.mat');
[AdpCoh_BWS_19] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_19, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_19); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_19.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_19.mat', 'AdpCoh_BWS_19', '-v7.3');
clear AdpCoh_BWS_19;clear SHP_BWS_19;

load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_21.mat');
[AdpCoh_BWS_21] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_21, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_21); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_21.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_21.mat', 'AdpCoh_BWS_21', '-v7.3');
clear AdpCoh_BWS_21;clear SHP_BWS_21;

load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_23.mat');
[AdpCoh_BWS_23] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_23, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_23); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_23.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_23.mat', 'AdpCoh_BWS_23', '-v7.3');
clear AdpCoh_BWS_23;clear SHP_BWS_23;

load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_25.mat');
[AdpCoh_BWS_25] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_25, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_25); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_25.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_25.mat', 'AdpCoh_BWS_25', '-v7.3');
clear AdpCoh_BWS_25;clear SHP_BWS_25;

load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_27.mat');
[AdpCoh_BWS_27] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_27, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_27); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_27.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_27.mat', 'AdpCoh_BWS_27', '-v7.3');
clear AdpCoh_BWS_27;clear SHP_BWS_27;

load('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_29.mat');
[AdpCoh_BWS_29] = AdpCohEst(powerstack.datastack, powerstack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_29, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_29); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_29.png');
save('/sar/guanshuao/Kumamoto/Data/not_connection/AdpCoh_BWS_29.mat', 'AdpCoh_BWS_29', '-v7.3');
clear AdpCoh_BWS_29;clear SHP_BWS_29;

% Boxcar coherence estimation 
[BoxCoh_3] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[3 3], 'average');
figure('Visible','off'); imagesc(BoxCoh_3); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_3.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_3.mat', 'BoxCoh_3', '-v7.3');
clear BoxCoh_3;

[BoxCoh_5] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[5 5], 'average');
figure('Visible','off'); imagesc(BoxCoh_5); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_5.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_5.mat', 'BoxCoh_5', '-v7.3');
clear BoxCoh_5;

[BoxCoh_7] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[7 7], 'average');
figure('Visible','off'); imagesc(BoxCoh_7); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_7.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_7.mat', 'BoxCoh_7', '-v7.3');
clear BoxCoh_7;

[BoxCoh_9] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[9 9], 'average');
figure('Visible','off'); imagesc(BoxCoh_9); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_9.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_9.mat', 'BoxCoh_9', '-v7.3');
clear BoxCoh_9;

[BoxCoh_11] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[11 11], 'average');
figure('Visible','off'); imagesc(BoxCoh_11); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_11.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_11.mat', 'BoxCoh_11', '-v7.3');
clear BoxCoh_11;

[BoxCoh_13] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[13 13], 'average');
figure('Visible','off'); imagesc(BoxCoh_13); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_13.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_13.mat', 'BoxCoh_13', '-v7.3');
clear BoxCoh_13;

[BoxCoh_15] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[15 15], 'average');
figure('Visible','off'); imagesc(BoxCoh_15); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_15.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_15.mat', 'BoxCoh_15', '-v7.3');
clear BoxCoh_15;

[BoxCoh_17] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[17 17], 'average');
figure('Visible','off'); imagesc(BoxCoh_17); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_17.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_17.mat', 'BoxCoh_17', '-v7.3');
clear BoxCoh_17;

[BoxCoh_19] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[19 19], 'average');
figure('Visible','off'); imagesc(BoxCoh_19); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_19.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_19.mat', 'BoxCoh_19', '-v7.3');
clear BoxCoh_19;

[BoxCoh_21] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[21 21], 'average');
figure('Visible','off'); imagesc(BoxCoh_21); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_21.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_21.mat', 'BoxCoh_21', '-v7.3');
clear BoxCoh_21;

[BoxCoh_23] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[23 23], 'average');
figure('Visible','off'); imagesc(BoxCoh_23); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_23.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_23.mat', 'BoxCoh_23', '-v7.3');
clear BoxCoh_23;

[BoxCoh_25] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[25 25], 'average');
figure('Visible','off'); imagesc(BoxCoh_25); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_25.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_25.mat', 'BoxCoh_25', '-v7.3');
clear BoxCoh_25;

[BoxCoh_27] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[27 27], 'average');
figure('Visible','off'); imagesc(BoxCoh_27); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_27.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_27.mat', 'BoxCoh_27', '-v7.3');
clear BoxCoh_27;

[BoxCoh_29] = BoxCohEst(powerstack.datastack,powerstack.filename,intfstack.datastack,intfstack.filename,[29 29], 'average');
figure('Visible','off'); imagesc(BoxCoh_29); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Kumamoto/Data/BoxCoh_29.png');
save('/sar/guanshuao/Kumamoto/Data/BoxCoh_29.mat', 'BoxCoh_29', '-v7.3');