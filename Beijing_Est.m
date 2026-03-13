clc;close all;clear;

NumWorkers=48;
Biascorr='n';%make bias correction



load('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_11.mat');
[AdpCoh_BWS_11] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_11, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_11); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_11.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_11.mat', 'AdpCoh_BWS_11', '-v7.3');
clear AdpCoh_BWS_11;clear SHP_BWS_11;

load('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_13.mat');
[AdpCoh_BWS_13] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_13, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_13); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_13.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_13.mat', 'AdpCoh_BWS_13', '-v7.3');
clear AdpCoh_BWS_13;clear SHP_BWS_13;

load('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_15.mat');
[AdpCoh_BWS_15] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_15, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_15); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_15.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_15.mat', 'AdpCoh_BWS_15', '-v7.3');
clear AdpCoh_BWS_15;clear SHP_BWS_15;

load('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_17.mat');
[AdpCoh_BWS_17] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_17, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_17); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_17.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_17.mat', 'AdpCoh_BWS_17', '-v7.3');
clear AdpCoh_BWS_17;clear SHP_BWS_17;

load('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_19.mat');
[AdpCoh_BWS_19] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_19, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_19); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_19.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_19.mat', 'AdpCoh_BWS_19', '-v7.3');
clear AdpCoh_BWS_19;clear SHP_BWS_19;

load('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_21.mat');
[AdpCoh_BWS_21] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_21, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_21); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_21.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/AdpCoh_BWS_21.mat', 'AdpCoh_BWS_21', '-v7.3');
clear AdpCoh_BWS_21;clear SHP_BWS_21;

load('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_11.mat');
[AdpCoh_BWS_11] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_11, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_11); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_11.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_11.mat', 'AdpCoh_BWS_11', '-v7.3');
clear AdpCoh_BWS_11;clear SHP_BWS_11;

load('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_13.mat');
[AdpCoh_BWS_13] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_13, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_13); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_13.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_13.mat', 'AdpCoh_BWS_13', '-v7.3');
clear AdpCoh_BWS_13;clear SHP_BWS_13;

load('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_15.mat');
[AdpCoh_BWS_15] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_15, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_15); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_15.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_15.mat', 'AdpCoh_BWS_15', '-v7.3');
clear AdpCoh_BWS_15;clear SHP_BWS_15;

load('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_17.mat');
[AdpCoh_BWS_17] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_17, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_17); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_17.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_17.mat', 'AdpCoh_BWS_17', '-v7.3');
clear AdpCoh_BWS_17;clear SHP_BWS_17;

load('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_19.mat');
[AdpCoh_BWS_19] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_19, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_19); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_19.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_19.mat', 'AdpCoh_BWS_19', '-v7.3');
clear AdpCoh_BWS_19;clear SHP_BWS_19;

load('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_21.mat');
[AdpCoh_BWS_21] = AdpCohEst_New('/sar/guanshuao/Beijing_Sentinel/SL/SLC', 'slc', 'cpxfloat32', 'b', '/sar/guanshuao/Beijing_Sentinel/SL/DIFF', 'int', 'cpxfloat32','b', 6400, SHP_BWS_21, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_21); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_21.png');
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/AdpCoh_BWS_21.mat', 'AdpCoh_BWS_21', '-v7.3');
clear AdpCoh_BWS_21;clear SHP_BWS_21;