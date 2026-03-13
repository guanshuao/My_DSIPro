clc;close all;clear;

slcstack = ImgRead('/sar/guanshuao/Beijing_Sentinel/SL/SLC','slc',6400,'cpxfloat32', 'b');
powerstack = slcstack; 
powerstack.datastack = abs(slcstack.datastack).^2;
clear slcstack;


CalWin = [21 21];
Alpha=0.05;
EstAgr='BWS'; 
NumWorkers=48; 

[SHP_BWS_21]=SHP_SelPoint(powerstack.datastack, CalWin, Alpha, EstAgr, NumWorkers, false);
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_21.mat', 'SHP_BWS_21', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_21.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_21.png'); 
delete(gcp('nocreate'));

SHP_BWS_19 = SHP_Resize(SHP_BWS_21, [19 19]);
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_19.mat', 'SHP_BWS_19', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_19.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_19.png');
clear SHP_BWS_21;

SHP_BWS_17 = SHP_Resize(SHP_BWS_19, [17 17]);
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_17.mat', 'SHP_BWS_17', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_17.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_17.png');
clear SHP_BWS_19;

SHP_BWS_15 = SHP_Resize(SHP_BWS_17, [15 15]);
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_15.mat', 'SHP_BWS_15', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_15.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_15.png');
clear SHP_BWS_17;

SHP_BWS_13 = SHP_Resize(SHP_BWS_15, [13 13]);
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_13.mat', 'SHP_BWS_13', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_13.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_13.png');
clear SHP_BWS_15;

SHP_BWS_11 = SHP_Resize(SHP_BWS_13, [11 11]);
save('/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_11.mat', 'SHP_BWS_11', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_11.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/not_connection/SHP_BWS_11.png');
clear SHP_BWS_13 SHP_BWS_11;


[SHP_BWS_21]=SHP_SelPoint(powerstack.datastack, CalWin, Alpha, EstAgr, NumWorkers, true);
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_21.mat', 'SHP_BWS_21', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_21.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_21.png'); 
delete(gcp('nocreate'));

SHP_BWS_19 = SHP_Resize(SHP_BWS_21, [19 19]);
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_19.mat', 'SHP_BWS_19', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_19.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_19.png');
clear SHP_BWS_21;

SHP_BWS_17 = SHP_Resize(SHP_BWS_19, [17 17]);
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_17.mat', 'SHP_BWS_17', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_17.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_17.png');
clear SHP_BWS_19;

SHP_BWS_15 = SHP_Resize(SHP_BWS_17, [15 15]);
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_15.mat', 'SHP_BWS_15', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_15.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_15.png');
clear SHP_BWS_17;

SHP_BWS_13 = SHP_Resize(SHP_BWS_15, [13 13]);
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_13.mat', 'SHP_BWS_13', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_13.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_13.png');
clear SHP_BWS_15;

SHP_BWS_11 = SHP_Resize(SHP_BWS_13, [11 11]);
save('/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_11.mat', 'SHP_BWS_11', '-v7.3');
figure('Visible','off'); imagesc(SHP_BWS_11.BroNum); colorbar; axis image; saveas(gcf,'/sar/guanshuao/Beijing_Sentinel/SL/connection/SHP_BWS_11.png');
clear SHP_BWS_13 SHP_BWS_11;