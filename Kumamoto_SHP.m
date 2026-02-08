clear;close all;clc;

[slcstack, ~] = readISCE2SLCStack('/titan/guanshuao/Kumamoto/process/dates_resampled', 'Subset', [201, 25200, 14001, 23000]);
mlistack = slcstack; mlistack.datastack = abs(slcstack.datastack).^2;
clear slcstack;


Alpha = 0.05;
EstAgr ='BWS'; 
NumBlocks = 60;

c = parcluster('local');
c.NumWorkers = NumBlocks;        % 允许的最大 worker 数
parpool(c, NumBlocks);           % 或 parpool('local', NumBlocks)

CalWin = [21 21]; 
[SHP_BWS_21]=SHP_SelPoint_Parallel(mlistack.datastack, CalWin, Alpha, EstAgr, NumBlocks);
save('/titan/guanshuao/Kumamoto/Data/SHP_BWS_21.mat', 'SHP_BWS_21', '-v7.3');
clear SHP_BWS_21;