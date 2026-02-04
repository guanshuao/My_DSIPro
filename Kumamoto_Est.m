clc;clear;close all;


datesDir = '/titan/guanshuao/Kumamoto/process/dates_resampled';
pairsDir = '/titan/guanshuao/Kumamoto/process/pairs';
subset = [201, 25200, 14001, 23000];



% 多窗口一次性估计，避免重复读盘
calWinList = [29 29; 27 27; 25 25; 23 23; 21 21; 19 19; 17 17; 15 15; 13 13; 11 11; 9 9; 7 7; 5 5; 3 3];
boxOutDirs = arrayfun(@(i) fullfile(pwd, sprintf('BoxCoh_%d', calWinList(i,1))), 1:size(calWinList,1), 'UniformOutput', false);
for i = 1:numel(boxOutDirs)
	if ~exist(boxOutDirs{i},'dir'), mkdir(boxOutDirs{i}); end
end

% 逐对读取并计算 Boxcar 相干性（每对数据只读一次）
BoxCohEst_StreamByPair(datesDir, pairsDir, 'Subset', subset, 'CalWin', calWinList, 'CohEstAgr', 'fast', 'OutDir', boxOutDirs);


%% 自适应相干性估计
load('/titan/guanshuao/Kumamoto/Data/SHP_BWS_21.mat');
whos;
Biascorr = 'n';
adpOutDir = '/titan/guanshuao/Kumamoto/Data/AdpCoh';
if ~exist(adpOutDir,'dir'), mkdir(adpOutDir); end
AdpCohEst_StreamByPair(datesDir, pairsDir, SHP_BWS_21, Biascorr, adpOutDir, 'Subset', subset);
