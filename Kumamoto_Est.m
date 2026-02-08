clc;close all;


datesDir = '/titan/guanshuao/Kumamoto/process/dates_resampled';
pairsDir = '/titan/guanshuao/Kumamoto/process/pairs';
subset = [201, 25200, 14001, 23000];


if exist('mlistack', 'var')
    display('mlistack exists, skip reading SLC stack');
else
    [slcstack, ~] = readISCE2SLCStack(datesDir, 'Subset', subset);
    mlistack = slcstack; mlistack.datastack = abs(slcstack.datastack).^2;
    clear slcstack;
    display('SLC stack reading completed');
end

if exist('intfstack', 'var')
    display('intfstack exists, skip reading interferogram stack');
else
    [intfstack, ~] = readISCE2intStack(pairsDir, 'Subset', subset);
    intfstack.datastack = intfstack.datastack ./ abs(intfstack.datastack);
    display('Interferogram stack reading completed');
end


delete(gcp('nocreate'));
NumWorkers=14;

% c = parcluster('local');
% c.NumWorkers = NumWorkers;        
% parpool(c, NumWorkers);
Biascorr='n';%make bias correction

if exist('SHP_BWS_11', 'var')
    display('SHP_BWS_11 exists, skip SHP reading');
else
    load('/titan/guanshuao/Kumamoto/Data/SHP_BWS_11.mat');
end



[AdpCoh_BWS_11] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_11, Biascorr, 'average', NumWorkers);
figure('Visible','off'); imagesc(AdpCoh_BWS_11); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/AdpCoh_BWS_11.png');
save('/titan/guanshuao/Kumamoto/Data/AdpCoh_BWS_11.mat', 'AdpCoh_BWS_11', '-v7.3');
clear AdpCoh_BWS_11;clear SHP_BWS_11;

delete(gcp('nocreate'));

% Boxcar coherence estimation 
[BoxCoh_3] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[3 3], 'average');
figure('Visible','off'); imagesc(BoxCoh_3); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_3.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_3.mat', 'BoxCoh_3', '-v7.3');
clear BoxCoh_3;

[BoxCoh_5] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[5 5], 'average');
figure('Visible','off'); imagesc(BoxCoh_5); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_5.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_5.mat', 'BoxCoh_5', '-v7.3');
clear BoxCoh_5;

[BoxCoh_7] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[7 7], 'average');
figure('Visible','off'); imagesc(BoxCoh_7); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_7.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_7.mat', 'BoxCoh_7', '-v7.3');
clear BoxCoh_7;

[BoxCoh_9] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[9 9], 'average');
figure('Visible','off'); imagesc(BoxCoh_9); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_9.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_9.mat', 'BoxCoh_9', '-v7.3');
clear BoxCoh_9;

[BoxCoh_11] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[11 11], 'average');
figure('Visible','off'); imagesc(BoxCoh_11); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_11.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_11.mat', 'BoxCoh_11', '-v7.3');
clear BoxCoh_11;

[BoxCoh_13] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[13 13], 'average');
figure('Visible','off'); imagesc(BoxCoh_13); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_13.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_13.mat', 'BoxCoh_13', '-v7.3');
clear BoxCoh_13;

[BoxCoh_15] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[15 15], 'average');
figure('Visible','off'); imagesc(BoxCoh_15); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_15.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_15.mat', 'BoxCoh_15', '-v7.3');
clear BoxCoh_15;

[BoxCoh_17] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[17 17], 'average');
figure('Visible','off'); imagesc(BoxCoh_17); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_17.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_17.mat', 'BoxCoh_17', '-v7.3');
clear BoxCoh_17;

[BoxCoh_19] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[19 19], 'average');
figure('Visible','off'); imagesc(BoxCoh_19); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_19.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_19.mat', 'BoxCoh_19', '-v7.3');
clear BoxCoh_19;

[BoxCoh_21] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[21 21], 'average');
figure('Visible','off'); imagesc(BoxCoh_21); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_21.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_21.mat', 'BoxCoh_21', '-v7.3');
clear BoxCoh_21;

[BoxCoh_23] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[23 23], 'average');
figure('Visible','off'); imagesc(BoxCoh_23); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_23.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_23.mat', 'BoxCoh_23', '-v7.3');
clear BoxCoh_23;

[BoxCoh_25] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[25 25], 'average');
figure('Visible','off'); imagesc(BoxCoh_25); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_25.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_25.mat', 'BoxCoh_25', '-v7.3');
clear BoxCoh_25;

[BoxCoh_27] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[27 27], 'average');
figure('Visible','off'); imagesc(BoxCoh_27); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_27.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_27.mat', 'BoxCoh_27', '-v7.3');
clear BoxCoh_27;

[BoxCoh_29] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[29 29], 'average');
figure('Visible','off'); imagesc(BoxCoh_29); colorbar; axis image; saveas(gcf,'/titan/guanshuao/Kumamoto/Data/BoxCoh_29.png');
save('/titan/guanshuao/Kumamoto/Data/BoxCoh_29.mat', 'BoxCoh_29', '-v7.3');
clear;