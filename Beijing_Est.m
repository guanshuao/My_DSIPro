clc;clear;close all;

% curpath   = pwd;
mlipath   = ['/titan/guanshuao/Beijing_Sentinel',filesep,'SL',filesep,'MLI'];
diffpath  = ['/titan/guanshuao/Beijing_Sentinel',filesep,'SL',filesep,'DIFF'];

%read parameter file
tag_files = dir([mlipath,filesep, '*.par']);
fname     = [mlipath,filesep,tag_files(1).name]; 
nlines    = readparam('azimuth_lines',fname); 

slcstack =ImgRead(mlipath,'slc',nlines,'cpxfloat32');
[~,~,num_slcstack] = size(slcstack.datastack);
for i = 1:num_slcstack
    slcstack.datastack(:,:,i) = flip(slcstack.datastack(:,:,i));
end
mlistack = slcstack; mlistack.datastack = abs(slcstack.datastack).^2;
clear slcstack;


intfstack = ImgRead(diffpath,'slc',nlines,'cpxfloat32');
intfstack.datastack = intfstack.datastack ./ abs(intfstack.datastack);
[~,~,num_intfstack] = size(intfstack.datastack);
for i = 1:num_intfstack
    intfstack.datastack(:,:,i) = flip(intfstack.datastack(:,:,i));
end

load('/titan/guanshuao/Beijing_Sentinel/SL/SHP_BWS_21.mat');

Biascorr='n';%make bias correction
[AdpCoh_BWS_21,~] = AdpCohEst(mlistack.datastack, mlistack.filename, intfstack.datastack, intfstack.filename, SHP_BWS_21, Biascorr, 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/AdpCoh_BWS_21.mat', 'AdpCoh_BWS_21', '-v7.3');
clear AdpCoh_BWS_21;clear SHP_BWS_21;

% Boxcar coherence estimation 
[BoxCoh_3,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[3 3], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_3.mat', 'BoxCoh_3', '-v7.3');
clear BoxCoh_3;

[BoxCoh_5,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[5 5], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_5.mat', 'BoxCoh_5', '-v7.3');
clear BoxCoh_5;

[BoxCoh_7,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[7 7], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_7.mat', 'BoxCoh_7', '-v7.3');
clear BoxCoh_7;

[BoxCoh_9,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[9 9], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_9.mat', 'BoxCoh_9', '-v7.3');
clear BoxCoh_9;

[BoxCoh_11,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[11 11], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_11.mat', 'BoxCoh_11', '-v7.3');
clear BoxCoh_11;

[BoxCoh_13,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[13 13], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_13.mat', 'BoxCoh_13', '-v7.3');
clear BoxCoh_13;

[BoxCoh_15,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[15 15], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_15.mat', 'BoxCoh_15', '-v7.3');
clear BoxCoh_15;

[BoxCoh_17,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[17 17], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_17.mat', 'BoxCoh_17', '-v7.3');
clear BoxCoh_17;

[BoxCoh_19,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[19 19], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_19.mat', 'BoxCoh_19', '-v7.3');
clear BoxCoh_19;

[BoxCoh_21,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[21 21], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_21.mat', 'BoxCoh_21', '-v7.3');
clear BoxCoh_21;

[BoxCoh_23,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[23 23], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_23.mat', 'BoxCoh_23', '-v7.3');
clear BoxCoh_23;

[BoxCoh_25,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[25 25], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_25.mat', 'BoxCoh_25', '-v7.3');
clear BoxCoh_25;

[BoxCoh_27,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[27 27], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_27.mat', 'BoxCoh_27', '-v7.3');
clear BoxCoh_27;

[BoxCoh_29,~] = BoxCohEst(mlistack.datastack,mlistack.filename,intfstack.datastack,intfstack.filename,[29 29], 'fast', 'average');
save('/titan/guanshuao/Beijing_Sentinel/SL/BoxCoh_29.mat', 'BoxCoh_29', '-v7.3');
clear BoxCoh_29;