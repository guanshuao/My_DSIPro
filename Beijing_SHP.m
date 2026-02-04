clc;close all;

% curpath   = pwd;
mlipath   = ['/titan/guanshuao/Beijing_Sentinel',filesep,'ML',filesep,'MLI'];
diffpath  = ['/titan/guanshuao/Beijing_Sentinel',filesep,'ML',filesep,'DIFF'];

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

Alpha=0.05;
EstAgr='BWS'; 
NumWorkers=40;
c = parcluster('local');
c.NumWorkers = NumWorkers;        
parpool(c, NumWorkers);      

CalWin = [21 21];
[SHP_BWS_21]=SHP_SelPoint_Parallel(mlistack.datastack, CalWin, Alpha, EstAgr, NumWorkers);
save('/titan/guanshuao/Beijing_Sentinel/ML/SHP_BWS_21.mat', 'SHP_BWS_21', '-v7.3');

SHP_BWS_19 = SHP_Resize(SHP_BWS_21, [19 19]);
save('/titan/guanshuao/Beijing_Sentinel/ML/SHP_BWS_19.mat', 'SHP_BWS_19', '-v7.3');

SHP_BWS_17 = SHP_Resize(SHP_BWS_21, [17 17]);
save('/titan/guanshuao/Beijing_Sentinel/ML/SHP_BWS_17.mat', 'SHP_BWS_17', '-v7.3');

SHP_BWS_15 = SHP_Resize(SHP_BWS_21, [15 15]);
save('/titan/guanshuao/Beijing_Sentinel/ML/SHP_BWS_15.mat', 'SHP_BWS_15', '-v7.3');

SHP_BWS_13 = SHP_Resize(SHP_BWS_21, [13 13]);
save('/titan/guanshuao/Beijing_Sentinel/ML/SHP_BWS_13.mat', 'SHP_BWS_13', '-v7.3');

SHP_BWS_11 = SHP_Resize(SHP_BWS_21, [11 11]);
save('/titan/guanshuao/Beijing_Sentinel/ML/SHP_BWS_11.mat', 'SHP_BWS_11', '-v7.3');
