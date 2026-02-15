close all;clc;

if exist('mlistack', 'var')
    display('mlistack exists, skip reading SLC stack');
else
    [slcstack, ~] = readISCE2SLCStack('/titan/guanshuao/Kumamoto/process/dates_resampled', 'Subset', [201, 25200, 14001, 23000]);
    mlistack = slcstack; mlistack.datastack = abs(slcstack.datastack).^2;
    clear slcstack;
    display('SLC stack reading completed');
end


Alpha = 0.05;
EstAgr ='BWS'; 
NumBlocks = 40;
Connection = false;
delete(gcp('nocreate'));


CalWin = [29 29]; 
[SHP_BWS_29]=SHP_SelPoint(mlistack.datastack, CalWin, Alpha, EstAgr, NumBlocks, Connection);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_29.mat', 'SHP_BWS_29', '-v7.3');

SHP_BWS_27 = SHP_Resize(SHP_BWS_29, [27, 27]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_27.mat', 'SHP_BWS_27', '-v7.3');
clear SHP_BWS_29;

SHP_BWS_25 = SHP_Resize(SHP_BWS_27, [25, 25]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_25.mat', 'SHP_BWS_25', '-v7.3');
clear SHP_BWS_27;

SHP_BWS_23 = SHP_Resize(SHP_BWS_25, [23, 23]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_23.mat', 'SHP_BWS_23', '-v7.3');
clear SHP_BWS_25;

SHP_BWS_21 = SHP_Resize(SHP_BWS_23, [21, 21]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_21.mat', 'SHP_BWS_21', '-v7.3');
clear SHP_BWS_23;

SHP_BWS_19 = SHP_Resize(SHP_BWS_21, [19, 19]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_19.mat', 'SHP_BWS_19', '-v7.3');
clear SHP_BWS_21;

SHP_BWS_17 = SHP_Resize(SHP_BWS_19, [17, 17]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_17.mat', 'SHP_BWS_17', '-v7.3');
clear SHP_BWS_19;

SHP_BWS_15 = SHP_Resize(SHP_BWS_17, [15, 15]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_15.mat', 'SHP_BWS_15', '-v7.3');
clear SHP_BWS_17;

SHP_BWS_13 = SHP_Resize(SHP_BWS_15, [13, 13]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_13.mat', 'SHP_BWS_13', '-v7.3');
clear SHP_BWS_15;

SHP_BWS_11 = SHP_Resize(SHP_BWS_13, [11, 11]);
save('/sar/guanshuao/Kumamoto/Data/not_connection/SHP_BWS_11.mat', 'SHP_BWS_11', '-v7.3');
clear SHP_BWS_13;