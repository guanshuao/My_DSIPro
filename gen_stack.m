clear;

data_dir = '/sar/guanshuao/Kumamoto/Data/not_connection';

Stack = struct();

Stack.Calwin = [15, 15];

Stack.level_num = 50;

load(fullfile(data_dir, 'AdpCoh_BWS_15.mat'));
Stack.AdpCoh = AdpCoh_BWS_15;
clear AdpCoh_BWS_15;

load(fullfile(data_dir, 'BoxCoh_3.mat'));
Stack.BoxCoh_3 = BoxCoh_3;
clear BoxCoh_3;

load(fullfile(data_dir, 'BoxCoh_5.mat'));
Stack.BoxCoh_5 = BoxCoh_5;
clear BoxCoh_5;

load(fullfile(data_dir, 'BoxCoh_7.mat'));
Stack.BoxCoh_7 = BoxCoh_7;
clear BoxCoh_7;

load(fullfile(data_dir, 'BoxCoh_9.mat'));
Stack.BoxCoh_9 = BoxCoh_9;
clear BoxCoh_9;

load(fullfile(data_dir, 'BoxCoh_11.mat'));
Stack.BoxCoh_11 = BoxCoh_11;
clear BoxCoh_11;

load(fullfile(data_dir, 'BoxCoh_13.mat'));
Stack.BoxCoh_13 = BoxCoh_13;
clear BoxCoh_13;

load(fullfile(data_dir, 'BoxCoh_15.mat'));
Stack.BoxCoh_15 = BoxCoh_15;
clear BoxCoh_15;

load(fullfile(data_dir, 'CV_15.mat'));
Stack.CV = CV_15;
clear CV_15;

load(fullfile(data_dir, 'RCS_15.mat'));
Stack.RCS = RCS_15;
clear RCS_15;

load(fullfile(data_dir, 'SHP_BWS_15.mat'));
Stack.SHPNum = SHP_BWS_15.BroNum;
Stack.Connection = SHP_BWS_15.Connection;
clear SHP_BWS_15;

Stack.TrueCoh = LUT_lookup('inverse', Stack.SHPNum, Stack.AdpCoh);

Stack.ModifiedCoh = LUT_lookup('forward', 225, Stack.TrueCoh);

find_best_win_est;

Stack.Best_Est = Best_Est;
Stack.Best_Win = Best_Win;

Stack.CV_level = gen_level_map(Stack.CV, 'float', Stack.level_num);

Stack.SHPNum_level = gen_level_map(Stack.SHPNum, 'integer', Stack.level_num);



% 保存结果
save(fullfile(data_dir, 'Stack.mat'), 'Stack', '-v7.3');
