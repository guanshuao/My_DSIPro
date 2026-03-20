close all; clear; clc;

data_dir = "Stack.mat";
disp('Loading Stack.mat...');
load(data_dir);
disp('Loaded Stack.mat...');

level_num = Stack.level_num;
levels = (1:level_num)';
Calwin = Stack.Calwin;
max_win = Calwin(1);
box_coh_field = sprintf('BoxCoh_%d', max_win);
if isfield(Stack, box_coh_field)
    BoxCoh_Max = Stack.(box_coh_field);
else
    BoxCoh_Max = NaN(size(Stack.TrueCoh));
end

fields = fieldnames(Stack);
box_coh_fields = fields(startsWith(fields, 'BoxCoh_'));
box_wins = zeros(length(box_coh_fields), 1);
for i = 1:length(box_coh_fields)
    tokens = regexp(box_coh_fields{i}, 'BoxCoh_(\d+)', 'tokens');
    box_wins(i) = str2double(tokens{1}{1});
end
[box_wins, sort_idx] = sort(box_wins);
box_coh_fields = box_coh_fields(sort_idx);

valid_mask = ~isnan(Stack.CV_level) & (Stack.CV_level >= 1 & Stack.CV_level <= level_num);
groups = Stack.CV_level(valid_mask);
groups = round(groups(:));

calc_mean = @(x) double(mean(x, 'omitnan'));

%% Fig 1
if isfield(Stack, 'AdpCph')
    adp_field = 'AdpCph';
elseif isfield(Stack, 'AdpCoh')
    adp_field = 'AdpCoh';
else
    adp_field = '';
end
has_modified = isfield(Stack, 'ModifiedCoh');

mean_TrueCoh = accumarray(groups, double(Stack.TrueCoh(valid_mask)), [level_num, 1], calc_mean, NaN);
mean_BoxCohMax = accumarray(groups, double(BoxCoh_Max(valid_mask)), [level_num, 1], calc_mean, NaN);
mean_BestEst = accumarray(groups, double(Stack.Best_Est(valid_mask)), [level_num, 1], calc_mean, NaN);
if ~isempty(adp_field)
    mean_AdpCoh = accumarray(groups, double(Stack.(adp_field)(valid_mask)), [level_num, 1], calc_mean, NaN);
else
    mean_AdpCoh = NaN(level_num, 1);
end
if has_modified
    mean_ModifiedCoh = accumarray(groups, double(Stack.ModifiedCoh(valid_mask)), [level_num, 1], calc_mean, NaN);
else
    mean_ModifiedCoh = NaN(level_num, 1);
end

T1 = table(levels, mean_TrueCoh, mean_BoxCohMax, mean_AdpCoh, mean_ModifiedCoh, mean_BestEst);
writetable(T1, 'fig1_cv_data.txt', 'Delimiter', '\t');

%% Fig 2
T2 = table(levels, mean_TrueCoh);
for j = 1:length(box_wins)
    tmp = double(Stack.(box_coh_fields{j})(valid_mask));
    T2.(sprintf('BoxCoh_%d', box_wins(j))) = accumarray(groups, tmp(:), [level_num, 1], calc_mean, NaN);
end
writetable(T2, 'fig2_cv_data.txt', 'Delimiter', '\t');

%% Fig 3
ratio_BoxAdp = mean_BoxCohMax ./ mean_AdpCoh;
diff_BoxAdp = mean_BoxCohMax - mean_AdpCoh;
if has_modified
    ratio_BoxMod = mean_BoxCohMax ./ mean_ModifiedCoh;
    ratio_AdpMod = mean_AdpCoh ./ mean_ModifiedCoh;
    diff_BoxMod = mean_BoxCohMax - mean_ModifiedCoh;
    diff_AdpMod = mean_AdpCoh - mean_ModifiedCoh;
else
    ratio_BoxMod = NaN(level_num, 1); ratio_AdpMod = NaN(level_num, 1);
    diff_BoxMod = NaN(level_num, 1); diff_AdpMod = NaN(level_num, 1);
end
T3 = table(levels, ratio_BoxAdp, ratio_BoxMod, ratio_AdpMod, diff_BoxAdp, diff_BoxMod, diff_AdpMod);
writetable(T3, 'fig3_cv_data.txt', 'Delimiter', '\t');

%% Fig 4, 5, 6 setup
unique_wins = unique(Stack.Best_Win(:));
unique_wins = unique_wins(~isnan(unique_wins) & unique_wins > 0);
unique_wins = sort(unique_wins);
num_wins = length(unique_wins);

%% Fig 4
win_labels = []; mean_cv_arr = []; q1_arr = []; q2_arr = []; q3_arr = []; min_arr = []; max_arr = [];
for i = 1:num_wins
    win_val = unique_wins(i);
    idx = (Stack.Best_Win == win_val);
    cv_data = Stack.CV_level(idx);
    cv_data(isnan(cv_data)) = []; % Remove NaNs before percentile calc
    
    win_labels = [win_labels; win_val];
    if isempty(cv_data)
        min_arr = [min_arr; NaN]; q1_arr = [q1_arr; NaN];
        q2_arr = [q2_arr; NaN]; q3_arr = [q3_arr; NaN];
        max_arr = [max_arr; NaN]; mean_cv_arr = [mean_cv_arr; NaN];
    else
        % For boxplots, use 1.5*IQR rule for whiskers instead of absolute min/max
        q1_val = prctile(cv_data, 25);
        q2_val = median(cv_data);
        q3_val = prctile(cv_data, 75);
        iqr_val = q3_val - q1_val;
        
        lower_bound = q1_val - 1.5 * iqr_val;
        upper_bound = q3_val + 1.5 * iqr_val;
        
        w_min = min(cv_data(cv_data >= lower_bound));
        w_max = max(cv_data(cv_data <= upper_bound));
        
        if isempty(w_min), w_min = min(cv_data); end
        if isempty(w_max), w_max = max(cv_data); end

        min_arr = [min_arr; w_min];
        q1_arr = [q1_arr; q1_val];
        q2_arr = [q2_arr; q2_val];
        q3_arr = [q3_arr; q3_val];
        max_arr = [max_arr; w_max];
        mean_cv_arr = [mean_cv_arr; mean(cv_data)];
    end
end
T4 = table(win_labels, min_arr, q1_arr, q2_arr, q3_arr, max_arr, mean_cv_arr);
writetable(T4, 'fig4_cv_boxplot.txt', 'Delimiter', '\t');

%% Fig 5
win_percent_matrix_cv = zeros(level_num, num_wins);
for i = 1:level_num
    level_idx = (Stack.CV_level == i);
    for j = 1:num_wins
        win_val = unique_wins(j);
        win_percent_matrix_cv(i, j) = sum(Stack.Best_Win(level_idx) == win_val);
    end
end
row_sums = sum(win_percent_matrix_cv, 2);
win_percent_matrix_cv = win_percent_matrix_cv ./ row_sums * 100;
win_percent_matrix_cv(isnan(win_percent_matrix_cv)) = 0;

T5 = table(levels);
for j = 1:num_wins
    T5.(sprintf('Win_%d', unique_wins(j))) = win_percent_matrix_cv(:, j);
end
writetable(T5, 'fig5_cv_stacked.txt', 'Delimiter', '\t');

%% Fig 6
mean_win_cv = zeros(level_num, 1);
for i = 1:level_num
    level_idx = (Stack.CV_level == i);
    best_wins = Stack.Best_Win(level_idx);
    best_wins = best_wins(~isnan(best_wins) & best_wins > 0);
    if ~isempty(best_wins)
        mean_win_cv(i) = mean(best_wins);
    else
        mean_win_cv(i) = NaN;
    end
end
T6 = table(levels, mean_win_cv);
writetable(T6, 'fig6_cv_mean_win.txt', 'Delimiter', '\t');


disp('Data exported to txt successfully (CV). MATLAB processing complete.');
