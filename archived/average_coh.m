clear;clc;

matpath = '/titan/guanshuao/Kumamoto/Data/AdpCoh_BWS_21';

% 获取该路径下的所有.mat文件
% 读取后检查矩阵的尺寸是否一致
% 若一致，将这些矩阵逐一读取，对其幅度进行平均（因为矩阵中的元素均为复数），保存为average.mat
files = dir(fullfile(matpath, '*.mat'));
numFiles = length(files);
if numFiles == 0
    error('No .mat files found in the specified directory.');
end

% 读取第一个文件以获取矩阵尺寸
firstData = load(fullfile(matpath, files(1).name));
fieldNames = fieldnames(firstData);
matrix = firstData.(fieldNames{1});
[m, n] = size(matrix);
sumMatrix = zeros(m, n);
count = 0;
for i = 1:numFiles
    data = load(fullfile(matpath, files(i).name));
    disp(['Processing file: ', files(i).name]);
    fieldNames = fieldnames(data);
    currentMatrix = data.(fieldNames{1});
    [curM, curN] = size(currentMatrix);
    if curM == m && curN == n
        sumMatrix = sumMatrix + abs(currentMatrix);
        count = count + 1;
        disp(count);
    else
        fprintf('Skipping file %s due to size mismatch.\n', files(i).name);
    end
end
if count == 0
    error('No matrices with matching dimensions were found.');
end
AdpCoh_21 = sumMatrix / count;
save(fullfile(matpath, 'AdpCoh_21.mat'), 'AdpCoh_21', '-v7.3');

disp('number of matrices averaged: ');
disp(count);
disp('Average matrix saved as AdpCoh_21.mat');