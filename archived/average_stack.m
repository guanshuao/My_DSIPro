clear; clc;

% Configuration
baseDir = '/titan/guanshuao/Kumamoto/process/dates_resampled';
subset = [201, 25200, 14001, 23000];
innerPath = fullfile('f1_2960', 's1');

% Find date folders (logic from readISCE2SLCStack.m)
if ~isfolder(baseDir)
    error('Base directory does not exist: %s', baseDir);
end

entries = dir(baseDir);
isDateDir = [entries.isdir];
names = {entries(isDateDir).name};
names = names(~ismember(names, {'.', '..'}));
% Filter for yymmdd format
dateNames = names(~cellfun(@isempty, regexp(names, '^\d{6}$', 'once')));

if isempty(dateNames)
    error('No date folders found under %s', baseDir);
end

% Construct SLC file paths
slcFiles = {};
for i = 1:numel(dateNames)
    d = dateNames{i};
    slcPath = fullfile(baseDir, d, innerPath, [d, '.slc']);
    if exist(slcPath, 'file')
        slcFiles{end+1} = slcPath; %#ok<AGROW>
    else
        warning('Skip %s: SLC not found at %s', d, slcPath);
    end
end

if isempty(slcFiles)
    error('No SLC files were found.');
end

fprintf('Found %d SLC files.\n', numel(slcFiles));

% Initialize accumulators
% Calculate size from subset
nRows = subset(2) - subset(1) + 1;
nCols = subset(4) - subset(3) + 1;

sum_amp = zeros(nRows, nCols, 'double');
sum_intensity = zeros(nRows, nCols, 'double');
nFiles = numel(slcFiles);

% Process files
for k = 1:nFiles
    slcPath = slcFiles{k};
    fprintf('Processing %d / %d : %s\n', k, nFiles, slcPath);
    
    % Read SLC
    try
        [slcData, ~] = readISCE2SLC(slcPath, 'Subset', subset);
    catch ME
        warning('Failed to read %s: %s. Skipping.', slcPath, ME.message);
        % If read fails, we should probably not count this file in the average
        % But since we pre-calculated nFiles, we need to adjust or handle it.
        % For simplicity, let's assume read success or error out. 
        % But to be robust, let's just fill with zeros or handle it.
        % Actually, if read fails, slcData is not assigned.
        % Let's just error out if read fails, as missing data in a stack is usually bad.
        rethrow(ME);
    end
    
    % Handle NaNs (treat as 0)
    slcData(isnan(slcData)) = 0;
    
    % Calculate amplitude
    amp = abs(slcData);
    
    % Accumulate amplitude
    sum_amp = sum_amp + double(amp);
    
    % Accumulate intensity (amplitude squared)
    sum_intensity = sum_intensity + double(amp.^2);
    
    % Clear temporary variables to save memory
    clear slcData amp;
end

% Calculate averages
avg_amp = sum_amp / nFiles;
avg_intensity = sum_intensity / nFiles;

% Save results
fprintf('Saving results...\n');
save('avg_amp.mat', 'avg_amp', '-v7.3');
save('avg_intensity.mat', 'avg_intensity', '-v7.3');

fprintf('Done. Processed %d files.\n', nFiles);

