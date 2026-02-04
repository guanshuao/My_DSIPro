clc;clear;close all;

curpath   = pwd;
mlipath   = [curpath,filesep,'ML',filesep,'MLI'];
diffpath  = [curpath,filesep,'ML',filesep,'DIFF'];

% Check if paths exist
if ~exist(mlipath, 'dir') || ~exist(diffpath, 'dir')
    error('MLI or DIFF directory not found. Please ensure you are in the correct directory.');
end

% Read parameter file to get nlines
tag_files = dir([mlipath,filesep, '*.par']);
if isempty(tag_files)
    error('No .par files found in %s', mlipath);
end
fname     = [mlipath,filesep,tag_files(1).name]; 
nlines    = readparam('azimuth_lines',fname); 

% Data input
% Reading intfstack to get dimensions and filenames
disp('Reading intfstack...');
intfstack = ImgRead(diffpath,'slc',nlines,'cpxfloat32');

[rows, cols, npages] = size(intfstack.datastack);
disp(['Data size: ', num2str(rows), ' x ', num2str(cols), ' x ', num2str(npages)]);

% Generate ones matrix
% Each element is 1+0i (complex float 32)
disp('Generating ones stack...');
ones_stack = complex(ones(rows, cols, npages, 'single'), 0);

% Prepare filenames
% The user requested: yyyymmdd_yyyymmdd.rslc, where both dates are the reference (master) date.
% intfstack.filename is usually N x 2 [Master, Slave] for interferograms.
% We will use [Master, Master] for the new filenames.
if size(intfstack.filename, 2) >= 1
    ref_dates = intfstack.filename(:, 1);
    ones_filenames = [ref_dates, ref_dates];
else
    error('Unexpected filename format in intfstack.');
end

% Output path
outpath = [curpath, filesep, 'ONES'];
if ~exist(outpath, 'dir')
    mkdir(outpath);
end

% Write the data
disp(['Writing output to ', outpath, '...']);
ImgWrite(ones_stack, ones_filenames, outpath, 'rslc', 'cpxfloat32', 'b');

disp('Done.');