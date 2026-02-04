function outfile = SingleWrite(image,filename,imgpath,suffixname,bkformat,machinefmt)
%SingleWrite Write a single image using fwritebkj helper.
%
%   outfile = SingleWrite(image,filename,imgpath,suffixname,bkformat,machinefmt)
%
%   Inputs:
%   - image:      2-D matrix to be written.
%   - filename:   Base name (numeric date, [date1 date2], or string without extension).
%   - imgpath:    Output directory.
%   - suffixname: Extension without dot, e.g. 'mli', 'diff'.
%   - bkformat:   Data format (see fwritebkj), e.g. 'float32', 'cpxfloat32'.
%   - machinefmt: Byte order, e.g. 'b' (big-endian). Defaults match ImgWrite.
%
%   Output:
%   - outfile:    Full path of the written file.
%
%   Example:
%   SingleWrite(data,20100311,'/path/to/MLI','mli','float32','b');

if nargin < 6
    machinefmt = 'b';
end

if nargin < 5
    bkformat = 'float32';
end

if nargin < 4
    help SingleWrite
    return;
end

% Normalize output directory (drop trailing separator).
if ~isempty(imgpath) && strcmp(imgpath(end),filesep)
    imgpath = imgpath(1:end-1);
end

% Create directory when missing.
if ~exist(imgpath,'dir')
    mkdir(imgpath);
else
    k = strfind(imgpath,filesep);
    if ~isempty(k)
        disp([imgpath(k(end)+1:end), ' directory already exists...']);
    end
end

% Build base filename.
if isnumeric(filename)
    filename = filename(:).';
    if numel(filename) == 1
        file_base = num2str(filename);
    elseif numel(filename) == 2
        file_base = [num2str(filename(1)),'_',num2str(filename(2))];
    else
        error('Numeric filename should be scalar or 1x2 for dual-date naming.');
    end
elseif isstring(filename) || ischar(filename)
    [~,file_base,~] = fileparts(char(filename));
else
    error('filename must be numeric or text.');
end

outfile = [imgpath,filesep,file_base,'.',suffixname];

tic;
fwritebkj(image,outfile,bkformat,machinefmt);
elapsed = toc;

fprintf('Writing Img to %s, time = %.0f sec\n', outfile, elapsed);

end
