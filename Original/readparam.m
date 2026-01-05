function value = readparam(param,fname)
%read specific parameter from parameter file, e.g., par file in GAMMA
%   Mi JIANG, Sun Yat-sen University,
if nargin < 1
    help readparam
    return;
end

fileID=fopen(fname);
params=textscan(fileID,'%s');
fclose(fileID);
params=params{1};

idx = strmatch(param,params);

if isempty(idx)
    error('empty value ...');
end

value = str2double(params(idx+1));

if isnan(value)
    error('the value does not exist...');
end

