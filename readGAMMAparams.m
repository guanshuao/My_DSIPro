function ImgParms = readGAMMAparams(fname)
% This function is used to read parameter list from a image file
% Using GAMMA as a example
%
%   Mi JIANG, Sun Yat-sen University,  

ImgParms = struct('created',date);

ImgParms.nlines  = readparam('azimuth_lines',fname);
ImgParms.nwidths = readparam('range_samples',fname);
ImgParms.rlks  = readparam('range_looks',fname);
ImgParms.azlks = readparam('azimuth_looks',fname);
ImgParms.rg_spacing = readparam('range_pixel_spacing',fname);
ImgParms.az_spacing = readparam('azimuth_pixel_spacing',fname);
ImgParms.inc = readparam('incidence_angle',fname);
