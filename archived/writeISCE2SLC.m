function writeISCE2SLC(slcData, slcFilename, varargin)
%WRITEISCE2SLC 将复数矩阵写入ISCE2 SLC格式文件
%
%   writeISCE2SLC(slcData, slcFilename)
%   writeISCE2SLC(slcData, slcFilename, 'Name', Value, ...)
%
%   输入参数:
%       slcData     - 复数矩阵 (length x width)
%       slcFilename - 输出SLC文件路径
%
%   可选Name-Value参数: 
%       'ByteOrder' - 字节序: 'l'小端(默认), 'b'大端
%       'WriteXML'  - 是否生成XML文件 (默认true)
%       'WriteVRT'  - 是否生成VRT文件 (默认true)
%
%   示例:
%       % 写入SLC文件及元数据
%       writeISCE2SLC(myComplexData, 'output.slc');
%
%   See also: readISCE2SLC

    p = inputParser;
    addRequired(p, 'slcData', @isnumeric);
    addRequired(p, 'slcFilename', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ByteOrder', 'l', @ischar);
    addParameter(p, 'WriteXML', true, @islogical);
    addParameter(p, 'WriteVRT', true, @islogical);
    parse(p, slcData, slcFilename, varargin{:});
    
    opts = p. Results;
    slcFilename = char(opts.slcFilename);
    
    % 确定字节序
    if ismember(lower(opts.ByteOrder), {'l', 'little'})
        machineFormat = 'ieee-le';
        byteOrderStr = 'l';
    else
        machineFormat = 'ieee-be';
        byteOrderStr = 'b';
    end
    
    [length_img, width] = size(slcData);
    
    % 确保数据为single精度复数
    slcData = single(slcData);
    
    % 打开文件写入
    fid = fopen(slcFilename, 'wb', machineFormat);
    if fid == -1
        error('无法创建文件: %s', slcFilename);
    end
    
    cleanupObj = onCleanup(@() fclose(fid));
    
    % 交错写入实部和虚部
    realPart = real(slcData)';
    imagPart = imag(slcData)';
    
    % 交错排列
    interleavedData = zeros(2, width * length_img, 'single');
    interleavedData(1, :) = realPart(: );
    interleavedData(2, :) = imagPart(:);
    
    fwrite(fid, interleavedData, 'float32');
    
    fprintf('成功写入SLC文件: %s\n', slcFilename);
    fprintf('  尺寸: %d 行 x %d 列\n', length_img, width);
    
    % 生成XML文件
    if opts.WriteXML
        xmlFilename = [slcFilename, '.xml'];
        writeISCE2XML(xmlFilename, width, length_img, 'CFLOAT', byteOrderStr);
        fprintf('  已生成XML文件: %s\n', xmlFilename);
    end
    
    % 生成VRT文件
    if opts. WriteVRT
        vrtFilename = [slcFilename, '.vrt'];
        writeISCE2VRT(vrtFilename, slcFilename, width, length_img);
        fprintf('  已生成VRT文件: %s\n', vrtFilename);
    end
end


function writeISCE2XML(xmlFilename, width, length_img, dataType, byteOrder)
%WRITEISCE2XML 生成ISCE2格式的XML元数据文件
    
    fid = fopen(xmlFilename, 'w');
    if fid == -1
        warning('无法创建XML文件: %s', xmlFilename);
        return;
    end
    
    fprintf(fid, '<imageFile>\n');
    fprintf(fid, '    <property name="WIDTH">\n');
    fprintf(fid, '        <value>%d</value>\n', width);
    fprintf(fid, '    </property>\n');
    fprintf(fid, '    <property name="LENGTH">\n');
    fprintf(fid, '        <value>%d</value>\n', length_img);
    fprintf(fid, '    </property>\n');
    fprintf(fid, '    <property name="NUMBER_BANDS">\n');
    fprintf(fid, '        <value>1</value>\n');
    fprintf(fid, '    </property>\n');
    fprintf(fid, '    <property name="DATA_TYPE">\n');
    fprintf(fid, '        <value>%s</value>\n', dataType);
    fprintf(fid, '    </property>\n');
    fprintf(fid, '    <property name="BYTE_ORDER">\n');
    fprintf(fid, '        <value>%s</value>\n', byteOrder);
    fprintf(fid, '    </property>\n');
    fprintf(fid, '    <property name="SCHEME">\n');
    fprintf(fid, '        <value>BIP</value>\n');
    fprintf(fid, '    </property>\n');
    fprintf(fid, '    <property name="IMAGE_TYPE">\n');
    fprintf(fid, '        <value>slc</value>\n');
    fprintf(fid, '    </property>\n');
    fprintf(fid, '</imageFile>\n');
    
    fclose(fid);
end


function writeISCE2VRT(vrtFilename, slcFilename, width, length_img)
%WRITEISCE2VRT 生成GDAL VRT格式文件
    
    [~, fname, ext] = fileparts(slcFilename);
    sourceFilename = [fname, ext];
    
    fid = fopen(vrtFilename, 'w');
    if fid == -1
        warning('无法创建VRT文件:  %s', vrtFilename);
        return;
    end
    
    fprintf(fid, '<VRTDataset rasterXSize="%d" rasterYSize="%d">\n', width, length_img);
    fprintf(fid, '    <VRTRasterBand dataType="CFloat32" band="1" subClass="VRTRawRasterBand">\n');
    fprintf(fid, '        <SourceFilename relativeToVRT="1">%s</SourceFilename>\n', sourceFilename);
    fprintf(fid, '        <ImageOffset>0</ImageOffset>\n');
    fprintf(fid, '        <PixelOffset>8</PixelOffset>\n');
    fprintf(fid, '        <LineOffset>%d</LineOffset>\n', width * 8);
    fprintf(fid, '        <ByteOrder>LSB</ByteOrder>\n');
    fprintf(fid, '    </VRTRasterBand>\n');
    fprintf(fid, '</VRTDataset>\n');
    
    fclose(fid);
end