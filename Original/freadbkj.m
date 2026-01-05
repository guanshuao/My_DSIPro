function [data, count] = freadbkj(infile, lines, bkformat, machinefmt ,r0, rN, c0, cN)
% FREADBK  --  Read binary data file.
%
%   DATA = FREADBK pops up a menu to select a filename, asks for the number
%   of lines and format of the input file, and reads in the data (double).
%
%   DATA = FREADBK(FILENAME) fopens the specified file with read permission and
%   returns a row vector DATA.  Format is assumed to be float32 format.
%
%   [DATA, COUNT] = FREADBK(FILENAME) optionally returns the number of elements
%   successfully read.
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES) includes an optional NUMLINES
%   argument for the number of lines in the file.  Matrix DATA is returned of
%   size NUMLINES rows and the appropriate number of columns.
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES, BKFORMAT) reads in a file in
%   the format BKFORMAT.  Which is one of the formats specified in FREAD (see
%   help fread), or one of these with a prepended 'cpx' for complex data types.
%   (DATA must be stored major row order, pixel interleaved):
%	'cpxfloat32'	complex floating point, 32 bits.
%	'cpxint16'	complex signed 16 bit short integers, ...
%	'cpx...'	...
%      ('mph'           format: use 'cpxfloat32')
%      ('hgt'           format: use 'freadhgt' function)
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES, BKFORMAT, MACHINEFMT) additionally 
%   specifies the order for reading or writing bytes or bits in the file using the machinefmt argument
%   MACHINEFMT:
%   'n' or 'native'         Your system byte ordering (default)
%   'b' or 'ieee-be'        Big-endian ordering (GAMMA software, for example)
%   'l' or 'ieee-le'        Little-endian ordering
%   's' or 'ieee-be.l64'    Big-endian ordering, 64-bit long data type
%   'a' or 'ieee-le.l64'    Little-endian ordering, 64-bit long data type
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES, BKFORMAT, r0, rN, c0, cN) only
%   read in the part between: row r0 to row rN, and between column c0 to cN.
%   Indexes start at 1, stop at numlines/width. r0=0 all rows, c0=0 means all
%   columns.
%
%   Examples:
%     To read in a file with 2 channels (complex), stored pixel interleaved, ie,
%     (row1: RE1 - IM1 - RE2 - IM2 - RE3 - ...)
%     (row2: RE1 - IM1 - RE2 - IM2 - RE3 - ...)
%     (      ...                              )
%     in float32 format, which has a height of 100 lines, and a width of 250
%     complex pixels, use:
%       D = freadbk('filename',100,'cpxfloat32');
%
%     To crop this file while reading, between rows 1:10, columns 101:110, use:
%       D = freadbk('filename',100,'cpxfloat32','n',1, 10, 101, 110);
%
%   See also FWRITEBK, FOPEN, FREAD, FWRITE, LOAD, SAVE, FPRINTF,
%   FSEEK, FTELL, FREADHGT, FWRITEHGT
%

% $Revision: 1.12 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 04-Mar-2000

%   ======================================================================
%   04/2009 MJ add machinefmt
%   ======================================================================

%%% Handle input. (varargin?)
%if (nargin >  3) error('working on this... not yet ok, sorry.'); end;
false=0; true=1; 
complextype=false; readwholefile=false;%		defaults

if (nargin <  8) cN = 0; end;%				fall through (not...)
if (nargin <  7) c0 = 0; end;%				fall through
if (nargin <  6) rN = 0; end;%				fall through 
if (nargin <  5) r0 = 0; readwholefile=true; end;%	fall through
if (nargin <  4) machinefmt = 'n'; end;%	fall through
if (nargin <  3) bkformat='float32'; end;%		fall through
if (nargin <  2) lines=1; end;%				fall through
if (nargin <  1)
    [infile, inpath] = uigetfile('*', 'Select binary input float file', 0,0);
    infile   = [inpath,infile];
    lines    = input('Enter number of lines in file: ');
    bkformat = input('Format (enter between single quotes): ');% make a gui...
end;
%%% Check if whole file should be read regardless of r0,rNc0,cN.
if (c0==1 & r0==1 & rN==lines)
  bytesperelem = datatypesize(bkformat);%	in bytes
  filesize     = fsize(infile);%		in bytes
  filewidth    = filesize/(bytesperelem*lines);%ok for complex
  if (cN==filewidth)
    readwholefile=true;
  end;
end
% if (readwholefile==true)
%   disp(['% Reading whole file: ', infile]);
% end;

%%% Check bkformat for complex type: 'cpx*'
if (~ischar(bkformat)) error('FREADBK: bkformat must be string.'); end;
if (~ischar(infile))   error('FREADBK: infile must be string.'); end;
if (length(bkformat)==3)
  if (bkformat=='mph') 
    disp('changing mph format to cpxfloat32');
    bkformat='cpxfloat32';
  end
  if (bkformat=='hgt') 
    error('please use freadhgt for hgt format.');
  end
end
%%% complex types defined as prepended 'cpx'
if (length(bkformat)>6)
  if (bkformat(1:3)=='cpx')
    complextype = true;
    bkformat=bkformat(4:length(bkformat));
  end;
end;

%%% Read from file in column vector.
fid = fopen(infile,'r',machinefmt);
if (fid<0)%						try one more time
  [infile, inpath] = uigetfile('*', 'Select binary input file', 0,0);
  infile = [inpath,infile];
  fid    = fopen(infile,'r',machinefmt);
  if (fid<0) error(ferror(fid)); end;
end;

if (readwholefile==true)
  [data,count]=fread(fid,bkformat);%		read data in column vector
				%		count is number of elements, not bytes...
else
  bytesperelem = datatypesize(bkformat);%	in bytes
  filesize     = fsize(fid);%			in bytes
  filewidth    = filesize/(bytesperelem*lines);%ok for complex
  if (r0 > rN) error('r0 > rN'); end;
  if (c0 > cN) error('c0 > cN'); end;
  if (c0==0) c0=1; cN=filewidth; end;
  if (r0==0) r0=1; rN=lines; end;
  if (r0 < 1)            error('r0 < 1'); end;
  if (rN > lines)        error('rN > lines'); end;
  if (c0 < 1)            error('c0 < 1'); end;
  if (cN > filewidth)    error('cN > width'); end;
  if (~isint(filewidth)) error('numlines file seems to be wrong(?)'); end;

  data    = [];
  count   = 0;
  offset  = bytesperelem*(c0-1);%		in bytes
  mywidth = cN-c0+1;%				number of elems 2b read
  lines   = rN-r0+1;%				number of lines in matrix
  if (complextype==true)
    mywidth=mywidth*2;%				correction for pix interleaved complex
    offset =offset*2;%				correction for pix interleaved complex
  end;
  for ii=r0-1:rN-1
    start   = ii*filewidth*bytesperelem+offset;
    status  = fseek(fid,start,'bof');%			beginning of file
    if (status ~= 0)
      disp(['r0 rN c0 cN: ',num2str(r0),' ',num2str(rN),...
			' ',num2str(c0),' ',num2str(cN)]);
      disp(['ii offset start: ',...
	     num2str(ii),' ',num2str(offset),' ',num2str(start)]);
      error(ferror(fid));
    end;
    [d cnt] = fread(fid,mywidth,bkformat);%		column vector
    count   = count+cnt;
    data    = [data;d];%				augment column vector
  end;
end;
fclose(fid);


if (complextype==true)
  % should be possible by something like, but for speed makes no big diff.
  %data = reshape(data,width,2*lines);%		transpose band interleaved if complex
  %data = complex(data(:,1:lines),data(:,lines+1:2*lines)).';
  data = complex(data(1:2:count),data(2:2:count));
  count=count/2;%					correction for complex types
end;

width = count/lines;%
data  = reshape(data,width,lines).';



%%% Export to type if requested 
% data = int8(data) ...


%%% EOF

function b = datatypesize(c)
% B=DATATYPESIZE(C)
% return number of bytes B based on string C.
% extended with cpx types for complex, see also freadbk.
% -1 is returned for unknown type
%
% help fread
%        MATLAB    C or Fortran     Description
%        'uchar'   'unsigned char'  unsigned character,  8 bits.
%        'schar'   'signed char'    signed character,  8 bits.
%        'int8'    'integer*1'      integer, 8 bits.
%        'int16'   'integer*2'      integer, 16 bits.
%        'int32'   'integer*4'      integer, 32 bits.
%        'int64'   'integer*8'      integer, 64 bits.
%        'uint8'   'integer*1'      unsigned integer, 8 bits.
%        'uint16'  'integer*2'      unsigned integer, 16 bits.
%        'uint32'  'integer*4'      unsigned integer, 32 bits.
%        'uint64'  'integer*8'      unsigned integer, 64 bits.
%        'single'  'real*4'         floating point, 32 bits.
%        'float32' 'real*4'         floating point, 32 bits.
%        'double'  'real*8'         floating point, 64 bits.
%        'float64' 'real*8'         floating point, 64 bits.
% 
%    The following platform dependent formats are also supported but
%    they are not guaranteed to be the same size on all platforms.
%        'char'    'char*1'         character,  8 bits (signed or unsigned).
%        'short'   'short'          integer,  16 bits.
%        'int'     'int'            integer,  32 bits.
%        'long'    'long'           integer,  32 or 64 bits.
%        'ushort'  'unsigned short' unsigned integer,  16 bits.
%        'uint'    'unsigned int'   unsigned integer,  32 bits.
%        'ulong'   'unsigned long'  unsigned integer,  32 bits or 64 bits.
%        'float'   'float'          floating point, 32 bits.
%
%
% See also FREAD, FREADBK, FSEEK, FTELL
%


% $Revision: 1.4 $  $Date: 2001/09/28 14:24:31 $
% Bert Kampes, 4/3/00

% %%% see fread
% if (nargin~=1 | ~ischar(c)) 
%     helphelp; 
%     return; 
% end


%%% Check complex type request.
cpx=1; % 2 for complex type
if (length(c)>3)
  if (c(1:3)=='cpx')
    c   = c(4:length(c));
    cpx = 2;
  end
end

switch (lower(c))
  case {'uchar','unsigned char','int8','integer*1'}, 	b=1;
  case {'int16','integer*2'}, 				b=2;
  case {'int32','integer*4'}, 				b=4;
  case {'int64','integer*8'}, 				b=8;
  case {'uint8','integer*1'}, 				b=1;
  case {'uint16','integer*2'}, 				b=2;
  case {'uint32','integer*4'}, 				b=4;
  case {'uint64','integer*8'}, 				b=8;
  case {'single','real*4'},				b=4;
  case {'float32','real*4'},				b=4;
  case {'double','real*8'},				b=8;
  case {'float64','real*8'},				b=8;
  case {'char','char*1'},				b=1;
  case {'short'},					b=2;
  case {'int'},						b=4;
  case {'long'},					b=4;
  case {'ushort'},					b=2;
  case {'uint'},					b=4;
  case {'ulong'},					b=4;
  case {'float'},					b=4;
  otherwise warning('Unknown size for type ...'), 	b=-1;
end

%%% Correct for complex type
if ( b>0 )
  b=b*cpx;
end

function numbytes = fsize(filename)
% FSIZE  return size of file in bytes.
%   b = FSIZE(FILENAME) return size in bytes of specified file.
%   b = FSIZE(FID)      return size in bytes of specified fid (file).
%   -1 is returned if specified file does not exist. 
%   (I found sometimes problems when file already was opened?)
%

% $Revision: 1.5 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 08-Mar-2000

% check if input is filename or fid
if (ischar(filename)) 
  tmpfid=fopen(filename,'r');
else
  tmpfid=filename;
  if (tmpfid<0)
    numbytes=-1;
  else
    oldpos=ftell(tmpfid);
  end
end

if (tmpfid<0)
  numbytes=-1;
else
  status=fseek(tmpfid,0,'eof');
  if (status==-1) error(ferror(tmpfid)); end;
  numbytes=ftell(tmpfid);
  if (numbytes==-1) error(ferror(tmpfid)); end;
  if (ischar(filename))
    fclose(tmpfid);%			close file if opened
  else
    fseek(tmpfid,oldpos,'bof');%	or reset file pointer
  end;
end;

%%% EOF
function res = isint(in)
% ISINT  True for integers.
%   ISINT(IN) returns 1 for those elements that are integers,
%   and 0 otherwise.  Valid for numbers in [-2147483648, 2147483647].
%   Mainly used for input checking in functions.
%
%   To test if any/all element of an array is an integer, use any or all:
%     if(any(isint(1:4))) disp('at least one element is an integer'); end;
%     if(all(isint(1:4))) disp('all elements are integers'); end;
%
%   Example:
%     a = [1 2 3 4 5.5 8.34];
%     isint(a)%  returns [1 1 1 1 0 0]
%
%   See also INT32
%

% $Revision: 1.3 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 15/03/00

%res = 0;
%if (int32(in)==in) res = 1; end;

%// BK 07-Aug-2001
%%%nico sneeuw: res=(rem(in,1));
res=(int32(in)==in);
