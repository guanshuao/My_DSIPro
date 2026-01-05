function mlishow(mli,sc,exp)
%this function is used to show intensity image
%mli: intensity image
%sc:  display scale factor (default=1.)
%exp: display exponent (default=.35)

%   Mi JIANG, Sun Yat-sen Univeristy,
if nargin < 3
    exp=1;
end

if nargin < 2
    sc=.35;
end

mli=abs(mli);
mli(isnan(mli))=0;
P = mli.^exp;
nv = numel(nonzeros(P));
P = sum(P(:));

if P==0
   av=1;
else
   av=P/nv; 
end

scale = sc*120/av;
IMG = mli.^exp*scale;
IMG(IMG>255) = 255;
IMG(IMG<1&IMG~=0)=1;
imagesc(IMG);colormap gray