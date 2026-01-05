function [Footprint]=SHP_FootPrint(imagename,SHP,PointIdx,imagesize_radius)
%This function is used to show footprint of SHPs for a given reference pixel
%
%   Usage:
%       [Footprint]=SHP_FootPrint(imagename,SHP,PointIdx,imagesize_radius)
%   
%
%   Inputs:
%   - imagename:        Basemap path with ras format
%   - SHP:              See script "SHP_SelPoint.m" for details
%   - PointIdx:         [row col] of basemap
%   - imagesize_radius: The radius of the background
%
%   Outputs:
%   - Footprint:        A point coordinate [row column]
%
%   Examples:
%   Select a reference point manually from original basemap, with a radius 25 pixels, use:
%    [Footprint]=SHP_FootPrint('/home/user/INSAR/COHEST/MLI/mli_ave.ras',SHP);
%   Compare the result with that from the other SHP set, e.g., BWS, use:
%    SHP_FootPrint('/home/user/INSAR/COHEST/MLI/mli_ave.ras',BWSSHP,Footprint);   
%
%
%   
%   This toolbox can be used only for research purposes, you should cite 
%   the aforementioned papers in any resulting publication.
%
%   Mi JIANG, Sun Yat-sen Univeristy,

% read base map
[mean_amp,amp_color]=imread(imagename);

if nargin < 4
    imagesize_radius=[25 25]; %[row col]
end

if nargin < 3
    disp('select a point from basemap...')
    clf;
    figure(1);image(mean_amp);colormap(amp_color);axis image
    PointIdx=fliplr(round(ginput(1)));    
end    

if nargin < 2
    help SHP_FootPrint
    return;
end
CalWin=SHP.CalWin;
[nlines,nwidths]=size(mean_amp);
IND = sub2ind([nlines,nwidths],PointIdx(1),PointIdx(2));
SHPS= reshape(SHP.PixelInd(:,IND), [CalWin(1),CalWin(2)]);
%Footprint coordinate in basemap
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;
[X,Y] = meshgrid(PointIdx(2)-RadiusCol:PointIdx(2)+RadiusCol,PointIdx(1)-RadiusRow:PointIdx(1)+RadiusRow);
footprintColREC = X(:);%footprint for REC
footprintRowREC = Y(:);
X=X(SHPS);
Y=Y(SHPS);
%remove footprint at image broundary
negvalue = X<=0|Y<=0|X>nwidths|Y>nlines;
X=X(~negvalue);
Y=Y(~negvalue);
%cut window
ShowRow = imagesize_radius(1);%image size
ShowCol = imagesize_radius(2);
CropWin = [PointIdx(1)-ShowRow,PointIdx(1)+ShowRow,PointIdx(2)-ShowCol,PointIdx(2)+ShowCol];
if CropWin(1) < 1
        CropWin(1) = 1;
end
if CropWin(2) > nlines
        CropWin(2) = nlines;
end
if CropWin(3) < 1;
        CropWin(3) = 1;
end
if CropWin(4) > nwidths
        CropWin(4) = nwidths;
end
Footprint = [PointIdx(1),PointIdx(2)];
%output
figure; colormap(amp_color);
subplot(1,2,2);image(CropWin(3):CropWin(4),CropWin(1):CropWin(2),mean_amp(CropWin(1):CropWin(2),CropWin(3):CropWin(4)));
hold on;plot(X,Y,'g.','MarkerSize',10);plot(PointIdx(2),PointIdx(1),'r.','MarkerSize',10);title('SHP'); axis image
subplot(1,2,1);image(CropWin(3):CropWin(4),CropWin(1):CropWin(2),mean_amp(CropWin(1):CropWin(2),CropWin(3):CropWin(4)));
hold on;plot(footprintColREC,footprintRowREC,'g.','MarkerSize',10);
plot(X,Y,'b.','MarkerSize',10); plot(PointIdx(2),PointIdx(1),'r.','MarkerSize',10);title('BoxCar');axis image
