%
%Download and Install the Image Processing Toolbox

clc;clear; close all;

stack_prefix='10.0mM zstack2_';

fibrin_fiber_overlap_agg = [];
long_aggr = [];
perc_porous = [];

for istack=1:82
fn = strcat(strcat(stack_prefix,int2str(istack)),'.tif');
I = imread(fn);
BW1=im2bw(I,graythresh(I));
BW2=imadjust(mat2gray(BW1));
BW3 = fibermetric(BW2, 30, 'StructureSensitivity', 15);
BW = BW3 > 0.00002;
figure, imshowpair(I,BW,'montage'), hold on;

[H,T,R] = hough(BW);
P  = houghpeaks(H,400,'threshold',ceil(0.4*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
%plot(x,y,'s','Color','blue');

%lines = houghlines(BW,T,R,P,'FillGap',3,'MinLength',35);
lines = houghlines(BW,T,R,P,'FillGap',15,'MinLength',40);
linelen = [];
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');

   % Plot beginnings and ends of lines
   %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','blue');
   
   linelen(k) = norm(lines(k).point1 - lines(k).point2);
end

% fibrin fiber overlap aggregation
% calc number intersections with other lines, average
numlineintersect = [];
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   numlineintersect(k) = 0;
   for m = 1:length(lines)
      xy_m = [lines(m).point1; lines(m).point2];
      % https://www.mathworks.com/matlabcentral/fileexchange/27205-fast-line-segment-intersection
      out = lineSegmentIntersect([xy(1,1) xy(1,2) xy(2,1) xy(2,2)],[xy_m(1,1) xy_m(1,2) xy_m(2,1) xy_m(2,2)]);
      % https://www.mathworks.com/matlabcentral/fileexchange/56835-lineintersection
      %[E, lambda, gamma, isConvex] = lineIntersection(xy(:,1),xy(:,2),xy_m(:,1),xy_m(:,2));
      %if not(isnan(E))
      if (out.intMatrixX>0 && out.intMatrixY>0)
          numlineintersect(k) = numlineintersect(k)+1;
      end
   end
end

% fibrin fiber overlap aggregation, average lineintersect
fibrin_fiber_overlap_agg(istack) = mean(numlineintersect)

% elongation, logitudinal aggregation = average length of line segments
long_aggr(istack) = mean(linelen)

% percent porosity
[s1,s2]=size(BW);
tarea=(s1*s2);
perc_porous(istack) = (tarea-bwarea(BW))/tarea

%imwrite(label2rgb(Pr_L),'Output.png')

fibrin_fiber_overlap_agg_mean = mean(fibrin_fiber_overlap_agg)
long_aggr_mean = mean(long_aggr)
perc_porous_mean = mean(perc_porous)

end;

%0.0mM zstack1_
%istack                          1           2       3       4           5       6           7           8
%fibrin_fiber_overlap_agg =   11.3987   12.2537   12.1214   11.9325   11.2103   11.7608   12.6884   12.9379
%long_aggr =                  72.2576   74.2070   74.4363   73.4357   72.9546   75.0652   74.9657   75.4506
%perc_porous =                 0.6854    0.6777    0.6759    0.6815    0.6782    0.6596    0.6516    0.6542