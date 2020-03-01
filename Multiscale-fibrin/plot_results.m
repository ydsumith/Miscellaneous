clc;
clear;

NX = 0;
NY = 0;

fid = fopen('navier_uu.data');
firstline = fgets(fid);
tempM = strsplit(firstline,'\t');
ROW = str2num(cell2mat(tempM(1,1)));
COL = str2num(cell2mat(tempM(1,2)));

data_uu = zeros(ROW,COL);
data_vv = zeros(ROW,COL);
data_w = zeros(ROW,COL);

tline = fgets(fid);
row_count = 1;
while ischar(tline)
    tempM = strsplit(tline,'\t');
    for i = 1:COL
      data_uu(row_count,i) = str2num(cell2mat(tempM(1,i)));
    end
    tline = fgets(fid);
    row_count = row_count + 1;
end
fclose(fid);

fid = fopen('navier_vv.data');
firstline = fgets(fid);
tline = fgets(fid);
row_count = 1;
while ischar(tline)
    tempM = strsplit(tline,'\t');
    for i = 1:COL
      data_vv(row_count,i) = str2num(cell2mat(tempM(1,i)));
    end
    tline = fgets(fid);
    row_count = row_count + 1;
end
fclose(fid);

fid = fopen('navier_w.data');
firstline = fgets(fid);
tline = fgets(fid);
row_count = 1;
while ischar(tline)
    tempM = strsplit(tline,'\t');
    for i = 1:COL
      data_w(row_count,i) = str2num(cell2mat(tempM(1,i)));
    end
    tline = fgets(fid);
    row_count = row_count + 1;
end
fclose(fid);

hold off;
axis equal;
axis([1 ROW 1 COL]);
quiver(flipud(rot90(data_uu)),  flipud(rot90(data_vv)), 0.4);
hold on;
contour(flipud(rot90(data_w)) ,100);

disp('done');
