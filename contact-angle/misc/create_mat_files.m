clear;
clc;
%load('data\water_thesis.mat');
[filename, pathname] = uigetfile({'*.xyz'},'File Selector');
filename = strcat(pathname,filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (ispc) %# Windows
    totl_lines_in_xyz = str2num( perl('readline.pl', filename) );
else
    error('...');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename);
tline = fgetl(fid);
NATOMS = ceil(str2double(tline));
NDATA_SETS = totl_lines_in_xyz/(NATOMS + 2); % 2 is to account for extra two lines in xyz file
fid = fopen(filename);
kope = ceil(NDATA_SETS * NATOMS);
data = zeros(kope,3);
COUNTER = 1;
disp('Reading the xyz file, Please wait...\n');
%tline = fgetl(fid);
while ischar(tline)
    chekk = tline(1:4);
    if(strcmp('  Ar',chekk)== 1)
        if(COUNTER == 1)
            data(1,1) = str2double(tline(8:20))/10;
            data(1,2) = str2double(tline(24:35))/10;
            data(1,3) = str2double(tline(39:48))/10;
        else
            data(COUNTER,1) = str2double(tline(8:20))/10;
            data(COUNTER,2) = str2double(tline(24:35))/10;
            data(COUNTER,3) = str2double(tline(39:48))/10;
        end
        COUNTER = COUNTER +1;
    end
    tline = fgetl(fid);
    if mod(COUNTER,NATOMS) == 0
       fprintf('%3.2f completed\n',100*COUNTER/kope); 
    end
end
fclose(fid);
disp('File read into COORDINATES');
