function [] = write_mol( data , TIMESTEP, NATOMS, BOX, gap,dump_file)

fid = fopen(dump_file,'a');

XSP = sprintf('ITEM: TIMESTEP\n');
fprintf(fid,XSP);

XSP = sprintf('%d\n', TIMESTEP);
fprintf(fid,XSP);

XSP = sprintf('ITEM: NUMBER OF ATOMS\n');
fprintf(fid,XSP);

XSP = sprintf('%d\n', NATOMS);
fprintf(fid,XSP);

XSP = sprintf('ITEM: BOX BOUNDS pp pp pp\n');
fprintf(fid,XSP);

xlo = min(data(:,2)) - gap;
xhi = max(data(:,2)) + gap;

XSP = sprintf('%f %f\n', xlo, xhi);
fprintf(fid,XSP);

XSP = sprintf('%f %f\n', BOX(2,1), BOX(2,2));
fprintf(fid,XSP);

XSP = sprintf('%f %f\n', BOX(3,1), BOX(3,2));
fprintf(fid,XSP);

XSP = sprintf('ITEM: ATOMS id mol x y z\n');
fprintf(fid,XSP);

for i=1:NATOMS
    % id mol x y z
    XSP = sprintf('%d %d %.2f %.2f %.2f\n', i, round(data(i,1)), data(i,2), data(i,3), data(i,4));
    fprintf(fid,XSP);
end

fclose(fid);
end

