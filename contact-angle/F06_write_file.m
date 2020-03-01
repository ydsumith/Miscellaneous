function F06_write_file(filename,xx,yy,zz,xy_COUNT,Low_L,Low_H,dx,dz,flag)
fid = fopen(filename,'w');
fprintf(fid, '%d\n', xy_COUNT); %+NATOMS
fprintf(fid, 'xyz-file-argon\n');
if yy ==0
    for navi=1:xy_COUNT
        XSP = sprintf('%6s%10.5f%10.5f%10.5f\n',flag,((xx(navi)*dx)+Low_L)*10,0.001*10,((zz(navi)*dz)+Low_H)*10);
        fprintf(fid,XSP);
    end
else
    for navi=1:xy_COUNT
        XSP = sprintf('%6s%10.5f%10.5f%10.5f\n',flag,((xx(navi)*dx)+Low_L)*10,yy(navi)*dz*10,((zz(navi)*dz)+Low_H)*10);
        fprintf(fid,XSP);
    end
end
fclose(fid);
end