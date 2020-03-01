clear;

load('singles\water-case5.mat'); % mat file contains only oxygen or center of mass of molecules
tic;
int_pol = 4; % don't change
disp('please wait...');
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% NATOMS = 4096;
Nw          = NATOMS; % no of atoms (read from mat file)
N           = NDATA_SETS; %no: of datasets  (read from mat file)
cellsize    = 0.1; %nm
Pt_thick    = 0.0; % single layer; use 0.46 nm for three layers
z_cutoff    = Pt_thick + 0.5; %nm ( 2*O-O)
cutoff      = ceil(z_cutoff/cellsize);
START_DATA  = 1;
END_DATA    = 1;
TOT_DATA = END_DATA- START_DATA +1;
mass        = 18.0; % mass of one molecule in amu (default 18 amu and 33 molecules/nm3)
d_theta = pi/5;
ang_theta = 0;
%ref_density = mass * 33.0;
create_file = 1; % 1 for yes, else no
INTPOL_TYPE = 2; % 0 for NGP, 1 for B-spline, 2 for Hardy
h_rc = 0.3;
drop_height = 7; % expected droplet height

h_rc2 = h_rc*h_rc;
A1 = 1.0/(h_rc2);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

[FACTORIAL_MAT] = F01_factorial(4);
RESULT = [];

%%
% loop through all data, fill x-z data
for i = START_DATA:END_DATA
    start_loc = (i-1)*Nw+1;
    end_loc = i * Nw;
    % copying x or y and z data for one time step
    X = data(start_loc:end_loc,1);
    Y = data(start_loc:end_loc,2);
    Z = data(start_loc:end_loc,3);
    L = max(X)+3;
    H = max(Z)+3;
    B = max(Y)+3;
    Low_L = 0;
    Low_H = 0;
    Norm_const = 2.0/( pi * h_rc2 * B);
    % find the grid dimensions
    NX = int32(L/cellsize);
    NZ = int32(H/cellsize);
    dx = cellsize;
    dz = cellsize;
    cell_vol = (dx*dz)*B;
    long_cell_density = mass/cell_vol;
    
    XZ_PLANE = zeros(NX,NZ);
    rot_count = 0;
    %%
    %--------------------------------------------
    %--
    % PASS #1 -- removing initial noise
    %--
    %--------------------------------------------
    MD_cutoffMAX = 1.2;
    MD_cutoffMIN = 0.0;
    [xy_COUNT,MODX,MODY,MODZ] = F02_find_MD(X,Y,Z,MD_cutoffMAX,MD_cutoffMIN,drop_height); % Mahalanobis distance
    X = MODX;
    Y = MODY;
    Z = MODZ;
    if create_file == 1
        F06_write_file('filtered.xyz',X,Y,Z,xy_COUNT,Low_L,Low_H,1,1,'OW');
    end
    %--------------------------------------------
    %--
    % PASS #2 -- smearing the molecules to a grid
    %--
    %--------------------------------------------
    for ang_theta = 0:d_theta:pi
        [MODX] = F05_rotate(X,Y,ang_theta);
        for j = 1:xy_COUNT
            ux = MODX(j)/dx; %dont round it
            uy = MODZ(j)/dz;
            switch INTPOL_TYPE
                case 0 % below three lines gives NGP assignment scheme
                    ROW = ceil(ux);
                    COL = ceil(uy);
                    if ROW > 0 && COL > 0 && ROW < NX && COL < NZ
                     XZ_PLANE (ROW,COL)= XZ_PLANE (ROW,COL)+1;
                    end
                case 1 % B-spline case
                    for k1 = 1: NX
                        u_x = ux-double(k1)+1.0;
                        if u_x < int_pol
                            if u_x > 0
                                wx = F03_get_M(u_x, 4,FACTORIAL_MAT);
                            else
                                continue
                            end
                        else
                            continue
                        end
                        for k2 = 1: NZ
                            u_y = uy-double(k2)+1.0;
                            if u_y < int_pol
                                if u_y > 0
                                    wy = F03_get_M(u_y, 4,FACTORIAL_MAT);
                                else
                                    continue
                                end
                            else
                                continue
                            end
                            XZ_PLANE(k1,k2) = XZ_PLANE(k1,k2) + wx * wy; % *1 is neglected
                        end
                    end
                case 2 % Hardy weight function
                    xi = MODX(j);
                    zi = MODZ(j);
                    for k1 = 1: NX
                        xp = double(k1)*dx;
                        rx = abs(xi - xp);
                        if rx > h_rc
                            continue;
                        end
                        for k2 = 1: NZ
                            zp = double(k2)*dz;
                            rz = abs(zi - zp);
                            if rz > h_rc
                                continue;
                            end
                            rsqr = rx*rx + rz*rz;
                            if rsqr < h_rc2
                                XZ_PLANE(k1,k2) = XZ_PLANE(k1,k2) + (1.0 - rsqr*A1) * Norm_const;
                            end
                        end
                    end
                    
                otherwise
                    disp('ERROR 18: check with author');
                    exit();
            end
        end
        rot_count = rot_count + 1;
    end
    XZ_PLANE = XZ_PLANE ./rot_count;
    if INTPOL_TYPE == 1
        XZ_PLANE = XZ_PLANE .* long_cell_density;
    end
    %plotting the density contour (optional)
    figure(1);
    %     subplot(1,2,1);
    [XX,YY] = meshgrid(1:NX,1:NZ);
    ZZ = XZ_PLANE;
    pcolor(ZZ);
    colormap('default');
    shading interp;
    %% %--------------------------------------------
    %--
    % PASS #3 -- Selecting only probable droplet contour pts
    %--
    %--------------------------------------------
    ref_density = max(max(XZ_PLANE));
    threshold_MAX = 0.2 * ref_density; % density threshold MAX 0.2
    threshold_MIN = 0.1 * ref_density; % density threshold MIN 0.02
    xy_COUNT = 1;
    NXNZ = NX*NZ;
    new_x = zeros(NXNZ,1);
    new_y = zeros(NXNZ,1);
    
    for ii = 1:NX
        for jj = 1:NZ
            if XZ_PLANE(ii,jj) > threshold_MIN
                if XZ_PLANE(ii,jj)< threshold_MAX
                    new_x(xy_COUNT) = double(ii);
                    new_y(xy_COUNT) = double(jj);
                    xy_COUNT = xy_COUNT +1;
                end
            end
        end
    end
    new_x(xy_COUNT:NXNZ) = [];
    new_y(xy_COUNT:NXNZ) = [];
    xy_COUNT = length(new_x);
    
%     [xy_COUNT,new_x,MODY,new_y] = F02_find_MD(new_x,0,new_y,MD_cutoffMAX,MD_cutoffMIN,drop_height/dz); % Mahalanobis distance
    %% %--------------------------------------------
    %--
    % PASS #4 -- Removing monolayer
    %--
    %--------------------------------------------
    [xc,yc,R] = F04_Landau_new(new_x,new_y); % initial fit to circle for cleaning
    if yc > R % drop not touching the platinum
        %do nothing
    else % drop touching platinum
        ii = 1;
        while  ii <= xy_COUNT
            if new_y(ii) <= cutoff
                new_x(ii) = [];
                new_y(ii) = [];
                xy_COUNT = xy_COUNT -1;
            else
                ii = ii + 1;
            end
        end
    end
    xy_COUNT = length(new_x);
    MD_cutoffMIN = 0;
    %[xy_COUNT,new_x,new_y] = F02_find_MD(new_x,new_y,MD_cutoffMAX,MD_cutoffMIN); % Mahalanobis distance
    xx = new_x;
    yy = new_y;
    %% %--------------------------------------------
    %--
    % PASS #? -- write filtered x,z to a file
    %--
    %--------------------------------------------
    if create_file == 1
        F06_write_file('final_data.xyz',xx,0,yy,xy_COUNT,Low_L,Low_H,dx,dz,'Ar');
    end
    %% %--------------------------------------------
    %--
    % PASS #5.c -- fit with circle
    %--
    %--------------------------------------------
    [xc,yc,R] = F04_Landau_new(xx,yy); % fitting to a circle [FINAL]
    if create_file == 1 % creating file of atoms representing fitted circle
        theta=0:pi/700:2*pi-pi/180;
        [sdsx, temp_cnt] = size(theta);
        xcircle = R*cos(theta')+xc;
        ycircle = R*sin(theta')+yc;
        F06_write_file('fitted_circle.xyz',xcircle,0,ycircle,temp_cnt,Low_L,Low_H,dx,dz,'Xe');
    end
    SStot = var(yy)*xy_COUNT; % var gives variance
    SSres = 0;
    fi = [];
    for ii = 1:xy_COUNT
        fi = [fi ; yc + real(sqrt(R*R - (xx(ii) - xc)^2))];
        SSres = SSres + (yy(ii) - fi(ii))^2;
    end
    R2value = (1 - SSres/SStot);
    % solving for contact angle
    y0 = cutoff;
    arg1 = R*R -(y0- yc)^2;
    if arg1 >= 0
        if yc < cutoff % Hydrophilic
            x0 = xc - sqrt(arg1);
            gradient = (xc - x0)/(y0 - yc);
            angle = atand(gradient);
            c = y0 - gradient*x0;
            lstart = x0 - 15;
            lend = x0 + 15;
            lxx = lstart:2:lend;    % tangent line
            lyy = gradient.*lxx + c;% tangent line
            lx = [lstart ; lend];   % cutoff line
            ly = [y0 ; y0];         % cutoff line
        else % Hydrophobic
            x01 = xc + sqrt(arg1);
            x02 = xc - sqrt(arg1);
            if x01 < xc
                x0 = x01;
            else
                x0 = x02;
            end
            gradient = (xc - x0)/(y0 - yc);
            angle = atand(gradient);
            if angle < 0
                angle = 180+angle;
            end
            c = y0 - gradient*x0;
            lstart = x0 - 15;
            lend = x0 + 15;
            lxx = lstart:2:lend;    % tangent line
            lyy = gradient.*lxx + c;% tangent line
            lx = [lstart ; lend];   % cutoff line
            ly = [y0 ; y0];         % cutoff line
        end
        
    else
        temp_str = sprintf('No solution for dataset number: %d',i);
        disp (temp_str);
        angle = -1;
    end
    RESULT = [RESULT; i angle R2value ];
    %end
    fprintf('Analysis progress %3.2f percentage\n',100*(i-START_DATA)/TOT_DATA);
end

disp('See RESULT for final results');

theta=0:pi/180:2*pi;
xcircle = R*cos(theta')+xc;
ycircle = R*sin(theta')+yc;
figure(2);
% subplot(1,2,2);
plot(xx,yy,'.',lxx,lyy,lx,ly,xcircle,ycircle,'LineWidth',2);
%axis([0,100,-85,100]);
axis equal;
disp('done');
toc;
