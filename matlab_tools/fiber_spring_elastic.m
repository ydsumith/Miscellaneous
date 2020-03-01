% Program by Sumith Yesudasan
% Fibrin Fiber dynamics simulation
clc; clear;
error_flag = 0;
tic;

EXT_FORCE = 1;
FRICTION = 1;
runfile = 'kscale_15';
%-----------------------------------------
%--simulation parameters--
timestep = 30; % ps
thermo_intvl = 1; % write thermo output every this many steps
dump_intvl = 100; % write dump files every this many steps
dump_chunk = 250000; % write a new dump file every 'dump_chunk' steps
pull_start = 1;

pull_rate = 0.36324; % LJ units of velocity (1 m/s = 0.36324)
K0_scale = 1 / 15; % multiply this number to K0
epsilon_LJ = 1.5;

%%
%----constants-------
T = 310; % K
kB = 1.38064852e-23; % (m2 kg s-2 K-1)
perm0 = 8.85418782e-12; % (m^-3 kg^-1 s^4 A^2)
rho_MD = 33.369; % no: of waters per nm^3
N_avogadro = 6.0221409e+23;
K0 = 0.0139626; % N/m or J/m2 of fiber--dont change this

%%
%----variables-------
NA = 20; % no: of beads in a protofibril element

mass_of_bead = 340000; % amu or g/mol
kangle = 1000; % in LJ units
theta0 = 180;

rcut = 10; % nm
L_fiber = 2000; % nm
skin = 50; % (nm) increase box size
radius = 51; % (nm) radius of fiber [Cell Biochem Biophys (2007) 49:165–181]
drad = 7; % (nm)

N = 10000; % no: of protofibrils to start with

B_attract = 1e-9; % (N) inter fiber friction
B_attract = B_attract / NA; % convert to udl



%%
%---initial calculations---

%--length of protofibril is 450 nm
bondlength = 450/(NA-1); % nm

gap = bondlength; % (nm) gap between protofibrils
sigma_nm = rcut; % (nm)
cell_size = rcut;
K0 = K0 * (NA-1) * K0_scale; % individual springs
dS = 1.1*drad; % (nm)

%%
%----- conversion to LJ units
sigma = sigma_nm * 1e-9; % nm -> m
epsilon = kB*T;

bondlength = bondlength / sigma_nm;
L_fiber = L_fiber / sigma_nm;
gap = gap / sigma_nm;
radius  = radius / sigma_nm;
drad = drad / sigma_nm;
dS = dS / sigma_nm;
rcut = rcut / sigma_nm;
rcut2 = rcut*rcut;
cell_size = cell_size / sigma_nm;
skin = skin / sigma_nm;

K0 = K0 * (sigma* sigma / epsilon); % LJ units
B_attract = B_attract * sigma / epsilon; % LJ units (force)
%pull_force = pull_force * sigma / epsilon; % LJ units (force)
T_ref = T * kB / epsilon;

Ncells = ceil(3*L_fiber / cell_size); % create 3 times larger box for extension
cell_size = 3*L_fiber / Ncells;
Ncells = ceil(3*L_fiber / cell_size);

atomlist = zeros(Ncells,100); % max 100 atoms in a cell
neighlist = zeros(Ncells,5); % max 5 neighbors for a cell

REFX = zeros(NA,1);
REFY = zeros(NA,1);
REFZ = zeros(NA,1);
MASS = zeros(NA,1);
length_proto = gap + (NA-1) *bondlength;
L = L_fiber + length_proto;

if L < 10*2*radius
    XSP = sprintf('Warning! L (%.2f) < 10D (%.2f)',L,10*2*radius);
    disp(XSP);
    disp('Try increasing L');
end

galaxy = zeros(N*NA, 5);
BOX = zeros(3,2);
BOX(1,1) = 0; % xlo
BOX(1,2) = 4*L; % xhi
BOX(2,1) = -radius - skin; % ylo
BOX(2,2) = radius + skin; % yhi
BOX(3,1) = -radius - skin; % zlo
BOX(3,2) = radius + skin; % zhi

for i = 1:NA
    REFX(i) = (i-1) * bondlength;
    MASS(i) = mass_of_bead;
end

mref = mass_of_bead / (1000*N_avogadro); % kg/bead
tau = sqrt(mref*sigma*sigma/epsilon); %time conversion constant
timestep_star = timestep * 1e-12 / tau; % LJ units
half_dt = timestep_star/2;
half_dt_sqr = half_dt*half_dt;

T_star = T * kB / epsilon;
MASS = MASS ./ (1000 * 6.0221409e+23); % convert to SI units
MASS = MASS ./ mref; % convert to reduced units

REF_MOL = zeros(NA,5);
REF_MOL(:,2) = REFX;
REF_MOL(:,3) = REFY;
REF_MOL(:,4) = REFZ;
%%
N = 1;
dtheta = (2*pi/10.0);
for irad = 0:drad:radius
    circum = 2*pi*irad;
    Ntemp = ceil(circum / dS);
    dStemp = circum / Ntemp;
    if Ntemp == 0
        dStemp = dS;
    end
    i = 1;
    for dcircum = 0:dStemp:circum-dStemp
        theta = dcircum/irad;
        if irad == 0
            theta = 0;
        end
        if mod(i,2)==0
            pitch = length_proto/2;
        else
            pitch = 0;
        end
        for jlen = 0:length_proto:L
            
            xpos = jlen + pitch + skin;
            ypos = irad*cos(theta);
            zpos = irad*sin(theta);
            galaxy((N-1)*NA+1:N*NA, 1) = N;
            galaxy((N-1)*NA+1:N*NA, 2) = REF_MOL(:, 2) + xpos;
            galaxy((N-1)*NA+1:N*NA, 3) = REF_MOL(:, 3) + ypos;
            galaxy((N-1)*NA+1:N*NA, 4) = REF_MOL(:, 4) + zpos;
            N = N + 1;
        end
        i = i + 1;
    end
end

N = N-1;
Natoms = NA*N;
sizer = length(galaxy(:,1));
galaxy(Natoms+1:sizer,:)=[];

%---create neighlist
neighlist(:,1) = -1;
for i = 1:Ncells-1
    ipos = (i - 0.5) * cell_size;
    for j = i + 1: Ncells
        jpos = (j - 0.5) * cell_size;
        abs_rx = abs(ipos - jpos);
        if abs_rx <= rcut
            for k = 1:5
                if neighlist(i,k) == -1
                    neighlist(i,k) = j;
                    neighlist(i,k+1) = -1;
                    break;
                else
                    continue;
                end
            end
        end
    end
end
%%
%---create bonds---
nbonds = (NA-1)*N;
bonds = zeros(nbonds,2);
ibond = 1;
for i = 1: Natoms-1
    if galaxy(i,1) == galaxy(i+1,1) % same molecule
        bonds(ibond,1) = i;
        bonds(ibond,2) = i+1;
        ibond = ibond + 1;
    end
end

%----create angles--
nangles = (NA-2)*N;
angles = zeros(nangles,3);
iangle = 1;
for i = 1: Natoms-2
    if galaxy(i,1) == galaxy(i+2,1) % same molecule
        angles(iangle,1) = i;
        angles(iangle,2) = i+1;
        angles(iangle,3) = i+2;
        iangle = iangle + 1;
    end
end

%---create end conditions---
galaxy(:,5) = 1; % setting type of bead
up_cut = max(galaxy(:,2)) - length_proto;
low_cut = min(galaxy(:,2)) + length_proto;
for i = 1: Natoms
    if galaxy(i,2) > up_cut
        galaxy(i,5) = 2;
    elseif galaxy(i,2) < low_cut
        galaxy(i,5) = 3;
    end
end

%---create force, velocity, displacement memories--
force = zeros(Natoms, 3);
acceleration = zeros(Natoms, 3);
velocity = zeros(Natoms, 3);
displacement = galaxy;
velocity(:,1) = rand(Natoms,1)./10;
velocity(:,1) = velocity(:,1) - mean(velocity(:,1));

%---write preliminary data---
% dump_file = sprintf('K0_%.2f_Battrac_%.2f_fiber.mol', K0, B_attract);
% log_file = sprintf('K0_%.2f_Battrac_%.2f_output.txt', K0, B_attract);
%
% fid = fopen(dump_file,'w');
% fclose(fid);
% fid = fopen(log_file,'w');
% XSP = sprintf('dx\tright_force\n');
% fprintf(fid,XSP);
% fclose(fid);
% step_val = 1;
% write_mol( displacement , step_val, Natoms, BOX, gap,dump_file);
%
% right = [];
% left = [];
% right_x = max(displacement(:,2)) - bondlength/2;
% left_x = min(displacement(:,2))  + bondlength/2;
% for i=1:Natoms
%     if displacement(i,2) > right_x
%         right = [right i];
%     elseif displacement(i,2) < left_x
%         left = [left i];
%     end
% end
% nright = length(right);
% nleft = length(left);
% ref_displacement_right = displacement(right(1),2);

%%
%---main simulation module
% disp('Please wait..simulation running...');
% for step_val = 1:total_steps
%
%     %----find displacement
%     displacement(:,2) = displacement(:,2) + timestep .* velocity(:,1) + half_dt_sqr .* acceleration(:,1);
%     displacement(:,3) = displacement(:,3) + timestep .* velocity(:,2) + half_dt_sqr .* acceleration(:,2);
%     displacement(:,4) = displacement(:,4) + timestep .* velocity(:,3) + half_dt_sqr .* acceleration(:,3);
%
%     %----find first half velocity
%     velocity(:,1)  = velocity(:,1)  + half_dt .* acceleration(:,1);
%     velocity(:,2)  = velocity(:,2)  + half_dt .* acceleration(:,2);
%     velocity(:,3)  = velocity(:,3)  + half_dt .* acceleration(:,3);
%
%     %----clear force
%     force = force * 0;
%
%     %----cell list updates
%     atomlist = atomlist .* 0 - 1;
%     for i = 1:Natoms
%         cellnum = ceil(displacement(i,2)/cell_size);
%         for k = 1:100
%             if atomlist(cellnum,k) == -1
%                 atomlist(cellnum,k) = i;
%                 atomlist(cellnum,k+1) = -1;
%                 break;
%             else
%                 continue;
%             end
%         end
%     end
%
%     %----find intra bond force----
%     for i=1:nbonds
%         iatom = bonds(i,1);
%         jatom = bonds(i,2);
%         rx = displacement(iatom,2) - displacement(jatom,2);
%         ry = displacement(iatom,3) - displacement(jatom,3);
%         rz = displacement(iatom,4) - displacement(jatom,4);
%         r = sqrt(rx*rx + ry*ry + rz*rz);
%         dr = r - bondlength;
%         rk = K0 * dr;
%         fbond = -2.0 * rk / r;
%         force (iatom,1) = force (iatom,1) + fbond * rx;
%         force (iatom,2) = force (iatom,2) + fbond * ry;
%         force (iatom,3) = force (iatom,3) + fbond * rz;
%
%         force (jatom,1) = force (jatom,1) - fbond * rx;
%         force (jatom,2) = force (jatom,2) - fbond * ry;
%         force (jatom,3) = force (jatom,3) - fbond * rz;
%     end
%
%     %----find inter protofibril friction--
%     if FRICTION == 1
%         %----using cell list
%         for icell=1: Ncells
%             for i = 1:99
%                 if atomlist(icell,i) == -1
%                     break;
%                 else
%                     iatom = atomlist(icell,i);
%                     imol = displacement(iatom,1);
%                     for j = i+1:100
%                         if atomlist(icell,j) == -1
%                             break;
%                         else
%                             jatom = atomlist(icell,j);
%                             jmol = displacement(jatom,1);
%                             if imol == jmol, continue; end % skip for same molecule
%                             rx = displacement(iatom,2) - displacement(jatom,2);
%                             ry = displacement(iatom,3) - displacement(jatom,3);
%                             rz = displacement(iatom,4) - displacement(jatom,4);
%                             rsqr = rx*rx + ry*ry + rz*rz;
%                             if rsqr > rcut2, continue; end
%
%                             r = sqrt(rsqr);
%                             wd = 1 - r/rcut;
%                             fpair = -B_attract * wd;
%
%                             force (iatom,1) = force (iatom,1) + fpair * rx;
%                             force (iatom,2) = force (iatom,2) + fpair * ry;
%                             force (iatom,3) = force (iatom,3) + fpair * rz;
%
%                             force (jatom,1) = force (jatom,1) - fpair * rx;
%                             force (jatom,2) = force (jatom,2) - fpair * ry;
%                             force (jatom,3) = force (jatom,3) - fpair * rz;
%                         end
%                     end
%                 end
%             end
%         end
%
%         for icell=1: Ncells
%             for k = 1:5
%                 if neighlist(icell,k) == -1
%                     break;
%                 else
%                     jcell = neighlist(icell,k);
%                     for i = 1:100
%                         if atomlist(icell,i) == -1
%                             break;
%                         else
%                             iatom =  atomlist(icell,i);
%                             imol = displacement(iatom,1);
%                             for j = 1:100
%                                 if atomlist(jcell,j) == -1
%                                     break;
%                                 else
%                                     jatom = atomlist(jcell,j);
%                                     jmol = displacement(jatom,1);
%                                     if imol == jmol, continue; end % skip for same molecule
%                                     rx = displacement(iatom,2) - displacement(jatom,2);
%                                     ry = displacement(iatom,3) - displacement(jatom,3);
%                                     rz = displacement(iatom,4) - displacement(jatom,4);
%                                     rsqr = rx*rx + ry*ry + rz*rz;
%                                     if rsqr > rcut2, continue; end
%
%                                     r = sqrt(rsqr);
%                                     dr = r - rcut;
%                                     wd = 1 - dr/rcut;
%                                     fpair = -B_attract * wd;
%
%                                     force (iatom,1) = force (iatom,1) + fpair * rx;
%                                     force (iatom,2) = force (iatom,2) + fpair * ry;
%                                     force (iatom,3) = force (iatom,3) + fpair * rz;
%
%                                     force (jatom,1) = force (jatom,1) - fpair * rx;
%                                     force (jatom,2) = force (jatom,2) - fpair * ry;
%                                     force (jatom,3) = force (jatom,3) - fpair * rz;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%
%     end
%     %---external force application
%     if EXT_FORCE == 1
%         if step_val == pull_start
%             ref_displacement_left = min(displacement(:,2));
%             ref_displacement_right = max(displacement(:,2));
%         end
%         if step_val > pull_start
%             for i = 1:nleft
%                 iatom = left(i);
%                 force(iatom,1) = 0;
%                 displacement(iatom,2) = ref_displacement_left;
%             end
%             for i = 1:nright
%                 iatom = right(i);
%                 displacement(iatom,2)= displacement(iatom,2) + pull_rate;
%             end
%         end
%     end
%
%     %----find acceleration---
%     acceleration = force ./ MASS(1);
%
%     %----find second half velocity
%     velocity(:,1)  = velocity(:,1)  + half_dt .* acceleration(:,1);
%     velocity(:,2)  = velocity(:,2)  + half_dt .* acceleration(:,2);
%     velocity(:,3)  = velocity(:,3)  + half_dt .* acceleration(:,3);
%
%     %-----temperature adjustment ----
%     vxx = sum(velocity(:,1).^2);
%     vyy = sum(velocity(:,2).^2);
%     vzz = sum(velocity(:,3).^2);
%     KE = vxx + vyy + vzz;
%
%     KE = KE * 0.5 * MASS(1);
%     T_BULK = 2* KE / (3 * Natoms);
%     scale =  sqrt(T_ref / T_BULK); %Velocity scaling
%     velocity = velocity .* scale;
%
%     %----write output---
%     if mod(step_val, dump_intvl)==0
%         write_mol( displacement , step_val, Natoms, BOX, gap,dump_file);
%         message = sprintf('%d percentage completed..', round(100*step_val/total_steps));
%         disp(message);
%         message = sprintf('Temperature = %f', T_BULK);
%         disp(message);
%
%         fid = fopen(log_file,'a');
%         left_force = 0;
%         for i = 1:nleft
%             iatom = left(i);
%             left_force = left_force + force(iatom,1);
%         end
%
%         right_force = 0;
%         dx = displacement(right(1),2) - ref_displacement_right;
%         for i = 1:nright
%             iatom = right(i);
%             right_force = right_force + force(iatom,1);
%         end
%         XSP = sprintf('%f\t%f\n', dx,  -right_force);
%         fprintf(fid,XSP);
%         fclose(fid);
%
%     end
% end
%%
[status,msg,msgID] = mkdir(runfile);
if strcmp(msgID,'MATLAB:MKDIR:DirectoryExists')==1, disp('directory already exists'); end
outfile_temp = fullfile(runfile,'lammps.data');
fid = fopen(outfile_temp,'w');
fprintf(fid, '# timestep = %f (fs), sigma = %f (nm)\n', tau*timestep_star*1e+15, sigma_nm);
fprintf(fid, '\n%d atoms \n', Natoms);
fprintf(fid, '%d bonds \n',nbonds);
fprintf(fid, '%d angles \n\n',nangles);

fprintf(fid, '%d atom types \n', 3);
fprintf(fid, '%d bond types \n', 1);
fprintf(fid, '%d angle types \n\n', 1);

XSP = sprintf('%f %f xlo xhi \n', BOX(1,1), BOX(1,2));
fprintf(fid, XSP);
XSP = sprintf('%f %f ylo yhi \n', BOX(2,1), BOX(2,2));
fprintf(fid, XSP);
XSP = sprintf('%f %f zlo zhi \n \n', BOX(3,1), BOX(3,2));
fprintf(fid, XSP);
fprintf(fid, 'Masses \n \n');

fprintf(fid, '%d %.4f \n', 1, MASS(1));
fprintf(fid, '%d %.4f \n', 2, MASS(1));
fprintf(fid, '%d %.4f \n\n', 3, MASS(1));

fprintf(fid, 'Atoms #atomid, molid, type, charge, x, y, z\n \n');
j = 1;
for i=1:Natoms
    %atom-ID molecule-ID atom-type q x y z
    XSP = sprintf('%d %d %d %.2f %.3f %.3f %.3f\n', i, ...
        displacement(i,1), displacement(i,5), 0, displacement(i,2), displacement(i,3), displacement(i,4));
    fprintf(fid,XSP);
    if mod(i,NA)== 0
        j = j + 1;
    end
end

if nbonds ~= 0
    fprintf(fid, '\nBonds \n \n');
    j = 1;
    for i = 1:nbonds
        XSP = sprintf('%d %d %d %d\n',i, 1, bonds(i,1), bonds(i,2));
        fprintf(fid,XSP);
    end
end

if nangles ~= 0
    fprintf(fid, '\n Angles \n \n');
    j = 1;
    ang_type = 1;
    for i = 1:nangles
        XSP = sprintf('%d %d %d %d %d\n',i, ang_type , angles(i,1), angles(i,2), angles(i,3));
        fprintf(fid,XSP);
    end
end
fclose(fid);
disp('Coordinates written to lammps.data');

%%
%---writing pair style table
Nr = 4999; % number of table points (r values)
r_start = 0.00001;
r_end = rcut;
d0 = 0.5 * rcut; % equm distance at beads interact
U_A = 0.01; % strength of the potential
U_B = 12; % spread of the function (min 10)
U_C = -2*U_A*U_B*(rcut - d0)* exp(-U_B*(rcut-d0)^2);
U_D = (1+2*rcut*U_B*(rcut-d0))*U_A*exp(-U_B*(rcut-d0)^2);
dr1 = (r_end-r_start)/Nr;
r = r_start:dr1:r_end;
potential = -U_A*exp(-U_B*((r-d0).^2)) + U_C*rcut + U_D;
force = -2*U_A*U_B*(r-d0).*exp(-U_B*(r-d0).^2) - U_C;
% plot(r,potential, r, force);

% fid = fopen('gauss.table','w');
% fprintf(fid, '# Fiber elastic - Sumith Yesudasan - potential for inter protofibril friction \n');
% fprintf(fid, 'GAUSS\n');
% fprintf(fid, 'N %d R %f %f\n\n', Nr+1, r_start, r_end);
% for i = 1: Nr+1
%     fprintf(fid, '%d %f %f %f\n', i, r(i), potential(i),force(i));
% end
% fprintf(fid, '\n');
% fclose(fid);

%%
%---writing the LAMMPS input script
disp('Now writing lammps input file to conf.in');
outfile_temp = fullfile(runfile,'conf.in');
fid = fopen(outfile_temp,'w');

fprintf(fid, '# Fiber elastic - Sumith Yesudasan \n');

%fprintf(fid, '\npackage omp 4\n');

fprintf(fid, '\nunits lj \n');
fprintf(fid, 'boundary p p p\n');
fprintf(fid, 'atom_style full\n\n');

fprintf(fid, 'pair_style lj/cut %.3f\n\n', rcut);
%fprintf(fid, 'pair_style table linear %d\n\n', Nr+1);
%fprintf(fid, 'pair_style morse %.3f\n\n', rcut);

fprintf(fid, 'pair_modify tail yes table 12\n');

fprintf(fid, 'special_bonds charmm \n');
fprintf(fid, 'bond_style harmonic \n');
fprintf(fid, 'angle_style harmonic \n\n');

fprintf(fid, 'read_data lammps.data \n');

fprintf(fid, 'comm_modify vel yes\n\n');

fprintf(fid, 'log res_log.log \n\n');

sigma_LJ = drad / 2^(1/6);
fprintf(fid, '# epsilon, sigma, cutoff1, cutoff2 \n');
fprintf(fid, 'pair_coeff * * %.3f %.3f\n\n',epsilon_LJ, sigma_LJ);

%fprintf(fid, '#pair_coeff * 3 morse.table ENTRY1 cutoff\n');
%fprintf(fid, 'pair_coeff * * gauss.table GAUSS %f\n\n',rcut);

% De = 100;
% alpha = 0.02;
% fprintf(fid, '#pair_coeff * * D0 alpha r0\n');
% fprintf(fid, 'pair_coeff * * %.3f %.3f %.3f\n\n',De, alpha, d0);

fprintf(fid, 'bond_coeff %d %.2f %.2f\n\n', 1, K0, bondlength);

fprintf(fid, 'angle_coeff %d %.2f %.2f\n\n', 1, kangle, theta0);

fprintf(fid, 'timestep %.8f \n', timestep_star);

neigh_skin = 3*(2*bondlength) - rcut; % required to calculate angle potential
fprintf(fid, 'neighbor %.2f bin\n',neigh_skin);

fprintf(fid, 'neigh_modify delay 0 every 1\n\n');

fprintf(fid, 'group fiber type 1\n');
fprintf(fid, 'group left type 3\n');
fprintf(fid, 'group ends type 2 3\n');
fprintf(fid, 'group right type 2\n\n');

%binsize = 5 / sigma_nm;
%fprintf(fid, 'compute density all chunk/atom bin/1d x lower %.3f\n\n',binsize);

fprintf(fid, 'compute LFX left reduce sum fx\n');
fprintf(fid, 'compute LFY left reduce sum fy\n');
fprintf(fid, 'compute LFZ left reduce sum fz\n\n');

fprintf(fid, 'compute RFX right reduce sum fx\n');
fprintf(fid, 'compute RFY right reduce sum fy\n');
fprintf(fid, 'compute RFZ right reduce sum fz\n\n');

fprintf(fid, 'compute midFX fiber reduce sum fx\n');
fprintf(fid, 'compute midFY fiber reduce sum fy\n');
fprintf(fid, 'compute midFZ fiber reduce sum fz\n\n');

fprintf(fid, 'compute extL left reduce min x\n');
fprintf(fid, 'compute extR right reduce max x\n\n');

% 2 = Nose Hoover
fprintf(fid, 'fix 1 all nvt temp %.5f %.5f %.3f\n\n',T_ref, T_ref, 100 * timestep_star );

fprintf(fid, 'variable sigma equal %10.6e\n', sigma);
fprintf(fid, 'variable epsilon equal %10.6e\n', epsilon);
fprintf(fid, 'variable kB equal %10.6e\n', kB);
fprintf(fid, 'variable tau equal %10.6e\n', tau);
fprintf(fid, 'variable time equal dt*step*v_tau\n');
fprintf(fid, 'variable P equal press*v_epsilon*1e-5/(v_sigma*v_sigma*v_sigma)\n');
fprintf(fid, 'variable T equal temp*v_epsilon/v_kB \n\n');

fprintf(fid, 'variable Ffactor equal 1e+12*v_epsilon/v_sigma\n'); % converting to pN
fprintf(fid, 'variable xfactor equal "1e+9 * v_sigma"\n'); % converting to nm
fprintf(fid, 'variable natoms equal "count(all)"\n');
fprintf(fid, 'variable extR equal "v_xfactor * c_extR"\n');
fprintf(fid, 'variable extL equal "v_xfactor * c_extL"\n');
fprintf(fid, 'variable LF equal "v_Ffactor * sqrt(c_LFX*c_LFX + c_LFY*c_LFY + c_LFZ*c_LFZ)/v_natoms"\n');
fprintf(fid, 'variable RF equal "v_Ffactor * sqrt(c_RFX*c_RFX + c_RFY*c_RFY + c_RFZ*c_RFZ)/v_natoms"\n');
fprintf(fid, 'variable midF equal "v_Ffactor * sqrt(c_midFX*c_midFX + c_midFY*c_midFY + c_midFZ*c_midFZ)/v_natoms"\n\n');

fprintf(fid, 'thermo %d\n', thermo_intvl);
fprintf(fid, '#extR (nm), force (pN)\n');
fprintf(fid, 'thermo_style custom step temp v_T pe etotal press v_P v_extL v_extR v_LF v_RF v_midF\n\n');


fprintf(fid, 'dump 1 all custom %d dump1.mol id mol type  x y z \n', dump_intvl);


fprintf(fid, 'velocity all create 1.0 312345\n');
fprintf(fid, 'run 0   # temperature may not be 300K\n');
fprintf(fid, 'velocity all scale 1.0   # now it should be\n\n');

%fprintf(fid, 'min_style cg\n');
%fprintf(fid, 'minimize 1.0e-6 1.0e-8 500 1000\n\n');

fprintf(fid, 'run 5000 \n\n');

fprintf(fid, 'unfix 1\n');
fprintf(fid, 'fix 1 fiber nvt temp %.5f %.5f %.3f\n',T_ref, T_ref, 100 * timestep_star );
fprintf(fid, 'fix 2 right move linear %.3f 0 0\n\n',pull_rate);

counter = 1;
expectd_strain  = 125; % in percentage
pull_per_step_nm = pull_rate*timestep_star*sigma_nm; % this much nm pull every step
pull_per_step = pull_per_step_nm / sigma_nm; % in LJ units
pull_nm_per_ns = pull_per_step_nm/(timestep/1000); % nm/ns
total_steps = round((expectd_strain/100) * (L_fiber/pull_per_step));

for i = 1:total_steps
    if mod(i, dump_chunk)==0
        counter = counter +1;
        fprintf(fid, 'undump 1\n');
        fprintf(fid, 'dump 1 fiber custom %d dump%d.mol id mol type  x y z \n', dump_intvl, counter);
        fprintf(fid, 'run %d\n\n',dump_chunk);
        
    end
    if (total_steps-i) < dump_chunk
        fprintf(fid, 'undump 1\n');
        fprintf(fid, 'dump 1 fiber custom %d dump%d.mol id mol type  x y z \n', dump_intvl, counter+1);
        fprintf(fid, 'run %d\n\n',total_steps-i);
        break;
    end
end

fprintf(fid, 'write_data final.data\n');
fclose(fid);
%------------------------------------------------------------------
%%
%--writing the submission script
write_submit_script(runfile);
%------------------------------------------------------------------
disp('Program completed ');
toc;
