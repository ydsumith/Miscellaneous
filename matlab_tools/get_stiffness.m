clc; clear;

error_flag = 0;
dx_cutoff = 35;
buffer_len = 20000;

list_of_files = dir(fullfile('*.res'));
N_FILES = length(list_of_files);

fid_steps = fopen('break_steps.txt','r');
if fid_steps  == -1,     error_flag = 1; end

if N_FILES == 0 || error_flag == 1
    disp('No files found..exiting');
else
    STEPS = zeros(N_FILES,1);
    break_strain = zeros(N_FILES,3);
    i = 0;
    while ~feof(fid_steps)
        i = i + 1;
        line = fgets(fid_steps); %# read line by line
        STEPS(i) = sscanf(line,'%d');
    end
    if i ~= N_FILES, disp('no: of files and break_steps.txt doesnt match'); end
    fclose(fid_steps);
    
    DATA = zeros(buffer_len, N_FILES * 2);
    for nfiler = 1:N_FILES
        current_file = list_of_files(nfiler).name;
        break_step = STEPS(nfiler);
        
        counter = 0;
        disp('Please wait..');
        fid = fopen(current_file,'r');
        i = 0;
         while ~feof(fid)
            line = fgets(fid); %# read line by line
            line_dat = sscanf(line,'%d %f %f %f');
            if length(line_dat) == 4
                i = i+1;
                dx = line_dat(2);
                
                if dx <= dx_cutoff
                    DATA(i,2*nfiler-1) = dx; % dx
                    DATA(i,2*nfiler)  = line_dat(3); % force
                    N = i;
                end
                
                if line_dat(1) == break_step
                    break_strain(nfiler,1) = line_dat(4); % break strain
                    break_strain(nfiler,2) = line_dat(3); % force
                    
                    %--estimate stiffness--
                    X = DATA(1:N,2*nfiler-1);
                    Y = DATA(1:N,2*nfiler);
                    p = polyfit(X,Y,1);
                    slope = p(1);
                     break_strain(nfiler,3) = slope; % force
                    break;
                end
                
                if mod(i,500) == 0
                    clc;
                    xsp = sprintf('current step is : %d, dx = %f, file = %s', line_dat(1), dx, current_file);
                    disp(xsp);
                end
            end
        end
        fclose(fid);
        
        disp(['Last file processed: ',current_file]);
    end
    
    %------------------------------------
    %---write output file----------------
    disp('writing output files');
    outfile = strcat('2.stiffness','.stiff');
    fid = fopen(outfile,'w');
    for j = 1:N_FILES
        A = list_of_files(j).name;
        B = strsplit(A,'.out');
        C = strjoin(cellstr(B(1)), '');
        fprintf(fid, 'dx\t%s\t',C);
    end
    fprintf(fid,'\n');
    for i = 1:buffer_len
        all_zeros = 1; % 1 yes, 0 no
        for j=1:N_FILES
            dx = DATA(i,2*j-1);
            force = DATA(i,2*j);
            if force == 0
                if i > 10
                    continue;
                end
            end
            all_zeros = 0; % no
            fprintf(fid, '%f\t%f\t', dx, force);
        end
        if all_zeros == 0
            fprintf(fid,'\n');
        else
            break; % stop
        end
    end
    fclose(fid);
    %-------writing break strain file -------
    outfile = strcat('1.break_strain','.strain');
    fid = fopen(outfile,'w');
    fprintf(fid, 'break_strain_(percent)\tforce_(pN)\tslope_(pN/nm)\n');
    for i = 1:N_FILES
        fprintf(fid, '%f\t%f\t%f\n',break_strain(i,1),break_strain(i,2),break_strain(i,3));
    end
    fclose(fid);
end
disp('Program completed successfully');