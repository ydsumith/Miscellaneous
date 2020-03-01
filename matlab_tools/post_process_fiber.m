clc; clear;

error_flag = 0;
dx_cutoff = 2000;
step_cutoff = 254900;

buffer_len = step_cutoff + 10;

list_of_files = dir(fullfile('*.out'));
N_FILES = length(list_of_files);

if N_FILES == 0
    disp('No files found..exiting');
else
    for nfiler = 1:N_FILES
        current_file = list_of_files(nfiler).name;
        
        if error_flag == 0 % no errors and proceed
            step = zeros(buffer_len,1);
            x_left = zeros(buffer_len,1);
            x_right = zeros(buffer_len,1);
            Fx_left = zeros(buffer_len,1);
            Fx_right = zeros(buffer_len,1);
            Fx_fiber = zeros(buffer_len,1);
            
            counter = 0;
            disp('Please wait..');
            fid = fopen(current_file,'r');
            i = 0;
            while ~feof(fid)
                line = fgets(fid); %# read line by line
                data = strsplit(line);
                if strcmp(data(1),'Step')==1
                    counter = counter +1;
                end
                if counter >= 3
                    data = sscanf(line,'%d %f %f %f %f %f %f %f %f %f %f %f');
                    if length(data) == 12
                        i = i+1;
                        step(i) = data(1);
                        x_left(i) = data(8);
                        x_right(i) = data(9);
                        Fx_left(i) = data(10);
                        Fx_right(i) = data(11);
                        Fx_fiber(i) =  data(12);
                        if length(x_right)>1
                            %                     if x_right(i) - x_right(1) > dx_cutoff
                            %                        break;
                            %                     end
                            if step(i)>step_cutoff
                                break;
                            end
                        end
                        if mod(step(i),500) == 0
                            clc;
                            xsp = sprintf('current step is : %d, dx = %f', step(i), x_right(i) - x_right(1));
                            disp(xsp);
                        end
                    end
                end
            end
            fclose(fid);
            
            %---write output file---
            outfile = strcat(current_file,'.res');
            
            N = length(Fx_fiber);
            dx = x_right - x_right(1);
            L0 = x_right(1) - x_left(1);
            strain = 100 * dx / L0;
            
            fid = fopen(outfile,'w');
            fprintf(fid,'Step\tdx(nm)\tforce(pN)\tstrain(%%)\n');
            for i = 1:N
                %         if dx(i) > dx_cutoff, break; end
                if dx(i) < -10, break; end
                fprintf(fid, '%d\t%f\t%f\t%.3f\n', step(i), dx(i), Fx_right(i), strain(i) );
            end
            fclose(fid);
            
            disp('Program completed successfully');
            disp(['Last file processed: ',current_file]);
        else
            disp('Nothing executed');
        end
    end
end
% [inputfile,path] = uigetfile('*.out');
% if isequal(inputfile,0)
%     disp('User selected Cancel');
%     error_flag = 1;
% else
%     disp(['User selected ', fullfile(path,inputfile)]);
% end



