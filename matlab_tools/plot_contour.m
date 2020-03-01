clc;clear;

% reads a x,y,Z data, convert to [X], [Y], [Z] and plot contour.

error_flag = 0; % 1 error, 0 no
[inputfile,path] = uigetfile('plot_data.txt','Select a file');

if isequal(inputfile,0)
    disp('User selected Cancel');    error_flag = 1;
else
    disp(['User selected ', fullfile(path,inputfile)]);
    
    answer = questdlg('Would you like to continue?', ...
        'Dessert Menu', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            disp('Proceeding..please wait');
        case 'No'
            error_flag = 1;
    end
end

input_full_file = fullfile(path,inputfile);
fid = fopen(input_full_file,'r');
if fid  == -1,     error_flag = 1; end

if error_flag == 0
    %---estimate total lines
    N = 0;
    while ~feof(fid)
        fgets(fid); %# read line
        N = N +1;
    end
    fclose(fid);
    %----read data
    galaxy = zeros(N,5);
    fid = fopen(input_full_file,'r');
    i = 0;
    while ~feof(fid)
        i = i + 1;
        line = fgets(fid); %# read line by line
        galaxy(i,:) = sscanf(line,'%f %f %f %f %f'); % change this depending on the columns
    end
    fclose(fid);
    %--prepare---
    xmin=10000;
    xmax = 0;
    ymin = 10000;
    ymax = 0;
    for i = 1:N
        if  galaxy(i,1) < xmin
            xmin =galaxy(i,1);
        end
        if galaxy(i,1) > xmax
            xmax = galaxy(i,1);
        end
        if  galaxy(i,2) < ymin
            ymin =galaxy(i,2);
        end
        if galaxy(i,2) > ymax
            ymax = galaxy(i,2);
        end
    end
    
    NX = 0;
    NY = 0;
    x = zeros(N,1)-1;
    y = x;
    
    for i = 1:N
        for j =1: N
            if x(j) == galaxy(i,1) %duplicate skip
                break;
            elseif x(j) == -1
                x(j) = galaxy(i,1);
                NX = NX + 1;
                break;
            end
        end
        
        for j =1: N
            if y(j) == galaxy(i,2) %duplicate skip
                break;
            elseif y(j) == -1
                y(j) = galaxy(i,2);
                NY = NY + 1;
                break;
            end
        end
    end
    
    x = x(1:NX);
    y = y(1:NY);
    
    [X,Y] = meshgrid(x,y);
    strain = zeros(NY,NX);
    Youngs = zeros(NY,NX); % MPa
    Stress = zeros(NY,NX); % MPa
    
    k = 0;
    for i = 1:NX
        for j = 1:NY
            k = k + 1;
            strain(j,i) = galaxy(k,3);
            Youngs(j,i) = galaxy(k,4);
            Stress(j,i) = 1000*galaxy(k,5);
        end
    end
    
    [xnew, ynew] = meshgrid(linspace(xmin,xmax,100),linspace(ymin,ymax,100));
    strain_new = interp2(X,Y,strain,xnew,ynew, 'spline');
    Youngs_new = interp2(X,Y,Youngs,xnew,ynew, 'spline');
    Stress_new = interp2(X,Y,Stress,xnew,ynew, 'spline');
    
    fontsize = 16;
    
    figure(1);
    [C, h] = contourf(xnew,ynew,strain_new);
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Break Strain %';
    xlabel('K_{scale}','fontsize',fontsize, 'Interpreter','tex');
    ylabel('\epsilon^*','fontsize',fontsize, 'Interpreter','tex');
    set(gca,'FontSize',fontsize)
    
    figure(2);
    [C, h] = contourf(xnew,ynew,Youngs_new);
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Elastic Modulus (MPa)';
    xlabel('K_{scale}','fontsize',fontsize, 'Interpreter','tex');
    ylabel('\epsilon^*','fontsize',fontsize, 'Interpreter','tex');
    set(gca,'FontSize',fontsize)
    
    figure(3);
    [C, h] = contourf(xnew,ynew,Stress_new);
    colormap(jet);
    c = colorbar;
    c.Label.String = 'Break Stress (kPa)';
    xlabel('K_{scale}','fontsize',fontsize, 'Interpreter','tex');
    ylabel('\epsilon^*','fontsize',fontsize, 'Interpreter','tex');
    set(gca,'FontSize',fontsize)
end

disp('completed');