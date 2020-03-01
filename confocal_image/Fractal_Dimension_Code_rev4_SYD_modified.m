clear;
clc;

error_flag = 0; % 1 error, 0 no
[filename,path] = uigetfile('*.tif','Select a file');

if isequal(filename,0)
    disp('User selected Cancel');    error_flag = 1;
else
    disp(['User selected ', fullfile(path,filename)]);
    inputfile = fullfile(path,filename);
    
    I = imread(inputfile);
    
    BW1 = im2bw(I,graythresh(I));
    BW2 = imadjust(mat2gray(BW1));
    BW3 = fibermetric(BW2, 30, 'StructureSensitivity', 15);
    BW = BW3 > 0.00002;
    
    %BW tubular filter above misses the most intense areas (HOT)
    HOT1=imadjust(I,[0.97 1.0],[]);
    HOT2=im2bw(HOT1,graythresh(HOT1));
    HOT3=imadjust(mat2gray(HOT2));
    HOT4=imboxfilt(HOT3,25);
    HOT = HOT4 > 0.6;
    
    %add the most intense (HOT) back to the tubular(BW) results
    BW = BW + HOT;
        
    % detect the edge of image 'p' using the Canny algorithm
    % this gives edge as 'e2'
    im = im2bw(BW, graythresh(BW));
    e  = edge(double(im));
    fi = imfill(im, 'holes');
    op = imerode(fi,strel('disk',4));
    e2 = edge(double(op));
    figure(2)
    imshow(e2)
    % once we have e2, set up a grid of blocks across the image
    % and scan each bloch too see if the edge occupies any of the blocks.
    % If a block is occupied then flag it and record it in boxCount --
    % store both size of blocks (numBlocks) and no of occupied boxes (boxCount)
    % in table()
    Nx = size(e2,1);
    Ny = size(e2,2);
    for numBlocks = 1:25
        
        sizeBlocks_x = floor(Nx./numBlocks);
        sizeBlocks_y = floor(Ny./numBlocks);
        
        flag = zeros(numBlocks,numBlocks);
        for i = 1:numBlocks
            for j = 1:numBlocks
                xStart = (i-1)*sizeBlocks_x + 1;
                xEnd   = i*sizeBlocks_x;
                
                yStart = (j-1)*sizeBlocks_y + 1;
                yEnd   = j*sizeBlocks_y;
                
                block = e2(xStart:xEnd, yStart:yEnd);
                
                flag(i,j) = any(block(:)); %mark this if ANY part of block is true
            end
        end
        boxCount = nnz(flag);
        table(numBlocks,1)    = numBlocks;
        table(numBlocks,2)    = boxCount;
    end
    % from the above table of discrete points, take a line of best fit and plot
    % the raw data (ro) and line of best fit (r-)
    x      = table(:,1);    % x is numBlocks
    y      = table(:,2);    % y is boxCount
    
    p       = polyfit(x,y,1);
    BestFit = polyval(p,x);
    
    figure(3)
    hold on
    grid on
    plot(x,y,       'ko','LineWidth',1)
    plot(x,BestFit, 'k-','LineWidth',2)
    xlabel('Number of blocks, N','FontSize',12)
    ylabel('Box Count, N(s)','FontSize',12)
    
    
    % calculate Hausdorff Dimension
    x2 = log(x);
    y2 = log(y);
    
    p2       = polyfit(x2,y2,1);
    BestFit2 = polyval(p2,x2);
    
    figure(4)
    hold on
    grid on
    plot(x2,y2,       'bo','LineWidth',1)
    plot(x2,BestFit2, 'b-','LineWidth',2)
    xlabel('Number of blocks, log N','FontSize',12)
    ylabel('Box Count, log N(s)'    ,'FontSize',12)
    
    HausdorffDimension = p2(:,1);
end

disp('Program completed');
