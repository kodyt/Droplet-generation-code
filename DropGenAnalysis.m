close all


%This script was written for the use of droplet & circle identification
%feel free to use it as needed, if you do have questions though, MATLAB
%help might be the best option to start with, search "imfindcircles" as a
%starting point, or let me know if there are issues and we can def improve
%upon it
%SCLP

%%
SRadius = 8;   % The minimum droplet/microsphere radius (in pixels) the script will seek out (and recognize)
LRadius = 14;   % The maximum droplet/microsphere radius (in pixels) the script will seek out (and recognize)
%Generally you want to know relative pixel size, if your radius range is
%too big accuracy for the script goes down. However, if you do not know the
%approximate pixel size, try a run from 10 to 1000 and identify approximate
%size and then update accordingly
% 

% TODO 
% Figure out how to get rid of manual radii 

% 0.91
SensitivityFactor = 0.94;   % The sensitivity of the circle-hunting command.  Lower values will 
                            % lead the script to be more selective, while higher values will 
                            % lead to it being less selective.




%% TODO %% 
% Change the senstivity factors and then make the average size of the
% particles
%% TODO %%

%%%%% START KODY's CHANGES %%%%%%
startingFolder = 'C:\Program Files\MATLAB';
if ~exist(startingFolder, 'dir')
  % If that folder doesn't exist, just start in the current folder.
  startingFolder = pwd;
end

% Choose the files you want to 
[file,path] = uigetfile('*.*', 'Select One or More Files', ...
    'MultiSelect', 'on');
if isequal(file,0)
   disp('User selected Cancel');
   return;
% else
%    disp(['User selected ', fullfile(path,file)]);
end


% TODO: Preallocate a vertical array
% Create a vertical array of all zeros
% Keep track how far you've indexed with a counter
% Append the new data
% Remove any extra zeros by taking size - counter and resize
% Analyze
% radiiBright = zeros(10000000, 1)

imgs = fullfile(path,file);

% Checking if it is one or more images
if ischar(imgs(1))
    imgs = {imgs};
    file = {file};
end

% TODO
% Create an output file and then clear it if it already exists
folderPath = fullfile(pwd, 'Output-Data');
if exist(folderPath, 'dir')
    rmdir(folderPath, 's');
end
mkdir(folderPath);
cd(folderPath); % This might be unnecessary?


radiiBright = 0;
% radiiBright = zeros(1000000,1);
% Preallocate and then with the counter from the temp, at the end resize
% the radiibright -= the total number of circles found

% Reads all of the files in the new directory
for k = 1:length(imgs)
    % fprintf(1, 'Now reading %s\n', fullFileName);
    RawImage = imread(string(imgs(k)));
    % imshow(imageArray);  % Display image.
    % drawnow; % Force display to update immediately.
    
    B = RawImage; % This is all your information for analysis
    % LofI = length(B); % Determining the length of the matrix, to know dimensions (normally all images are squares so one dimension is fine)
    %%
    red = B(:,:,1); % Red channel data
    % green = B(:,:,2); % Green channel data
    % blue = B(:,:,3); % Blue channel data
    a = zeros(size(B, 1), size(B, 2)); % Makes an all zero matrix to blot out channels that are not of interest, i.e they become 0

    %%
    %Producing matrices that are only composed of a single channel's worth of data
    
    just_red = cat(3, red, a, a);     % Creates a new set of data only composed of red channel data
    % just_green = cat(3, a, green, a); % Creates a new set of data only composed of green channel data
    % just_blue = cat(3, a, a, blue);  % Creates a new set of data only composed of blue channel data
    % back_to_original_img = cat(3, red, green, blue); % Back to your original image such that you can evaluate
    
    %Preparing image analysis and output
    
    L = bwlabel(red); %labels connected components in a 2D array/
    s = regionprops(L, 'PixelIdxList', 'PixelList', 'Area', 'Centroid', 'FilledArea'); %measures properties of image (i.e. property definitions) should not need to change
    
    [centersBright, tempRadiiBright] = imfindcircles(just_red,[SRadius LRadius],'ObjectPolarity','bright','Sensitivity',SensitivityFactor);%Algorithm that identifies the circles
    
    % TODO 
    % With preallocation, instead of adding to the end of the array have it
    % just index into the section
    % Could hold the current index of the latest value and then replace it with the tempRadiibright and then add the latest value
    radiiBright = [radiiBright; tempRadiiBright];
    
    figure()
    %% Individual photos 
    subplot(2,1,2), imshow(RawImage)
    viscircles(centersBright, tempRadiiBright,'EdgeColor','r') %places circles around your objects, should help visualize the specific beads or droplets identified
    
    %% Individual Photo Data
    DiaBrightAdj = tempRadiiBright*(1100/410)*2; 
    subplot(2,1,1), histogram(DiaBrightAdj)
    title('Droplet size distribution','FontSize', 14) 
    xlabel('Diameter (um)', 'FontSize', 14) 
    ylabel('Number of beads', 'FontSize', 14) 
    countDia = size(DiaBrightAdj,1);
    DiaMean = mean(DiaBrightAdj);    
    DiaSTD = std(DiaBrightAdj);

    A = [countDia; DiaMean; DiaSTD];
    DataInfo = {'Samples Counted: %.3f';...
               'Mean Diameter: %.3f';...
               'St.Dev: %.3f'};
    str = compose(DataInfo, A);
    dim = [0.55 0.6 0.3 0.3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');

    sgtitle(string(file(k))) % Title of the individual photo

    %% TODO 
    % Output the figure to the output file and DO NOT print
end
    
% The value inside of the parentheses below is the um / pixels
DiaBrightAdj = radiiBright*(1100/410)*2; % Gives you diameter of the beads/droplets... NEED TO MAKE SURE CORRECT MULTIPLIER IS USED  

figure,histogram(DiaBrightAdj) % Gives you a histogram of the data
title('Droplet size distribution','FontSize', 14) % title of your plot
xlabel('Diameter (um)', 'FontSize', 14) % x-axis label of your plot
ylabel('Number of beads', 'FontSize', 14) % y-axis label of your plot 
countDia = size(DiaBrightAdj,1); %Data that may be useful to have present
DiaMean = mean(DiaBrightAdj);    %Data that may be useful to have present
DiaSTD = std(DiaBrightAdj);      %Data that may be useful to have present

%Information that will be placed on the histogram, useful if you want to
%visualize output, if you want to eliminate comment out
A = [countDia; DiaMean; DiaSTD];
DataInfo = {'Samples Counted: %.3f';...
           'Mean Diameter: %.3f';...
           'St.Dev: %.3f'};
str = compose(DataInfo, A);
dim = [0.55 0.6 0.3 0.3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');


%%% END KODY'S CHANGES %%%

%%
%This portion below is only necessary to visualize, the different versions
%of the image you may have created, names should be descriptive enough to use
%I would suggest to comment most or all of them out if you are not using them

% figure, imshow(B), title('Original image')
% figure, imshow(just_red), title('Red channel')
% figure, imshow(just_green), title('Green channel')
% figure, imshow(just_blue), title('Blue channel')
% figure, imshow(back_to_original_img), title('Back to original image')


%% MUlTIPLIER SECTION
%If using the Celestron microscope images make sure to put in the
%appropriate pixel to micron conversion, same goes for nikon, etc... 
%multipliers are in micron/pixel

%Celestron 
%For 4X use 1.6447 as multiplier
%For 10X use 0.6964 as multiplier

%Nikon
%For 10x use 1.24 as multiplier
%For 20x use 0.62 as multiplier
%

%EVOS LP Lab
%2X objective = 1302 pixels per 5618um
%4X objective = 1216 pixels per 2620um
%10X objective = 1724 pixels per 1500um
%10X multiplier would be 1500um/1724pix
%20X objective = 1793 pixels per 800um