%% --- clean previous data ---
clearvars -except i data k N peri; clc; 

%% --- read images into a data-structre  ---
% ask the pictures path from the user, save it into "selectedPath" variable
selectedPath = uigetdir('.');
% take the pictures (.jpg suffix) from selectedPath, save it into "dirOutput" struct
cd(selectedPath);
dirOutput = dir('*.jpg');
% extract data for the variables above, from the fields in "dirOutput" struct
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);
% read the data of the first picture into the variable "Pic1"
Pic1 = imread(fileNames{1});
% preallocate array - create zeros 4-D uint8 variable, with the size of Pic1 * numFrames
imageSeries = zeros([size(Pic1) numFrames],class(Pic1));
% fill the first numerical pixels values of Pic1
imageSeries(:,:,:,1) = Pic1;
% fill the rest numerical pixels values into the imageSeries struct
for p = 2:numFrames
    imageSeries(:,:,:,p) = imread(fileNames{p});
end

%% --- cropping ---
% crop the images menually
img = imageSeries(:,:,:,1); % keep the first image
% initialize new data-structure for the images after cropping
croppedImageSeries = [];
% initialize a "flag" for the drawing loop. if it's 1=>draw again, if it's 0=>crop
draw = 1;
% drawing loop
while draw == 1
clear croppedImageSeries; % clean previous versions of cropped images
figure(1) % open visual box
imshow(img) % show the first image in the visual box
h = msgbox('Draw a rectengle to crop the image'); % open another visual box with intructions
pause(2) %in seconds
% let the user create a rectangle on the screen, save the rectangle positions (coordinates) in r1
r1 = drawrectangle('Label','','Color',[1 0 0]);
%crop all the images by r1 positions, keep the cropped images in "croppedImageSeries" data-structure
for k = 1:numFrames
    croppedImageSeries(:,:,:,k) = imcrop(imageSeries(:,:,:,k),r1.Position);
end
% display the cropped images to the user
figure(1)
title('Cropped images')
delete(h);
montage(croppedImageSeries) % Display multiple image frames as rectangular montage
% ask for user feedback
selection = questdlg('Would you like to continue?', ...
	'Confirm cropped image series', ...
	'Confirm','Crop again','Crop again');
% if the approve, stop the loop. else, do it all over again
if strcmp(selection,'Confirm')
    draw = 0;
end
end
close(figure(1))

%% --- preprocessing ---
% Short explanation for the software developer:
% Write in the Matlab command line : size(croppedImageSeries). you should
% get 4 numbers. for example: 154   277     3    24. it means:
% height   width   RGB*   number of images in the structe.
% RGB = red, green, blue. for each image, each pixel have those 3 colors, expressed by numerical value.
% in the croppedImageSeries case, the value in index 1 = red values, index 2 = green values, index 3 = blue values. 

for k = 1:numFrames
% split channels - green channel (index 2 of the third dimension - as i explained above)
% greenChannelSeries keeps only the green values of the pictures
greenChannelSeries(:,:,k) = croppedImageSeries(:,:,2,k);

% ### J = imnoise(greenChannelSeries(:,:,k),'salt & pepper',0.02);
% Gaussian filter - Smooth Image. 2 = std
guass(:,:,k) = imgaussfilt(greenChannelSeries(:,:,k),2);
% ### figure(1) % open visual box
% ### imshow(guass(:,:,k)) % show the first image in the visual box
% ### medfilter for anomalies (as bubbles)
% ### guass(:,:,k) = medfilt2(J);

% ### figure(2) % open visual box
% ### imshow(guass(:,:,k)) % show the first image in the visual box

% find threshold
OtsuLevel(:,k) = graythresh(guass(:,:,k));
% apply threshold - "bw" conatins pixel values of "1" to the dark areas, "0" to the green areas
bw(:,:,k) = imcomplement(imbinarize(guass(:,:,k), OtsuLevel(:,k)));

% Enframe the gut (the unchecked area), 8 pixels away
radius = 8;
decomposition = 0;
% strel command create Morphological structuring element
se = strel('disk', radius, decomposition);
% bw data-structure now conatins the "1" pixels minus 8 pixels
bw(:,:,k) = imerode(bw(:,:,k), se);
end

%% --- mask the inverted image with the thresholded image ---
writer_obj = VideoWriter('check.avi');
writer_obj.FrameRate = 5;
open(writer_obj);

figure
for k = 1:numFrames
%medium mask
Medium_Bw_uint8(:,:,k) = uint8(bw(:,:,k)); % convert logical to integer
% bring the "1" pixels ( - 8 eight pixels) their real numeric value, keep the 0 pixels as 0
Medium_Masked(:,:,k) = greenChannelSeries(:,:,k) .* Medium_Bw_uint8(:,:,k);
if (k==1)
    lowTemp = min(min(Medium_Masked(:,:,k)));
    highTemp = max(max(Medium_Masked(:,:,k)));
end
tmpLOW = min(min(Medium_Masked(:,:,k)));
if (tmpLOW < lowTemp)
    lowTemp = tmpLOW;
end
tmpHIGH = max(max(Medium_Masked(:,:,k)));
if (tmpHIGH > highTemp)
    highTemp = tmpHIGH;
end
thermalImage = rescale(Medium_Masked(:,:,k), lowTemp, 100*highTemp);
image(thermalImage)
writeVideo(writer_obj, getframe);
axis equal
axis off
title('image ', k);
pause(0.1)
end
close(writer_obj);



%% --- measure grayscale values ---
for k = 1:numFrames
%medium measurements
MediumStats (k)= regionprops(Medium_Bw_uint8(:,:,k), Medium_Masked(:,:,k), 'MeanIntensity'); %mean intensity
MediumMI(k) = MediumStats(k).MeanIntensity % MediumMI values are used in the main script
end

%% --- find perimeter ---
for k = 1:numFrames 
findper(:,:,k) = bwperim(bw(:,:,k)); % perimeter
per(:,:,k) = greenChannelSeries(:,:,k) + 255 * uint8(findper(:,:,k)); % convert,scale, and overlay perimeter on image
end

