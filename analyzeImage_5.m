% Research Script
% SupersonicImagine
%
% Matlab Script that uses extractDicomEmbeddedData function to access
% embedded data within a DICOM file and displays the elastography map
%extractEmbeddedData is a part of the dcmtk library that must be under the
%same folder as this code.
% contact thomas.frappart@supersonicimagine.com for assistance

dicomFileIn = 'images/I0000036';
[SQPTHeader, SQPTRegionHeader, SQPTRegionParameter, SQPTRegionData] = extractDicomEmbeddedData(dicomFileIn);

% Axis
depth   = SQPTRegionHeader.regionTopLeftOrdmm + linspace(0, SQPTRegionHeader.regionHeightmm, SQPTRegionHeader.Samples);
lateral = SQPTRegionHeader.regionTopLeftAbsmm + linspace(0, SQPTRegionHeader.regionWidthmm, SQPTRegionHeader.Lines);

figure,subplot(121)
imagesc(lateral, depth, SQPTRegionData.dataSWE)
colormap(SQPTRegionParameter.LUT.data/255)
axis image; colorbar;
xlabel('Lateral (mm)');ylabel('Depth (mm)');
title('Elasto Box in kPa');

subplot(122)
imagesc(lateral, depth, SQPTRegionData.qualityMap)
colormap(SQPTRegionParameter.LUT.data/255);
axis image; colorbar;
xlabel('Lateral (mm)');ylabel('Depth (mm)');
title('Quality map');
% code till here is from researchSWE.m code


[I,cmap] = dicomread(dicomFileIn); 
figure, 
subplot(121),imshow(I)

% ROI of colorbar. I think the location of colorbar is fixed, so you just
% need to set the values
ROI_colorbar = I(141:360,1386:1403,:);
%subplot(1,3,2),imshow(ROI_colorbar)

% color range of the colorbar in the original image
% 3 columns, each for color channel R, G, and B
% convert the color to MATLAB style, i.e. 0 to 1 instead of 0 to 255
RGB_range = squeeze(double(ROI_colorbar(:,9,:)))./255;

% SWE range, I assume the deep blue is 0 and the red is 10
% we make a linear scale between the minimum and maximum value of SWE
% the number of the division need to be the same as the RGB range division
linear_scale_SWE = linspace(0,10,size(RGB_range,1));
% our RGB range starts with red (the first row) and it should denote SWE
% value of 10, since we make linear scale from 0 to 10, we need to flip the
% array
linear_scale_SWE = fliplr(linear_scale_SWE);

% increase the red color resolution
% find the first color that has full red, i.e the value in R channel is 255
% or 1 if the scale is 0 to 1
idx_R = find(RGB_range(:,1)==1,1,'first');
% interpolate only the red color range by 2
new_scale_R = interp1(1:idx_R,RGB_range(1:idx_R,1),1:0.5:idx_R);
% update the color range, first create a new variable
RGB_range_new = RGB_range;
% remove the previous color range and replace it with the new red color range
RGB_range_new(1:idx_R-1,:) = [];
RGB_range_new=[[new_scale_R',zeros(length(new_scale_R),1),zeros(length(new_scale_R),1)];RGB_range_new];
% do the same thing with the SWE scale (we want to match it with the color so it need to have the same division)
new_scale_SWE_R = interp1(1:idx_R,linear_scale_SWE(1:idx_R),1:0.5:idx_R);
linear_scale_SWE(1:idx_R-1) = [];
linear_scale_SWE=[new_scale_SWE_R,linear_scale_SWE];

% round the value to make the matching process easier
RGB_range_new_rounded = round(RGB_range_new.*10000)./10000;

% find the region of interest using masking method
% loop through the pixel
% the idea is, the ROI is colorful and the background + B mode is black and
% white image. So if the color of the pixel does not match with the color
% in the color bar, it means that it is part of the background image (we
% label it with 0 since it is not our region of interest)
for pixel_row = 1: size(I,1)
    for pixel_col = 1: 1380
        % take the color of the pixel in that position, convert it to
        % MATLAB color ranged from 0 to 1
        pixel_color = double(squeeze(I(pixel_row,pixel_col,:)))./255;
        % round the value to make the matching process easier
        pixel_color = round(pixel_color.*10000)./10000;
        % is the color of the pixel in the color bar?
        idx_color_match = (RGB_range_new_rounded(:,1)==pixel_color(1)) & ...
            (RGB_range_new_rounded(:,2)==pixel_color(2)) & ...
            (RGB_range_new_rounded(:,3)==pixel_color(3));
        % if not in the color bar, it means that the pixel is part of the
        % background or the B mode, label it with 0, otherwise it is the
        % region of interest that we need
        if sum(idx_color_match)==0
            mask_ROI(pixel_row,pixel_col)=0;
        else
            mask_ROI(pixel_row,pixel_col)=1;
        end
    end
end

% find the edge of the ROI
left_edge = find(sum(mask_ROI)>0,1,'first');
right_edge = left_edge-1 + find(sum(mask_ROI(:,left_edge:left_edge+300))>0,1,'last');
top_edge = find(sum(mask_ROI,2)>0,1,'first');
low_edge = top_edge-1 + find(sum(mask_ROI(top_edge:top_edge+300,:),2)>0,1,'last');

% Extract the region of interest
ROI = I(top_edge:low_edge,left_edge:right_edge,:);
subplot(1,2,2),imshow(ROI);

%figure,imshow(ROI);
[r,c] = size(ROI(:,:,1));

%algorithm for image segmentation
I = rgb2lab(ROI);

% seperate color channels
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);

% La*b* color thresholding for red/high, green/medium, blue/low stiffness muscle, it is a more
% robust way of seperating the colors 
% Define thresholds for channel 1 based on histogram settings
redchannel1Min = 0.000;
redchannel1Max = 100.000;

% Define thresholds for channel 2 based on histogram settings
redchannel2Min = 5.000;
redchannel2Max = 80.002;

% Define thresholds for channel 3 based on histogram settings
redchannel3Min = 5.000;
redchannel3Max = 94.410;

% Create mask based on chosen histogram thresholds
red_sliderBW = (R >= redchannel1Min ) & (R <= redchannel1Max) & ...
    (G >= redchannel2Min ) & (G <= redchannel2Max) & ...
    (B >= redchannel3Min ) & (B <= redchannel3Max);
BW_R = red_sliderBW;
% Initialize output masked image based on input image.
redmask = ROI;

% Set background pixels where BW is false to zero.
redmask(repmat(~BW_R,[1 1 3])) = 0;

% Defining the parameters to get green color.
% Define thresholds for channel 1 based on histogram settings
greenchannel1Min = 0.000;
greenchannel1Max = 100.000;

% Define thresholds for channel 2 based on histogram settings
greenchannel2Min = -60.000;
greenchannel2Max = -19.000;

% Define thresholds for channel 3 based on histogram settings
greenchannel3Min = 5.000;
greenchannel3Max = 70.000;

% Create mask based on chosen histogram thresholds
green_sliderBW = (R >= greenchannel1Min ) & (R <= greenchannel1Max) & ...
    (G >= greenchannel2Min ) & (G <= greenchannel2Max) & ...
    (B >= greenchannel3Min ) & (B <= greenchannel3Max);
BW_G = green_sliderBW;

% Initialize output masked image based on input image.
greenmask = ROI;

% Set background pixels where BW is false to zero.
greenmask(repmat(~BW_G,[1 1 3])) = 0;

% Define parameters to get blue color
% Define thresholds for channel 1 based on histogram settings
bluechannel1Min = 0.000;
bluechannel1Max = 100.000;

% Define thresholds for channel 2 based on histogram settings
bluechannel2Min = -60.000; 
bluechannel2Max = 80.000;

% Define thresholds for channel 3 based on histogram settings
bluechannel3Min = -108.000;
bluechannel3Max = -20.000; %-20.000   

% Create mask based on chosen histogram thresholds
blue_sliderBW = (R >= bluechannel1Min ) & (R <= bluechannel1Max) & ...
    (G >= bluechannel2Min ) & (G <= bluechannel2Max) & ...
    (B >= bluechannel3Min ) & (B <= bluechannel3Max);
BW_B = blue_sliderBW;

% Initialize output masked image based on input image.
bluemask = ROI;

% Set background pixels where BW is false to zero.
bluemask(repmat(~BW_B,[1 1 3])) = 0;
figure,
subplot(131), imagesc(lateral, depth,redmask),
axis image;
xlabel('Lateral (mm)');ylabel('Depth (mm)');
subplot(132), imagesc(lateral, depth,greenmask),
axis image; 
xlabel('Lateral (mm)');ylabel('Depth (mm)');
subplot(133), imagesc(lateral, depth,bluemask)
axis image; 
xlabel('Lateral (mm)');ylabel('Depth (mm)');

% Match the pixel color with the SWE value
% Loop for each pixel
for pixel_row = 1: size(ROI,1)
    for pixel_col = 1: size(ROI,2)
        % take the color of the pixel in that position, convert it to
        % MATLAB color ranged from 0 to 1
        pixel_color = double(squeeze(ROI(pixel_row,pixel_col,:)))./255;
        % round the value to make the matching process easier
        pixel_color = round(pixel_color.*10000)./10000;
        % is the color of the pixel in the color bar? If so, give me the
        % index number of the color
        idx_color_match = (RGB_range_new_rounded(:,1)==pixel_color(1)) & ...
            (RGB_range_new_rounded(:,2)==pixel_color(2)) & ...
            (RGB_range_new_rounded(:,3)==pixel_color(3));
        % if it is not in the color bar (maybe the missing SWE data point),
        % label the SWE as nan
        % if it is in the color bar, using the index, find the match SWE
        % value.
        if sum(idx_color_match)==0
            SWE(pixel_row,pixel_col)=nan;
        else
            SWE(pixel_row,pixel_col) = linear_scale_SWE(idx_color_match);
        end
    end
end

% create figure of the ROI and the SWE value on each pixel
% (this require a lot of computer power, it will take a while until you see
% the tiny numbers on the ROI)
figure,
subplot(121),imcontour(SWE,'ShowText','on');
subplot(122),imagesc(lateral, depth, ROI)
axis image;
xlabel('Lateral (mm)');ylabel('Depth (mm)');