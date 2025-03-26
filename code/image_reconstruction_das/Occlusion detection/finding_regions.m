clc; close all; clear;

% Load image
image1 = open("image_3_angles.mat");
grayImg = image1.grayImg;

% env = h.data(:,:,f);
% im = 20*log10(grayImg./max(grayImg(:)));
% im = (im-min(im(:)))/(max(im(:))-min(im(:)));
grayImg(grayImg < 0) = 0;  % Set negative values to 0
grayImg(grayImg > 1) = 1;  % Set values greater than 1 to 1
% grayImg = (grayImg-min(grayImg(:)))/(max(grayImg(:))-min(grayImg(:)));
% grayImg = im;
figure();
imshow(grayImg);


grayImg = 1-grayImg;
figure();
imshow(grayImg);
title('Beamformed data')
% saveas(gcf, 'Results/Get_ROI/sim/original.png');


% spatialSigma = 0.1;  % Controls the spatial spread
% intensitySigma = 0.1;  % Controls intensity similarity
% filteredImg = imbilatfilt(grayImg, intensitySigma, spatialSigma);
% 
% figure();
% imshow(filteredImg);

filteredImg2 = medfilt2(grayImg, [16 16]);  % Apply a 3x3 median filter


figure();
imshow(filteredImg2);
title('Median-filtered Beamformed data')
% saveas(gcf, 'Results/Get_ROI/sim/median.png');

filteredImg5 = exp(-8*(1-filteredImg2));

figure();
imshow(filteredImg5);
title('Exponential decay-1')
% saveas(gcf, 'Results/Get_ROI/sim/exp1.png');

filteredImg6 = medfilt2(filteredImg5, [20 20]);
figure();
imshow(filteredImg6);
title('Exponential decay-1 median filter')
saveas(gcf, 'Results/Get_ROI/sim/exp1_med.png');

% filteredImg5 = exp(1*(1-filteredImg2));
% 
% figure();
% imshow(filteredImg5);
% title('Exponential decay-1')
% saveas(gcf, 'Results/Get_ROI/sim/sim1.png');
% 
% filteredImg6 = medfilt2(filteredImg5, [15 15]);
% figure();
% imshow(filteredImg6);
% title('simonential decay-1 median filter')
% saveas(gcf, 'Results/Get_ROI/sim/sim1_med.png');

filteredImg5 = exp(-2*(1-filteredImg6));
figure();
imshow(filteredImg5);
title('Exponential decay-2')
saveas(gcf, 'Results/Get_ROI/sim/exp2.png');

filteredImg6 = medfilt2(filteredImg5, [20 20]);
figure();
imshow(filteredImg6);
title('Exponential decay-2 median filter')
saveas(gcf, 'Results/Get_ROI/sim/exp2_med.png');

binaryMask = filteredImg6 > 0.48;

binaryMask = imclose(binaryMask, strel('disk', 3)); % Close small gaps
binaryMask = imopen(binaryMask, strel('disk', 2));  % Remove small noise

figure();
imshow(binaryMask);
title('Initial Binary Mask')
saveas(gcf, 'Results/Get_ROI/sim/bin_init.png');
% filteredImg4 = medfilt2(binaryMask, [15 15]);  % Apply a 3x3 median filter
% figure();
% imshow(filteredImg4);

filteredImg7 = medfilt2(binaryMask, [30 30]); 
figure(); 
imshow(filteredImg7);

filteredImg7 = imclose(filteredImg7, strel('disk', 3)); % Close small gaps
filteredImg7 = imopen(filteredImg7, strel('disk', 2));  % Remove small noise
figure(); 
imshow(filteredImg7);
title('Median filter Binary Mask')
saveas(gcf, 'Results/Get_ROI/sim/bin_fin.png');

filteredImg7 = medfilt2(filteredImg7, [30 30]); 
figure(); 
imshow(filteredImg7);
title('Median filter Binary Mask')
saveas(gcf, 'Results/Get_ROI/sim/bin_fin.png');
% filteredImg4 = medfilt2(filteredImg4, [13 13]);  % Apply a 3x3 median filter
% figure();
% imshow(filteredImg4);
% 
% laplacianFilter = fspecial('laplacian', 0.2); % Laplacian filter with default alpha
% filteredImg3 = imfilter(filteredImg2, laplacianFilter, 'replicate');
% 
% figure();
% imshow(filteredImg3);
% Apply adaptive thresholding
% Apply median filtering to reduce noise
% I_filt = medfilt2(grayImg, [5, 5]);

% Apply adaptive thresholding
bw = imbinarize(grayImg, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.4);

% Perform morphological opening to remove small noise
bw = imopen(bw, strel('disk', 5));

% Fill holes inside objects
bw = imfill(bw, 'holes');

% Identify connected components
stats = regionprops(bw, 'Area', 'Eccentricity', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength');

% Create an empty mask
filtered_bw = false(size(bw));

% Select only elliptical objects
for i = 1:numel(stats)
    if stats(i).Eccentricity > 0.7 && stats(i).MajorAxisLength / stats(i).MinorAxisLength > 1.2 ...
            && stats(i).Area > 500
        bbox = round(stats(i).BoundingBox);
        filtered_bw(bbox(2):(bbox(2)+bbox(4)), bbox(1):(bbox(1)+bbox(3))) = 1;
    end
end

% Display final result
imshow(filtered_bw);


% Reshape the image into a vector for clustering
% [m, n] = size(grayImg); % Get image dimensions
% pixelValues = double(grayImg(:)); % Convert grayscale image to a 1D vector
% 
% % Perform K-means clustering on ALL pixels
% numClusters = 2; 
% [idx, C] = kmeans(pixelValues, numClusters, 'Replicates', 3);
% 
% % Ensure idx has correct dimensions
% if length(idx) ~= m * n
%     error('Unexpected size of idx. Expected %d but got %d.', m*n, length(idx));
% end
% 
% % Reshape clustered labels back into the image shape
% segmentedImg = reshape(idx, m, n); % Corrected reshaping
% 
% % Convert to binary mask (assuming the darker regions are the ROIs)
% binaryMask = segmentedImg == mode(segmentedImg(:)); % Select dominant low-intensity cluster
% 
% % Morphological processing to refine detections
% binaryMask = imclose(binaryMask, strel('disk', 5)); % Close small gaps
% binaryMask = imopen(binaryMask, strel('disk', 3));  % Remove small noise
% 
% % Find connected components
% CC = bwconncomp(binaryMask);
% stats = regionprops(CC, 'BoundingBox', 'Area');
% 
% % Filter out small or irrelevant regions
% minArea = 200; % Adjust based on expected ROI size
% filteredROIs = [];
% for i = 1:length(stats)
%     if stats(i).Area > minArea
%         filteredROIs = [filteredROIs; stats(i).BoundingBox];
%     end
% end
% 
% % Display detected regions
% figure, imshow(grayImg), hold on;
% for i = 1:size(filteredROIs,1)
%     rectangle('Position', filteredROIs(i,:), 'EdgeColor', 'r', 'LineWidth', 2);
% end
% title('Detected Regions of Interest');
% hold off;
% 
% % Compute vertical gradient manually
% Compute vertical gradient
vertgrad = zeros(size(grayImg));
vertgrad(2:end, :) = diff(double(grayImg), 1, 1); % First derivative along rows

% Apply thresholding
threshold = 0.03;
binaryMask = vertgrad > threshold;

% Morphological processing to refine detections
binaryMask = imclose(binaryMask, strel('disk', 3)); % Close small gaps
binaryMask = imopen(binaryMask, strel('disk', 2));  % Remove small noise

% Initialize mobile mask
mobile = zeros(size(grayImg));

% Sliding window approach to check vertical continuity
for i = 1:size(grayImg,1)-10
    for j = 1:size(grayImg,2)  % Fixed incorrect loop
        k = 0;
        for p = i:i+10
            if binaryMask(p, j) == 0 % Fixed condition
                k = k + 1;
            end
        end
        if k >= 10
            mobile(i, j) = 1;
        end
    end
end

% Display result
imshow(mobile);
title('Mobile Mask');

% Display gradient and histogram
% figure(); imshow(vertgrad, []);
% figure(); histogram(vertgrad);

% Apply thresholding
threshold = 0.03;
binaryMask = vertgrad > threshold;

% Morphological processing to refine detections
binaryMask = imclose(binaryMask, strel('disk', 3)); % Close small gaps
binaryMask = imopen(binaryMask, strel('disk', 2));  % Remove small noise

mask = zeros(size(binaryMask));
linmask = zeros(size(mask, 2), 1);
linmask = linmask+10;

for i=5:size(mask,1)-5
    
end



% Remove small noise adaptively
stats = regionprops(binaryMask, 'Area');
areaValues = [stats.Area]; % Extract numerical areas
areaThreshold = prctile(areaValues, 30); % Dynamic threshold
cleanedImage = bwareaopen(binaryMask, round(areaThreshold)); 

% Label connected components
CC = bwconncomp(cleanedImage);
stats = regionprops(CC, 'Area', 'Eccentricity', 'Solidity');

% Create a new mask keeping only circular objects
finalMask = false(size(cleanedImage));
for i = 1:length(stats)
    if stats(i).Eccentricity < 0.7 && stats(i).Solidity > 0.9 % More robust circularity check
        finalMask(CC.PixelIdxList{i}) = 1;
    end
end

% Fill holes in detected circles
finalMask = imfill(finalMask, 'holes');

% Invert to set circles as 0 and rest as 1
finalMask = ~finalMask;

% Display results
figure;
subplot(1,2,1); imshow(binaryMask); title('Binary Mask After Thresholding');
subplot(1,2,2); imshow(finalMask); title('Final Mask with Circles as 0');

figure(); imshow(binaryMask);

% Remove small noise using area opening
cleanedImage = bwareaopen(binaryMask, 500); % Adjust 500 based on noise size

% Label connected components
CC = bwconncomp(cleanedImage);
stats = regionprops(CC, 'Area', 'Eccentricity');

% Create a new mask keeping only circular objects
finalMask = false(size(cleanedImage));
for i = 1:length(stats)
    if stats(i).Eccentricity < 0.7 % Keep circular shapes
        finalMask(CC.PixelIdxList{i}) = 1;
    end
end

% Invert to set circles as 0 and rest as 1
finalMask = ~finalMask;

% Display results
figure;
subplot(1,2,1); imshow(binaryImage); title('Original Binary Image');
subplot(1,2,2); imshow(finalMask); title('Cleaned Mask with Circles as 0');

% Morphological operations to refine detection
binaryMask = imclose(binaryMask, strel('disk', 3)); % Close small gaps
binaryMask = imopen(binaryMask, strel('disk', 2));  % Remove noise

% Find connected components
CC = bwconncomp(binaryMask);
stats = regionprops(CC, 'BoundingBox', 'Area');

% Filter small regions
minArea = 100; % Adjust based on expected ROI size
filteredROIs = [];
for i = 1:length(stats)
    if stats(i).Area > minArea
        filteredROIs = [filteredROIs; stats(i).BoundingBox];
    end
end

% Display results
figure, imshow(grayImg), hold on;
for i = 1:size(filteredROIs,1)
    rectangle('Position', filteredROIs(i,:), 'EdgeColor', 'r', 'LineWidth', 2);
end
title('Detected Regions of Interest');
hold off;


%-- Setting axis limits (mm)
x_lim = [min(scan.x_matrix(:)) max(scan.x_matrix(:))]*1e3; 
z_lim = [min(scan.z_matrix(:)) max(scan.z_matrix(:))]*1e3; 
x = scan.x_matrix;
z = scan.z_matrix;

padding = 1;
maskOcclusion = cell(1,9);
maskInterference = cell(1,9);

for k=1:length(pht.occlusionDiameter)

    r = pht.occlusionDiameter(k) / 2;
    xc = pht.occlusionCenterX(k);
    zc = pht.occlusionCenterZ(k);
    maskOcclusion{k} = ( ((x-xc).^2 + (z-zc).^2) <= r^2);
    
    x_dash = pht.occlusionCenterX;
    x_dash(k) = []; 
    z_dash = pht.occlusionCenterZ;
    z_dash(k) = [];

    maskInterference{k} = false(size(x));

    % Apply mask for all 8 circles
    for i = 1:length(x_dash)
        maskInterference{k} = maskInterference{k} | ...
            ((x - x_dash(i)).^2 + (z - z_dash(i)).^2 <= r^2);
    end

    figure();
    imshow(maskOcclusion{k});
    title(sprintf('Occlusion %2d', k));
    % filename = fullfile('Results/Count_Weight_Anglewise3/', num2str(k), '/','mask_occlusion.jpg');
    % saveas(gcf, ['Results/Count_Weight_Anglewise3/', num2str(k), '/','mask_occlusion.jpg'])

    figure();
    imshow(maskInterference{k});
    title('Occlusion %2d', k);
    % saveas(gcf, ['Results/Count_Weight_Anglewise3/', num2str(k), '/','mask_interference.jpg'])
end

% counts = cell(1,9);
% weights = cell(1,9);
angles = linspace(-16, 16, 75);
count_angles_all = cell(1,2);
weight_angles_all = cell(1,2);

count_angles = zeros(size(angles));
weight_angles = zeros(size(angles));
[height, width] = size(maskOcclusion{i});
count_all = cell(1,75);
weight_all = cell(1,75);

for i=1:2
    % Define angles for projection
    % Starting points (bottom edge of the image)
    x_start = 1:1:width; % Project lines at intervals along the bottom
    
    % figure; imshow(maskOcclusion{i}); hold on;
    % figure; imshow(maskInterference{i}); hold on;
    
    k=1;
    % p=1;
    % Loop through each starting point and angle
    for theta=angles
        count = zeros(size(x_start));
        weight = ones(size(x_start));
        % count = 0;
        rad = deg2rad(theta);
        % i=1;
        for j = 1:length(x_start)
            % Convert angle to radians
            % weight = 1;
            % Project a long enough line upwards
            for inc=1:height
                x_end = x_start(j) + inc * tan(rad);
                y_end = inc;
                
                dx = x_end - floor(x_end);
                % Bilinear interpolation for maskInterference
                if(x_end>=1 && x_end<=width)
                    Q11 = maskInterference{i}(y_end, floor(x_end));
                    Q21 = maskInterference{i}(y_end, ceil(x_end));
                    interp_maskInterf = Q11 * (1 - dx) + Q21 * dx;
                

                    % Reduce weight based on interpolated interference
                    weight(j) = weight(j) * (1 - interp_maskInterf * 0.5);
    
                    % Bilinear interpolation for maskOcclusion
                    Q11 = maskOcclusion{i}(y_end, floor(x_end));
                    Q21 = maskOcclusion{i}(y_end, ceil(x_end));
                    interp_maskOccl = Q11 * (1 - dx) + Q21 * dx;
    
                    % Accumulate occlusion count using interpolated value
                    count(j) = count(j) + weight(j) * interp_maskOccl;
                end
                
                % if(maskInterference{i}(y_end,x_end))
                %     weight = weight*0.5;
                % end
                % 
                % if(maskOcclusion{i}(y_end,x_end))
                %     count = count + weight; 
                % end
            end

            % i = i+1;
        end
        % if mod(k,10) == 0
        %     figure();
        %     plot(x_start, count);
        %     title(sprintf('Count vs x for \\theta = %.2f°, occlusion = %2d ', theta, i))
        %     xlabel('horzontal distance (transducer distance)')
        %     ylabel('Count on transducer')
        %     saveas(gcf, ['Results/Count_Weight_Anglewise3/', num2str(i), '/', num2str(theta), '_count.jpg'])
        % 
        %     figure();
        %     plot(x_start, weight);
        %     title(sprintf('Weight vs x for \\theta = %.2f°, occlusion = %2d', theta, i))
        %     xlabel('horzontal distance (transducer distance)')
        %     ylabel('Weight on transducer')
        %     saveas(gcf, ['Results/Count_Weight_Anglewise3/', num2str(i), '/', num2str(theta), '_weight.jpg'])
        % end

        count_angles(k) = sum(count);
        count_all{k} = count;

        weight_angles(k) = sum(weight);
        weight_all{k} = weight;
        k = k+1;
    end
    count_angles_all{i} = count_angles;
    weight_angles_all{i} = weight_angles;
    counts{i} = count_all;
    weights{i} = weight_all;
end

temp_count = zeros(size(angles));
temp_weight = zeros(size(angles));
for h=1:2
    temp_count = temp_count+count_angles_all{h};
    temp_weight = temp_weight+weight_angles_all{h};
    
    figure();
    plot(linspace(-16,16,75), count_angles_all{h});
    title('Total Count vs Angle for Occlusion %2d', h);
    % saveas(gcf, ['Results/Count_Weight_Anglewise3/', num2str(h), '/','count_angle_occlusion.jpg'])

    figure();
    plot(linspace(-16,16,75), weight_angles_all{h});
    title('Total Weight vs Angle for Occlusion %2d', h);
    % saveas(gcf, ['Results/Count_Weight_Anglewise3/', num2str(h), '/','cweight_angle_occlusion.jpg'])
end

figure();
plot(linspace(-16,16,75), temp_count);
title('Total Count vs Angle for All Occlusion');
% saveas(gcf, 'Results/Count_Weight_Anglewise3/Total_count_angle_occlusion.jpg')

figure();
plot(linspace(-16,16,75), temp_weight);
title('Total Weight vs Angle for All Occlusion');
% saveas(gcf, 'Results/Count_Weight_Anglewise3/Total_weight_angle_occlusion.jpg')
