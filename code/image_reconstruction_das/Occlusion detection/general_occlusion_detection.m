% Load the mask image
clear;
close all;
clc;
addpath(genpath('../../src'));


%-- Parameters
acquisition_type = 1;       %-- 1 = simulation || 2 = experiments
phantom_type = 2;           %-- 1 = resolution & distorsion || 2 = contrast & speckle quality
data_type = 1;              %-- 1 = IQ || 2 = RF
flag_display = 0;

%-- Parsing parameter choices
switch acquisition_type    
    case 1
        acquisition = 'simulation';
        acqui = 'simu';
        flag_simu = 1;
    case 2
        acquisition = 'experiments';
        acqui = 'expe';
        flag_simu = 0;
    otherwise       %-- Do deal with bad values
        acquisition = 'simulation';
        acqui = 'simu';
        flag_simu = 1;
end
switch phantom_type    
    case 1
        phantom = 'resolution_distorsion';
    case 2
        phantom = 'contrast_speckle';
    otherwise       %-- Do deal with bad values
        phantom = 'resolution';
end
switch data_type    
    case 1
        data = 'iq';
    case 2
        data = 'rf';
    otherwise       %-- Do deal with bad values
        data = 'iq';        
end
%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_pht = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);
pht = us_phantom();
pht.read_file(path_pht);

pw_indices{1} = [38 19 57];
flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);

%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
image1 = das_iq_original(scan,dataset,pw_indices); %path_scan, path_pht, flag_simu, flag_display);
% tools.exec_evaluation_contrast_speckle(path_scan,path_pht,'Results/Paper_sampling_strategy/Reference/all.hdf5',flag_simu,flag_display,'Results/Paper_sampling_strategy/Reference/all.txt');
% contrast = tools.contrast_score(scan,pht, image,15);
disp('Reconstruction Done')

grayImg = image1.data;

grayImg = (grayImg-min(grayImg(:)))/(max(grayImg(:))-min(grayImg(:)));
% figure();
% imshow(grayImg);
grayImg = 1-grayImg;
figure();
imshow(grayImg);
% saveas(gcf, 'Results/Get_ROI/sim/original.png');


spatialSigma = 0.1;  % Controls the spatial spread
intensitySigma = 0.1;  % Controls intensity similarity
filteredImg = imbilatfilt(grayImg, intensitySigma, spatialSigma);

% figure();
% imshow(filteredImg);

filteredImg2 = medfilt2(grayImg, [13 13]);  % Apply a 3x3 median filter


% figure();
% imshow(filteredImg2);

filteredImg5 = exp(-15*(1-grayImg));

figure();
imshow(filteredImg5);

filteredImg6 = medfilt2(filteredImg5, [13 13]);
figure();
imshow(filteredImg6);

filteredImg9 = exp(-3*(1-filteredImg6));
figure();
imshow(filteredImg9);

binaryMask = filteredImg9 > 0.4;

binaryMask = imclose(binaryMask, strel('disk', 3)); % Close small gaps
binaryMask = imopen(binaryMask, strel('disk', 2));  % Remove small noise

figure();
imshow(binaryMask);

filteredImg4 = medfilt2(binaryMask, [15 15]);  % Apply a 3x3 median filter
figure();
imshow(filteredImg4);

filteredImg7 = medfilt2(filteredImg4, [30 30]); 
figure(); 
imshow(filteredImg7);
% Compute vertical gradient using Sobel filter
% sobelFilter = fspecial('sobel'); % Get Sobel filter
% verticalGradient = imfilter(data, sobelFilter', 'replicate');
% Compute vertical gradient manually
vertgrad = zeros(size(grayImg));
vertgrad(2:end,:) = diff(double(grayImg), 1, 1); % First derivative along rows

% Normalize correctly
% vertgrad = (vertgrad - min(vertgrad(:))) / (max(vertgrad(:)) - min(vertgrad(:)));

% Apply thresholding
threshold = 0.04;
binaryMask = vertgrad > threshold;
% Apply Gaussian smoothing to reduce noise
% smoothedImg = imgaussfilt(double(grayImg), 2); % Sigma=2 for smoothing

% Compute Laplacian to find regions with intensity changes
% laplacianImg = del2(smoothedImg);

% Normalize for consistent thresholding
% laplacianImg = (laplacianImg - min(laplacianImg(:))) / (max(laplacianImg(:)) - min(laplacianImg(:)));

% Apply adaptive thresholding to detect regions of interest
% threshold = graythresh(laplacianImg); % Otsu's method
% binaryMask = laplacianImg > threshold;

% Morphological processing to refine detections
binaryMask = imclose(binaryMask, strel('disk', 3)); % Close small gaps
binaryMask = imopen(binaryMask, strel('disk', 2));  % Remove small noise
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
