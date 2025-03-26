% Load the mask image
clear;
close all;
clc;
addpath(genpath('../src'));


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

    % figure();
    % imshow(maskOcclusion{k});
    % title(sprintf('Occlusion %2d', k));
    % % filename = fullfile('Results/Count_Weight_Anglewise/', num2str(k), '/','mask_occlusion.jpg');
    % saveas(gcf, ['Results/Count_Weight_Anglewise/', num2str(k), '/','mask_occlusion.jpg'])
    % 
    % figure();
    % imshow(maskInterference{k});
    % title('Occlusion %2d', k);
    % saveas(gcf, ['Results/Count_Weight_Anglewise/', num2str(k), '/','mask_interference.jpg'])
end

% counts = cell(1,9);
% weights = cell(1,9);
angles = linspace(-16, 16, 8);
% count_angles_all = cell(1,9);
% weight_angles_all = cell(1,9);

% count_angles = zeros(size(angles));
% weight_angles = zeros(size(angles));
[height, width] = size(maskOcclusion{i});
% count_all = cell(1,75);
% weight_all = cell(1,75);
x_start = 1:1:width; 
count=zeros(size(angles));

for theta=angles
    k=1;
    rad = deg2rad(theta);
    for j = 1:length(x_start)
        weight = ones(9);
        for inc=1:height
            x_end = x_start(j) + inc * tan(rad);
            y_end = inc;
            dx = x_end - floor(x_end);
            % Bilinear interpolation for maskInterference
            if(x_end>=1 && x_end<=width)
                for i=1:9  
                    Q11 = maskInterference{i}(y_end, floor(x_end));
                    Q21 = maskInterference{i}(y_end, ceil(x_end));
                    interp_maskInterf = Q11 * (1 - dx) + Q21 * dx;
                    weight(i) = weight(i) * (1 - interp_maskInterf * 0.5);
    
                    Q11 = maskOcclusion{i}(y_end, floor(x_end));
                    Q21 = maskOcclusion{i}(y_end, ceil(x_end));
                    interp_maskOccl = Q11 * (1 - dx) + Q21 * dx;

                    count(k) = count(k) + weight(i) * interp_maskOccl;
                end
            end
        end
    end
    k = k+1;
end
% 
%         count_angles(k) = count;
%         % count_all{k} = count;
% 
%         % weight_angles(k) = sum(weight);
%         % weight_all{k} = weight;
%         k = k+1;
%     end
%     count_angles_all{i} = count_angles;
%     % weight_angles_all{i} = weight_angles;
%     % counts{i} = count_all;
%     % weights{i} = weight_all;
% end

temp_count = zeros(size(angles));
temp_weight = zeros(size(angles));
for h=1:9
    temp_count = temp_count+count_angles_all{h};
    temp_weight = temp_weight+weight_angles_all{h};
    
    figure();
    plot(linspace(-16,16,75), count_angles_all{h});
    title('Total Count vs Angle for Occlusion %2d', k);
    saveas(gcf, ['Results/Count_Weight_Anglewise/', num2str(k), '/','count_angle_occlusion.jpg'])

    figure();
    plot(linspace(-16,16,75), count_weights_all{h});
    title('Total Count vs Angle for Occlusion %2d', k);
    saveas(gcf, ['Results/Count_Weight_Anglewise/', num2str(k), '/','cweight_angle_occlusion.jpg'])
end

figure();
plot(linspace(-16,16,75), temp_count);
title('Total Count vs Angle for All Occlusion');
saveas(gcf, 'Results/Count_Weight_Anglewise/Total_count_angle_occlusion.jpg')

figure();
plot(linspace(-16,16,75), temp_count);
title('Total Count vs Angle for All Occlusion');
saveas(gcf, 'Results/Count_Weight_Anglewise/Total_weight_angle_occlusion.jpg')
