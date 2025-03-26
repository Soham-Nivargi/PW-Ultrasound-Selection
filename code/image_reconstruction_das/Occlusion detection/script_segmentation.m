%-- Script to be used as an example to manipulate the provided dataset

%-- After choosing the specific configuration through acquisition_type, 
%-- phantom_type and data_type parameters, this script allows reconstructing
%-- images for evaluation for different choices of steered plane waves involved
%-- in the compounding scheme (specified by the pw_indices parameter)

%-- The implemented method (to be used as example) corresponds to the standard Delay
%-- And Sum (DAS) technique with apodization in reception

%-- Authors: Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)
%--          Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)

%-- $Date: 2016/03/01 $  


clear all;
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

K = 2;
%-- Create path to load corresponding files
path_dataset = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_dataset_',data,'.hdf5'];
path_scan = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_scan.hdf5'];
path_pht = ['../../database/',acquisition,'/',phantom,'/',phantom,'_',acqui,'_phantom.hdf5'];


%-- Read the corresponding dataset and the region where to reconstruct the image
dataset = us_dataset();
dataset.read_file(path_dataset);
scan = linear_scan();
scan.read_file(path_scan);

% arr = cell(1, 1);
% for i = 1:75
%     arr{1} = [38 i]; % Each cell contains the index i
%     if i==38
%         continue
%     end
    %-- Indices of plane waves to be used for each reconstruction
    % pw_indices{1} = [38 25 26 27 28 49 50 51 48]; %overall
    % pw_indices{1} = [38 21 23 28 22 52 51 45 49];%regionwise
    % pw_indices{1} = [35,43,28,51,38,21,29,18,46];
pw_indices{1} = [38 34 42 22 54 30 46 26 50]; %uniform
% pw_indices{1} = [38 30 22 14 8 45 52 60 68];
% pw_indices{2} = 45;
% pw_indices{3} = round(linspace(1,dataset.firings,11));
% pw_indices{1} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves
flag_simu = num2str(flag_simu);
flag_display = num2str(flag_display);

path_reconstruted_img = ['../../reconstructed_image/',acquisition,'/',phantom,'/inv_gaussian_image_',phantom,'_',acqui,'_img_from_',data,'.hdf5'];
% disp(path_reconstruted_img);
% path_image = ['../../reconstructed_image/',acquisition,'/',phantom,'/double_image_',phantom,'_',acqui,'_img_from_',data, '.jpg'];
%-- Reconstruct Bmode images for each pw_indices
disp(['Starting image reconstruction from ',acquisition,' for ',phantom,' using ',data,' dataset'])
% switch data_type    
%     case 1
%         % disp(arr{i});
%         image1 = das_iq_windowed(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display);
%         % image2 = das_iq(scan,dataset,pw_indices);         
%     case 2
%         image = das_rf(scan,dataset,pw_indices);
%     otherwise       %-- Do deal with bad values
%         image = das_iq(scan,dataset,pw_indices);       
% end
% disp('Reconstruction Done')
% disp(['Result saved in "',path_reconstruted_img,'"'])

assert(~isempty(dataset.modulation_frequency)&&dataset.modulation_frequency~=0,'The supplied dataset is not IQ');
pht = us_phantom();
pht.read_file(path_pht);

%-- receive apodization: 
%-- dynamically expanding receive aperture with Tukey 25% apodization
rx_f_number = 1.75;
rx_aperture = scan.z/rx_f_number;
rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey50');

%-- angular apodization -> no apodization
angular_apodization = ones(scan.pixels,dataset.firings);

%-- beamforming loop
beamformed_data = zeros(scan.pixels,length(pw_indices));
time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
w0 = 2*pi*dataset.modulation_frequency;
rows = numel(scan.z_axis);
cols = numel(scan.x_axis);

reg_image = zeros([rows cols  1]);
temp_data = zeros(scan.pixels,1);

pw = pw_indices{1}(1);

transmit_delay = scan.z*cos(dataset.angles(pw))+scan.x*sin(dataset.angles(pw));
for nrx=1:dataset.channels
    %-- receive delay
    receive_delay = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
    %-- total delay
    delay = (transmit_delay+receive_delay)/dataset.c0;
    %-- phase shift
    phase_shift = exp(1i.*w0*(delay-2*scan.z/dataset.c0));

    %-- beamformed data
    beamformed_data(:,1) = beamformed_data(:,1)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
    temp_data(:,1) = temp_data(:,1)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
end
envelope_temp_data = reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);

temp_us = us_image('DAS-IQ beamforming');
temp_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
temp_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
temp_us.algorithm = 'Delay-and-Sum (IQ version)';
temp_us.scan = scan;
temp_us.number_plane_waves = 1;
temp_us.data = abs(envelope_temp_data);
temp_us.transmit_f_number = 0;
temp_us.receive_f_number = rx_f_number;
temp_us.transmit_apodization_window = 'none';
temp_us.receive_apodization_window = 'Tukey 50%';
% dynamic_range = 60;
% temp_us.show(dynamic_range);
% disp(gca);

% Correct definition of centres
centres = {[74,181], [74, 345], [74,507], [196,181], [196, 345], [196,507], [318,181], [318, 345], [318,507]};

% pixel_to_mm_x = 38.0466 / 387; % mm per pixel in x-direction
% pixel_to_mm_y = 44.9462 / 609; % mm per pixel in y-direction


figure();
imshow(envelope_temp_data);
hold on;
for i = 1:length(centres)
    plot(centres{i}(1), centres{i}(2), 'r*', 'MarkerSize', 10);
end
hold off;

% saveas(gcf, 'your_image.jpg');

% image = imread('your_image.jpg'); % Load the image
% figure();
% imshow(image);
% % imshow(image);
% if size(image, 3) > 1
%     image = rgb2gray(image); % Convert to grayscale if necessary
% end
% 
% bandwidth = 100; % Adjust based on target size
% max_iter = 100;
% tol = 1e-3;
% 
% % Mean Shift Algorithm to find target centers
% % Inputs:
% % - image: 2D grayscale image (e.g., ultrasound image)
% % - bandwidth: Bandwidth of the kernel (size of the search window)
% % - max_iter: Maximum number of iterations for convergence
% % - tol: Convergence tolerance (stopping criterion)
% % Outputs:
% % - centers: Nx2 array of target centers (x, y)
% 
% % Normalize the image for consistent intensity scale
% image = double(image);
% image = image / max(image(:));
% 
% % Find high-intensity points as initial seeds
% threshold = 0.5; % Adjust this based on your target intensity
% [y, x] = find(image > threshold);
% seeds = [x, y];
% 
% % Initialize variables
% centers = [];
% used = false(size(seeds, 1), 1); % Tracks if a seed is used
% 
% % Mean Shift Algorithm
% for i = 1:size(seeds, 1)
%     if used(i)
%         continue; % Skip already processed seeds
%     end
% 
%     % Initialize the current point
%     current_point = seeds(i, :);
%     for iter = 1:max_iter
%         % Extract a local region (search window)
%         % Ensure indices are integers
%         x_min = max(floor(current_point(1) - bandwidth), 1);
%         x_max = min(ceil(current_point(1) + bandwidth), size(image, 2));
%         y_min = max(floor(current_point(2) - bandwidth), 1);
%         y_max = min(ceil(current_point(2) + bandwidth), size(image, 1));
% 
%         % Compute the kernel weights (Gaussian)
%         [X, Y] = meshgrid(x_min:x_max, y_min:y_max);
%         distances = (X - current_point(1)).^2 + (Y - current_point(2)).^2;
%         kernel = exp(-distances / (2 * bandwidth^2));
% 
%         % Weighted mean shift
%         weights = kernel .* image(y_min:y_max, x_min:x_max);
%         new_x = sum(sum(X .* weights)) / sum(weights(:));
%         new_y = sum(sum(Y .* weights)) / sum(weights(:));
% 
%         % Check for convergence
%         shift_distance = sqrt((new_x - current_point(1))^2 + (new_y - current_point(2))^2);
%         if shift_distance < tol
%             break;
%         end
% 
%         % Update the current point
%         current_point = [new_x, new_y];
%     end
% 
%     % Add the converged point to centers
%     centers = [centers; current_point];
% 
%     % Mark nearby seeds as used
%     for j = 1:size(seeds, 1)
%         if ~used(j) && sqrt(sum((seeds(j, :) - current_point).^2)) < bandwidth
%             used(j) = true;
%         end
%     end
% end
% 
% 
% % Display the image and centers
% imshow(image, []);
% hold on;
% plot(centers(:, 1), centers(:, 2), 'r*', 'MarkerSize', 10);
% hold off;

% Create a Tukey window
% Create a Tukey window
win = tukeywin(floor(rows/3), 0.25) * tukeywin(floor(cols/3)-1, 0.25)';
figure();
imshow(win, []); % Display the Tukey window

% Initialize the window
window = ones([numel(scan.z_axis), numel(scan.x_axis)]);

% Iterate over centres
for i = 1:length(centres)
    % Calculate the center and bounds for the Tukey window
    center_y = centres{i}(1);
    center_x = centres{i}(2);
    

    half_rows = floor(size(win, 1) / 2);
    half_cols = floor(size(win, 2) / 2);

    x_min = max(1, center_x - half_rows);
    x_max = min(size(window, 1), center_x + half_rows);
    y_min = max(1, center_y - half_cols)+1;
    y_max = min(size(window, 2), center_y + half_cols);

    % Calculate the corresponding indices in the Tukey window
    win_x_min = 1 + max(0, half_rows - center_x + 1);
    win_x_max = size(win, 1) - max(0, center_x + half_rows - size(window, 1));
    win_y_min = 1 + max(0, half_cols - center_y + 1);
    win_y_max = size(win, 2) - max(0, center_y + half_cols - size(window, 2));

    % Assign the Tukey window to the main window
    window(x_min:x_max, y_min:y_max) = win(win_x_min:win_x_max, win_y_min:win_y_max);
end

% Display the resulting window
figure();
imshow(window, []);


% padded_win = padarray(win, [rows-floor(rows/7) cols-floor(cols/6)], 0);
% coordinates = {[74,181], [74, 345], [74,507], [196,181], [196, 345], [196,507], [318,181], [318, 345], [318,507]};
% w = cell(1,9);
% for i=1:9
%     w{i} = padded_win(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1);
%     w{i} = w{i} / max(w{i}(:));
% end




win_top = tukeywin(265+15, 0.25) * tukeywin(132+10, 0.25)';
padded_win_top = padarray(win_top, [rows-floor(280/2) cols-floor(142/2)], 0);

win_other = tukeywin(172+21, 0.25) * tukeywin(132+10, 0.25)';
padded_win_other = padarray(win_other, [rows-floor(193/2) cols-floor(142/2)], 0);
% figure();
% imshow(padded_win_other);
% coordinates = {, [188, 132], , , [196,345], , , , };
coordinates = {[192,132], [192,343], [192,515], [68,132], [68,343], [68,515], [316, 132], [316,343], [316,515]};
w = cell(1,9);
for i=1:9
    if mod(i, 3) == 1
        w{i} = padded_win_top(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1);
    else
        w{i} = padded_win_other(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1); 
    end
    w{i} = w{i} / max(w{i}(:));
    figure();
    imshow(w{i});
end

% win = tukeywin(floor(rows/3.5), 0.25) * tukeywin(floor(cols/3)-1, 0.25)';
% padded_win = padarray(win, [rows-floor(rows/7) cols-floor(cols/6)], 0);
% coordinates = {[74,181], [74, 345], [74,507], [196,181], [196, 345], [196,507], [318,181], [318, 345], [318,507]};
% w = cell(1,9);
% for i=1:9
%     w{i} = padded_win(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1);
%     w{i} = w{i} / max(w{i}(:));
% end


wb = waitbar(0,'DAS beamforming');

for f=1:length(pw_indices)
    reg_image = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);
    j=1;
    % disp(size(im1));
    waitbar(f/length(pw_indices),wb,sprintf('DAS-IQ beamforming %0.0f%%',f/length(pw_indices)*100));
    contrast_reg = zeros([1 9]);
    % disp(contrast_reg);
    for pw=pw_indices{f}
        temp_data = zeros(scan.pixels,length(pw_indices));
        window_data = zeros(scan.pixels,length(pw_indices));
        %-- transmit delay
        transmit_delay = scan.z*cos(dataset.angles(pw))+scan.x*sin(dataset.angles(pw));
        for nrx=1:dataset.channels
            %-- receive delay
            receive_delay = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
            %-- total delay
            delay = (transmit_delay+receive_delay)/dataset.c0;
            %-- phase shift
            phase_shift = exp(1i.*w0*(delay-2*scan.z/dataset.c0));

            %-- beamformed data
            beamformed_data(:,f) = beamformed_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            temp_data(:,f) = temp_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            window_data(:,f) = window_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx);%.*interp1(time_vector,zeros(:,nrx,pw),delay,'spline',0);
        end
        envelope_temp_data = reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);
        % figure();
        % imshow(envelope_temp_data, []);
        % title('OG image');

        temp_us = us_image('DAS-IQ beamforming');
        temp_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
        temp_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
        temp_us.algorithm = 'Delay-and-Sum (IQ version)';
        temp_us.scan = scan;
        temp_us.number_plane_waves = 1;
        temp_us.data = abs(envelope_temp_data);
        temp_us.transmit_f_number = 0;
        temp_us.receive_f_number = rx_f_number;
        temp_us.transmit_apodization_window = 'none';
        temp_us.receive_apodization_window = 'Tukey 50%';

        contrast_curr = tools.contrast_score(scan, pht, temp_us, flag_simu, flag_display, 1);
        disp(contrast_curr);
        disp(contrast_reg);
        
        weights = contrast_curr - contrast_reg;

        % weights = contrast_curr;
        disp(weights);
        
        % if min(weights)<0
        %     weights = weights - min(weights);
        % end
        % weights = weights/max(weights);

        weights = 1 ./ (1 + exp(-0.5*weights));
        disp(weights);
        
        window = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);
        for n=1:9
            window = window + weights(n)*w{n};
        end

        window = window/max(window(:));
        window = (3+window)/4;
        figure();
        imshow(window);
        saveas(gcf, ['Results/uniform_simulation/new_window/window_', num2str(j), '.jpg']);

        % reg_image = reg_image.*(1-window) + envelope_temp_data.*(window);
        if j==1
            reg_image = envelope_temp_data;
        else
            reg_image = reg_image.*(2-window) + envelope_temp_data.*(window);
            % reg_image = reg_image.*((4-window)/5) + envelope_temp_data.*((2+window)/5);
            % reg_image = reg_image.*((2-window)/3) + envelope_temp_data.*((1+window)/3);                
        end
   
        dynamic_range = 60;
        filename1 = ['Results/uniform_simulation/new_window/iteration', num2str(j), '.jpg'];

        reg_us = us_image('DAS-IQ beamforming');
        reg_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
        reg_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
        reg_us.algorithm = 'Delay-and-Sum (IQ version)';
        reg_us.scan = scan;
        reg_us.number_plane_waves = j;
        reg_us.data = abs(reg_image);
        reg_us.transmit_f_number = 0;
        reg_us.receive_f_number = rx_f_number;
        reg_us.transmit_apodization_window = 'none';
        reg_us.receive_apodization_window = 'Tukey 50%';
        reg_us.show(dynamic_range);
        saveas(gcf, filename1);

        contrast_reg = tools.contrast_score(scan, pht, reg_us, flag_simu, flag_display, j);

        disp(contrast_reg);
        % if contrast reg <0 define another window where 

        j = j+1; 
        clc;
        disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])   
    end
end

dynamic_range = 60;
envelope_final_temp = abs(reg_image);
windowed_path_jpg = ['Results/uniform_simulation/new_window/windowed', '.jpg'];
close(wb);
reg_us = us_image('DAS-IQ beamforming');
reg_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
reg_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
reg_us.algorithm = 'Delay-and-Sum (IQ version)';
reg_us.scan = scan;
reg_us.number_plane_waves = 9;
reg_us.data = envelope_final_temp;
% image.data = filtered_img;
reg_us.transmit_f_number = 0;
reg_us.receive_f_number = rx_f_number;
reg_us.transmit_apodization_window = 'none';
reg_us.receive_apodization_window = 'Tukey 50%';
reg_us.show(dynamic_range);
saveas(gcf, windowed_path_jpg);

windowed_path = ['Results/uniform_simulation/new_window/windowed', '.hdf5'];
reg_us.write_file(windowed_path);
%-- reshape
envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  length(pw_indices)]));

unwindowed_path_jpg = ['Results/uniform_simulation/new_window/unwindowed', '.jpg'];
%-- declare an us_image object to store the beamformed data
image = us_image('DAS-IQ beamforming');
image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
image.algorithm = 'Delay-and-Sum (IQ version)';
image.scan = scan;
image.number_plane_waves = cellfun('length',pw_indices);
image.data = envelope_beamformed_data;
% image.data = filtered_img;
image.transmit_f_number = 0;
image.receive_f_number = rx_f_number;
image.transmit_apodization_window = 'none';
image.receive_apodization_window = 'Tukey 50%';

image.show(dynamic_range);
saveas(gcf, unwindowed_path_jpg);

unwindowed_path = ['Results/uniform_simulation/new_window/unwindowed', '.hdf5'];
image.write_file(unwindowed_path);

% contrast_unwindowed = tools.contrast_score(path_scan, path_pht, image, flag_simu, flag_display, 9);
% 
% contrast_windowed = tools.contrast_score(path_scan, path_pht, reg_us, flag_simu, flag_display, 9);
% disp(contrast_windowed);
% disp(contrast_unwindowed);
% disp(contrast_windowed-contrast_unwindowed);

 