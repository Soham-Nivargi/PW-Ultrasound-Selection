function image = das_iq(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display, phantom_type)

    %-- Function which implements the conventional Delay And Sum (DAS) beamform technique with apodization in reception
    %-- The corresponding code is dedicated to the reconstrucion of dataset (rawdata) saved in IQ format

    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    assert(~isempty(dataset.modulation_frequency)&&dataset.modulation_frequency~=0,'The supplied dataset is not IQ');

    %-- select the plane waves that will be used in each frame
    if nargin < 3
        pw_indices{1} = 1:dataset.firings;
    end

    %-- receive apodization: 
    %-- dynamically expanding receive aperture with Tukey 25% apodization
    rx_f_number = 1.75;
    rx_aperture = scan.z/rx_f_number;
    rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey50');

    %-- angular apodization -> no apodization
    angular_apodization = ones(scan.pixels,dataset.firings);


    % [rows, ~] = size(scan.x);
    % [cols, ~] = size(scan.z);
    % mid_row = floor(rows / 2);
    % mid_col = floor(cols / 2);

    % Define Gaussian windows with different parameters for each region
    % sigma_tl = 5; weight_tl = 1.5;
    % sigma_bl = 7; weight_bl = 1.2;
    % sigma_tr = 5; weight_tr = 1.3;
    % sigma_br = 7; weight_br = 1.1;

    % Generate Gaussian windows
    % gauss_tl = fspecial('gaussian', [mid_row, mid_col], sigma_tl) * weight_tl;
    % gauss_bl = fspecial('gaussian', [rows - mid_row, mid_col], sigma_bl) * weight_bl;
    % gauss_tr = fspecial('gaussian', [mid_row, cols - mid_col], sigma_tr) * weight_tr;
    % gauss_br = fspecial('gaussian', [rows - mid_row, cols - mid_col], sigma_br) * weight_br;
    

    %-- beamforming loop
    beamformed_data = zeros(scan.pixels,length(pw_indices));
    i = 0;
    j = 1;
    
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    wb = waitbar(0,'DAS beamforming');

    for f=1:length(pw_indices)
        reg_image = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);
        % disp(size(im1));
        waitbar(f/length(pw_indices),wb,sprintf('DAS-IQ beamforming %0.0f%%',f/length(pw_indices)*100));
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

                % Apply the Gaussian window to the plane wave data
                % windowed_data = dataset.data(:, nrx, pw);
                % windowed_data(1:mid_row, 1:mid_col) = windowed_data(1:mid_row, 1:mid_col) .* gauss_tl;
                % windowed_data(mid_row + 1:end, 1:mid_col) = windowed_data(mid_row + 1:end, 1:mid_col) .* gauss_bl;
                % windowed_data(1:mid_row, mid_col + 1:end) = windowed_data(1:mid_row, mid_col + 1:end) .* gauss_tr;
                % windowed_data(mid_row + 1:end, mid_col + 1:end) = windowed_data(mid_row + 1:end, mid_col + 1:end) .* gauss_br;
                % filtered_img = zeros(scan.pixels,length(pw_indices));
                % filtered_img = phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                % sigma = 2; % Standard deviation of the Gaussian filter, adjust as needed
                % filtered_image = imgaussfilt(filtered_img, sigma);
                % Display the original and filtered images
                % og = dataset.data(:,nrx,pw);
                % figure;
                % subplot(1, 2, 1);
                % imshow(og);
                % title('Original Image');

                % subplot(1, 2, 2);
                % imshow(filtered_img);
                % title(['Gaussian Filtered Image (\sigma = ', num2str(sigma), ')']);

                %-- beamformed data
                beamformed_data(:,f) = beamformed_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                temp_data(:,f) = temp_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                window_data(:,f) = window_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx);%.*interp1(time_vector,zeros(:,nrx,pw),delay,'spline',0);
               
                % temp_data(:,f) = phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                % env_temp_data = abs(reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  length(pw_indices)]));
                % % dynamic_range = 60;
                % figure;
                % subplot(1, 2, 1);
                % imshow(env_temp_data);
                % title('Original Image');
            end
            envelope_temp_data = reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);
            % envelope_window_data = abs(reshape(window_data,[numel(scan.z_axis) numel(scan.x_axis)  1]));
            
            rows = numel(scan.z_axis);
            cols = numel(scan.x_axis);
            % figure();
            % imshow(envelope_temp_data, []);
            % title('OG image');
            % for j=1:9
            %     [x, y] = ginput(1);    % Get one point's coordinates from the user
            % 
            %     % Display the clicked coordinates
            %     fprintf('You clicked at X: %.2f, Y: %.2f\n', x, y);
            % end

            
            
            % alpha = 1;
            % if pw==26
            %     cy = 74;  
            %     cx = 181;
            % elseif pw==27
            %     cy = 74;  
            %     cx = 345;
            % elseif pw==38
            %     cy = 74;  
            %     cx = 507;
            % elseif pw==28
            %     cy = 196;
            %     cx = 181;
            % elseif pw==50
            %     cy = 196;
            %     cx = 345;
            % elseif pw==25
            %     cy = 196;
            %     cx = 507;
            % elseif pw==49
            %     cy = 318;
            %     cx = 507;
            % elseif pw==48
            %     cy = 318;
            %     cx = 345;
            % else
            %     cy = 318;
            %     cx = 181;
            % end

            % if pw==21
            %     cy = 74;  
            %     cx = 181;
            % elseif pw==23
            %     cy = 74;  
            %     cx = 345;
            % elseif pw==28
            %     cy = 74;  
            %     cx = 507;
            % elseif pw==38
            %     cy = 196;
            %     cx = 181;
            % elseif pw==52
            %     cy = 196;
            %     cx = 345;
            % elseif pw==22
            %     cy = 196;
            %     cx = 507;
            % elseif pw==51
            %     cy = 318;
            %     cx = 507;
            % elseif pw==45
            %     cy = 318;
            %     cx = 345;
            % else
            %     cy = 318;
            %     cx = 181;
            % end

            if pw < 23
                cy = 74;  
                cx = 181;
            elseif pw<27
                cy = 74;  
                cx = 345;
            elseif pw<31
                cy = 74;  
                cx = 507;
            elseif pw<35
                cy = 196;
                cx = 181;
            elseif pw<40
                cy = 196;
                cx = 345;
            elseif pw<44
                cy = 196;
                cx = 507;
            elseif pw<48
                cy = 318;
                cx = 507;
            elseif pw<52
                cy = 318;
                cx = 345;
            else
                cy = 318;
                cx = 181;
            end
            % 
            % Create 2D Gaussian manually centered in each region
            % [X, Y] = deal(zeros(rows, cols));
            % % Manually calculate the Gaussian weights
            % for i = 1:rows
            %     for j = 1:cols
            %         X(i, j) = (i - cx)^2;  % Distance squared in the x direction
            %         Y(i, j) = (j - cy)^2;  % Distance squared in the y direction
            %     end
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Gaussian
            % Gaussian function: w = exp(-((x - cx)^2 + (y - cy)^2) / (2 * sigma^2))

            % sigma = 50;  % Standard deviation of Gaussian, adjust as needed
            % 
            % 
            % 
            % w = exp(-(X + Y) / (2 * sigma^2));
            % w = w / max(w(:));  % Normalize the window
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tukey
            % % Create a Tukey window along both x and y dimensions
            % tukey_window_x = create_tukey_window(cols, cy, alpha);
            % tukey_window_y = create_tukey_window(rows, cx, alpha);
            w = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);

            % temp = tukeywin(bottom_marg-top_marg+1, 0.25) * tukeywin(right_marg-left_marg+1, 0.25)';
            temp = tukeywin(floor(rows/3), 0.25) * tukeywin(floor(cols/3)-1, 0.25)';

            % figure();
            % imshow(temp);
            % title('temp');

            padded_temp = padarray(temp, [rows-floor(rows/6) cols-floor(cols/6)], 0);
            % figure();
            % imshow(padded_temp);
            % title('temp');
            w = padded_temp(rows-cx:2*rows-cx-1, cols-cy:2*cols-cy-1);

            % figure();
            % imshow(w);
            % title('w');

            % w(top_marg:bottom_marg, left_marg:right_marg) = temp(10:250-1, 15:135-1);
            % disp(size(w));


            % Combine the 1D Tukey windows to form a 2D window
            % w = tukey_window_y * tukey_window_x';  % Outer product
            % 
            % Normalize the Tukey window
            w = w / max(w(:));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % w = w+1;
            % w = w / max(w(:));
            % figure();
            % imshow(w);
            % title('Window');
            % saveas(gcf, ['oct30/window/one/not_averaged_',num2str(j),'.jpg']);
            % 
            figure();
            imshow((1+w)/3);
            title('Window Averaged');
            % saveas(gcf, ['oct30/window/one/averaged_',num2str(j),'.jpg']);
            

            disp(size(w));
            % Blend the region into the final image using the manually computed Gaussian window
            reg_image = reg_image.*((2-w)/3) + envelope_temp_data.*((1+w)/3); %envelope_temp_data .*(w/2)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HANN
            % x = linspace(0, 1, cols);  % Horizontal direction (left to right)
            % y = linspace(0, 1, rows);  % Vertical direction (top to bottom)
            % w1_x = 0.5 * (1 + cos(pi * x));  % Weight for the left side (img1, img3 more)
            % w2_x = 1 - w1_x;  % Weight for the right side (img2, img4 more)
            % 
            % % Vertical weight functions for top-bottom blending
            % w1_y = 0.5 * (1 + cos(pi * y'));  % Weight for the top side (img1, img2 more)
            % w2_y = 1 - w1_y;  % Weight for the bottom side (img3, img4 more)
            % 
            % 
            % disp(pw);
            % if pw<17
            %     w = w2_y * w1_x; % Bottom-left 
            % elseif pw<38
            %     w = w1_y * w1_x; % Top-left
            % elseif pw<55
            %     w = w1_y * w2_x;  % Top-right
            % else
            %     w = w2_y * w2_x;  % Bottom-right
            % 
            % end
            % disp(size(w));
            % w = w/max(w(:));
            % disp(size(envelope_temp_data))
            % 
            % reg_image = reg_image.*(1-w) + envelope_temp_data .* w;
            % % w1 = w1_y * w1_x;  % Top-left (img1)
            % % w2 = w1_y * w2_x;  % Top-right (img2)
            % % w3 = w2_y * w1_x;  % Bottom-left (img3)
            % % w4 = w2_y * w2_x;  % Bottom-right (img4)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filename1 = ['tukey_regionwise_3/half_window_', num2str(i), '.jpg'];
            filename2 = ['oct30/simulation/regioned/one/iteration', num2str(j), '.jpg'];
            % figure();
            % % subplot(1,2,1);
            % % imshow(abs(reg_image), []);
            % % title('Weighted sum image');
            % % subplot(1,2,2);
            % imshow(w, []);
            % title('Window');
            % saveas(gcf, filename1);


            dynamic_range = 60;
            reg_us = us_image('DAS-IQ beamforming');
            reg_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
            reg_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
            reg_us.algorithm = 'Delay-and-Sum (IQ version)';
            reg_us.scan = scan;
            reg_us.number_plane_waves = j;
            reg_us.data = abs(reg_image);
            % image.data = filtered_img;
            reg_us.transmit_f_number = 0;
            reg_us.receive_f_number = rx_f_number;
            reg_us.transmit_apodization_window = 'none';
            reg_us.receive_apodization_window = 'Tukey 50%';
            reg_us.show(dynamic_range);
            saveas(gcf, filename2);

            j = j+1;

            % % disp(rows);
            % % disp(cols);
            % halfRows = floor(rows / 2);
            % halfCols = floor(cols / 2);
            % sigma = 100;
            % 
            % gaussTL = fspecial('gaussian', [halfRows halfCols], sigma);
            % 
            % gaussTL = gaussTL - min(gaussTL(:));  % Shift minimum value to 0
            % gaussTL = gaussTL / max(gaussTL(:));  % Normalize to [0, 1]
            % gaussTL = 0.5 + 0.5 * gaussTL;  % Rescale to [0.5, 1]
            % 
            % % Top-right region Gaussian window (flipped horizontally)
            % % gaussTR = fliplr(gaussTL);
            % % 
            % % % Bottom-left region Gaussian window (flipped vertically)
            % % gaussBL = flipud(gaussTL);
            % % 
            % % % Bottom-right region Gaussian window (flipped both horizontally and vertically)
            % % gaussBR = flipud(fliplr(gaussTL));
            % 
            % % Combine regions into one window
            % gaussian_window = zeros(rows, cols);
            % gaussian_window(1:halfRows, 1:halfCols) = gaussTL;  % Top-left
            % % gaussian_window(1:halfRows, halfCols+1:end) = gaussTR;  % Top-right
            % % gaussian_window(halfRows+1:end, 1:halfCols) = gaussBL;  % Bottom-left
            % % gaussian_window(halfRows+1:end, halfCols+1:end) = gaussBR;  % Bottom-right
            % 
            % 
            % 
            % 
            % figure();
            % subplot(1, 2, 1);
            % imshow(envelope_temp_data, []);
            % title('PW Image');
            % subplot(1,2,2);
            % imshow(gaussian_window, []);
            % title('Window');
            % 
            % windowed_data = envelope_temp_data .* gaussian_window;
            % 
            % figure();
            % subplot(1, 2, 1);
            % imshow(envelope_temp_data, []);
            % title('PW Image');
            % subplot(1,2,2);
            % imshow(windowed_data, []);
            % title('Window image');

            % im1 = us_image('DAS-IQ beamforming');
            % im1.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
            % im1.affiliation = 'Norwegian University of Science and Technology (NTNU)';
            % im1.algorithm = 'Delay-and-Sum (IQ version)';
            % im1.scan = scan;
            % im1.number_plane_waves = 1;
            % im1.data = envelope_temp_data;
            % % image.data = filtered_img;
            % im1.transmit_f_number = 0;
            % im1.receive_f_number = rx_f_number;
            % im1.transmit_apodization_window = 'none';
            % im1.receive_apodization_window = 'Tukey 25%';
            % 
            % im2 = us_image('DAS-IQ beamforming');
            % im2.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
            % im2.affiliation = 'Norwegian University of Science and Technology (NTNU)';
            % im2.algorithm = 'Delay-and-Sum (IQ version)';
            % im2.scan = scan;
            % im2.number_plane_waves = 1;
            % im2.data = windowed_data;
            % % image.data = filtered_img;
            % im2.transmit_f_number = 0;
            % im2.receive_f_number = rx_f_number;
            % im2.transmit_apodization_window = 'none';
            % im2.receive_apodization_window = 'Tukey 25%';

            % 
            % path_reconstructed_img1 = ['../../reconstructed_image/','/no_window','_K_',num2str(K),'.hdf5'];
            % path_reconstructed_img2 = ['../../reconstructed_image/','/yes_window','_K_',num2str(K),'.hdf5'];
            % im1.write_file(path_reconstructed_img1);
            % im2.write_file(path_reconstructed_img2);
            % 
            % path_op1 = ['../../evaluation/','no_window','_K_',num2str(K), '.txt'];
            % % path_op2 = ['../../evaluation/','yes_window','_K_',num2str(K), '.txt'];
            % 
            % switch phantom_type    
            %     case 1 	%-- evaluating resolution and distorsion
            %         tools.exec_evaluation_resolution_distorsion(path_scan,path_pht,path_reconstructed_img1,flag_simu,flag_display,path_op1);
            %     case 2 	%-- evaluating contrast and speckle quality
            %         tools.exec_evaluation_contrast_speckle(path_scan,path_pht,path_reconstructed_img1,flag_simu,flag_display,path_op1);
            %     otherwise       %-- Do deal with bad values
            %         tools.exec_evaluation_resolution(path_scan,path_pht,path_reconstructed_img1,flag_simu,flag_display,path_op1);
            % end
            % disp('Evaluation Done')
            % disp(['Result saved in "',path_op1,'"'])
            % 
            % 
            % % switch phantom_type    
            % %     case 1 	%-- evaluating resolution and distorsion
            % %         tools.exec_evaluation_resolution_distorsion(path_scan,path_pht,path_reconstructed_img2,flag_simu,flag_display,path_op2);
            % %     case 2 	%-- evaluating contrast and speckle quality
            % %         tools.exec_evaluation_contrast_speckle(path_scan,path_pht,path_reconstructed_img2,flag_simu,flag_display,path_op2);
            % %     otherwise       %-- Do deal with bad values
            % %         tools.exec_evaluation_resolution(path_scan,path_pht,path_reconstructed_img2,flag_simu,flag_display,path_op2);
            % % end
            % % disp('Evaluation Done')
            % % disp(['Result saved in "',path_op2,'"'])
            % 
            % K = K+1;
            % disp(K);
            
            % im1.show(dynamic_range);
            % 
            % im2.show(dynamic_range);
            
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])   
        end
    end

    dynamic_range = 60;
    envelope_final_temp = abs(reg_image);
    windowed_path_jpg = ['oct30/simulation/regioned/one/windowed', '.jpg'];
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
    
    windowed_path = ['oct30/simulation/regioned/one/windowed', '.hdf5'];
    reg_us.write_file(windowed_path);
    %-- reshape
    envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  length(pw_indices)]));

    % sigma = 2; % Standard deviation of the Gaussian filter, adjust as needed
    % filtered_img = imgaussfilt(envelope_beamformed_data, sigma);
    
    % Display the original and filtered images
    % figure;
    % subplot(1, 2, 1);
    % imshow(envelope_beamformed_data);
    % title('Original Image');
    % 
    % subplot(1, 2, 2);
    % imshow(filtered_img);
    % title(['Gaussian Filtered Image (\sigma = ', num2str(sigma), ')']);

    unwindowed_path_jpg = ['oct30/simulation/regioned/one/not_windowed', '.jpg'];
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

    unwindowed_path = ['oct30/simulation/regioned/one/not_windowed', '.hdf5'];
    image.write_file(unwindowed_path);

end


