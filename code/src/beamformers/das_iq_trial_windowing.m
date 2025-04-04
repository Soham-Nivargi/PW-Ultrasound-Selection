function reg_us = das_iq_trial_windowing(scan,pht,dataset,pw_indices)

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


    %-- beamforming loop
    beamformed_data = zeros(scan.pixels,length(pw_indices));
    rows = numel(scan.z_axis);
    cols = numel(scan.x_axis);
    %----------------------------------  
    %Window 21, Window 22
    % win_top = tukeywin(265+15, 0.25) * tukeywin(132+10, 0.25)';
    % padded_win_top = padarray(win_top, [rows-floor(280/2) cols-floor(142/2)], 0);
    % win_other = tukeywin(172+21, 0.25) * tukeywin(132+10, 0.25)';
    % padded_win_other = padarray(win_other, [rows-floor(193/2) cols-floor(142/2)], 0);
    % coordinates = {[192,132], [192,343], [192,515], [68,132], [68,343], [68,515], [316, 132], [316,343], [316,515]};
    % w = cell(1,9);
    % for i=1:9
    %     if mod(i, 3) == 1
    %         w{i} = padded_win_top(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1);
    %     else
    %         w{i} = padded_win_other(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1); 
    %     end
    %     w{i} = w{i} / max(w{i}(:));
    %     % figure();
    %     % imshow(w{i});
    % end
    
    %------------------------------------------
    % Window 11, Window 12
    % [X, Y] = meshgrid(1:cols, 1:rows);
    % coordinates = {[196,181], [196, 345], [196,507], [74,181], [74, 345], [74,507], [318,181], [318, 345], [318,507]};    
    % win = tukeywin(floor(rows/3), 0.25) * tukeywin(floor(cols/3), 0.25)';
    % w = cell(1,9);
    % mask = cell(1,9);
    % mask_all = zeros(size(X));
    % for i=1:9
    %     mask{i} = abs(X-coordinates{i}(1))<=floor(cols/6) & abs(Y-coordinates{i}(2))<=floor(rows/6);
    %     % w{i} = padded_win(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1);
    %     w{i} = zeros(size(X));
    %     if size(win, 1) == sum(mask{i}(:,coordinates{i}(1)))
    %         if size(win, 2) == sum(mask{i}(coordinates{i}(2),:))
    %             w{i}(mask{i}==1) = win;
    %         else
    %             if coordinates{i}(2) >=  size(X, 2)/2
    %                 win1 = win(:, 1:end-(size(win, 2)-sum(mask{i}(coordinates{i}(2),:))));
    %             else
    %                 win1 = win(:, size(win, 2)-sum(mask{i}(coordinates{i}(2),:)):end);
    %             end
    %             w{i}(mask{i}==1) = win1;
    %         end
    %     else
    %         if size(win, 2) == sum(mask{i}(coordinates{i}(2),:))
    %             if coordinates{i}(2) >=  size(X, 1)/2
    %                 win1 = win(1:end-(size(win, 1)-sum(mask{i}(:,coordinates{i}(1)))), :);
    %             else
    %                 win1 = win((size(win, 1)-sum(mask{i}(:,coordinates{i}(1)))):end, :);
    %             end
    %             w{i}(mask{i}==1) = win1;
    %         else
    %             if coordinates{i}(1) >=  size(X, 1)/2
    %                 win1 = win(1:end-(size(win, 1)-sum(mask{i}(coordinates{i}(1),:))), :);
    %             else
    %                 win1 = win(size(win, 1)-sum(mask{i}(coordinates{i}(1),:)):end, :);
    %             end
    %             if coordinates{i}(2) >=  size(X, 2)/2
    %                 win2 = win1(:, 1:end-(size(win, 2)-sum(mask{i}(coordinates{i}(2),:))));
    %             else
    %                 win2 = win1(:, size(win, 2)-sum(mask{i}(coordinates{i}(2),:)):end);
    %             end
    %             w{i}(mask{i}==1) = win2;
    %         end
    %     end
    %     figure();
    %     imshow(w{i});
    % 
    %     mask_all = mask_all | mask{i};
    % end
    
    %-------------------------------------
    % Window 31, Window 32, Window 33, Window 34
    % [X, Y] = meshgrid(1:cols, 1:rows);
    % 
    % w = cell(1,length(pht.occlusionDiameter));
    % mask_all = zeros(size(X));
    % x = scan.x_matrix;
    % z = scan.z_matrix;
    % mask = cell(1,length(pht.occlusionDiameter));
    % 
    % for k=1:length(pht.occlusionDiameter)
    % 
    %     r = pht.occlusionDiameter(k) / 2;
    %     xc = pht.occlusionCenterX(k);
    %     zc = pht.occlusionCenterZ(k);
    %     mask{k} = ( (abs(x-xc) <= sqrt(2)*r) & abs(z-zc) <= sqrt(2)*r);
    % 
    % 
    %     figure();
    %     imshow(mask{k});
    %     title(sprintf('Occlusion %2d', k));
    % 
    %     mask_all = mask_all | mask{k};
    % end
    % mask_all = ~mask_all;
    
    %------------------------------------
    %Window 13, Window 14, Window 15, Window 16, Window 17, Window 18
    coordinates = {[196,181], [196, 345], [196,507], [74,181], [74, 345], [74,507], [318,181], [318, 345], [318,507]};
    [X, Y] = meshgrid(1:cols, 1:rows);

    w = cell(1,length(pht.occlusionDiameter));
    mask_all = zeros(size(X));
    x = scan.x_matrix;
    z = scan.z_matrix;
    mask = cell(1,length(pht.occlusionDiameter));

    for k=1:length(pht.occlusionDiameter)

        r = pht.occlusionDiameter(k) / 2;
        xc = pht.occlusionCenterX(k);
        zc = pht.occlusionCenterZ(k);
        mask{k} = ( (abs(x-xc) <= sqrt(2)*r) & abs(z-zc) <= sqrt(2)*r);


        figure();
        imshow(mask{k});
        title(sprintf('Occlusion %2d', k));
    end

    for i=1:9
        win = tukeywin(sum(mask{i}(:,coordinates{i}(1))), 0.25) * tukeywin(sum(mask{i}(coordinates{i}(2),:)), 0.25)';
        w{i} = zeros(size(X));
        w{i}(mask{i}==1) = win;    
        figure();
        imshow(w{i});

        mask_all = mask_all | mask{i};
    end

    mask_all = ~mask_all;

    %--------------------------------------
    %Window 41, Window 42
    % coordinates = {[196,181], [196, 345], [196,507], [74,181], [74, 345], [74,507], [318,181], [318, 345], [318,507]};
    % weight_matrix = tools.generate_weight_matrix(rows, cols, coordinates, 50, 'euclidean');

    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    wb = waitbar(0,'DAS beamforming');
    temp_us = us_image('DAS-IQ beamforming');
    temp_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    temp_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    temp_us.algorithm = 'Delay-and-Sum (IQ version)';
    temp_us.scan = scan;
    temp_us.transmit_f_number = 0;
    temp_us.receive_f_number = rx_f_number;
    temp_us.transmit_apodization_window = 'none';
    temp_us.receive_apodization_window = 'Tukey 50%';

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
            end
            envelope_temp_data = reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);
            
            temp_us.number_plane_waves = 1;
            temp_us.data = abs(envelope_temp_data);
            contrast_curr = tools.contrast_score(scan, pht, temp_us, 1);

            % temp_us.number_plane_waves = j;
            % temp_us.data = abs(envelope_temp_data + reg_image);
            % contrast_curr = tools.contrast_score(scan, pht, temp_us, j);
            
            % disp(contrast_curr);
            % disp(contrast_reg);
            % 
            weights = contrast_curr - contrast_reg;
            weights = 1 ./ (1 + exp(-6*weights));
            window = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);
            for n=1:9
                % window = window + weights(n)*weight_matrix(:,:,n);
                window = window + weights(n)*w{n};
                % window = window + weights(n)*mask{n};
            end
            window = window + mean(weights)*mask_all;
            % window = window + 1*mask_all;

            window = min(window, 1);
            % window = (3+window)/4;

            % Define the Gaussian blur filter
            sigma = 20;  % Adjust the standard deviation for blur intensity
            filter_size = 25;  % Adjust kernel size for smoothness
            h = fspecial('gaussian', filter_size, sigma);
            
            % Apply the filter
            img_blurred = imfilter(window, h, 'replicate');
            
            % Display original and blurred images
            figure;
            subplot(1,2,1), imshow(window), title('Original Image');
            subplot(1,2,2), imshow(img_blurred), title('Blurred Image');

            figure();
            imshow(window);
            % saveas(gcf, ['Results weighting strategy 1/simulation/MEAN-FILTER/Window-41 Uniform/TwoMinus/window_', num2str(j), '.jpg']);

            % reg_image = reg_image.*(1-window) + envelope_temp_data.*(window);
            if j==1
                reg_image = envelope_temp_data;
            else
                % reg_image = reg_image.*(1.5-window) + envelope_temp_data.*(0.5+window);
                reg_image = reg_image.*(2-window) + envelope_temp_data.*(window);             
            end

            dynamic_range = 60;
            filename1 = ['Results weighting strategy 1/simulation/MEAN-FILTER/Window-41 Uniform/TwoMinus/iteration', num2str(j), '.jpg'];

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
            % saveas(gcf, filename1);

            contrast_reg = tools.contrast_score(scan, pht, reg_us, j);
            j = j+1; 
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])   
        end
    end

    dynamic_range = 60;
    envelope_final_temp = abs(reg_image);
    windowed_path_jpg = ['Results weighting strategy 1/simulation/MEAN-FILTER/Window-41 Uniform/TwoMinus/windowed_15', '.jpg'];
    close(wb);
    reg_us = us_image('DAS-IQ beamforming');
    reg_us.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    reg_us.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    reg_us.algorithm = 'Delay-and-Sum (IQ version)';
    reg_us.scan = scan;
    reg_us.number_plane_waves = length(pw_indices{f});
    reg_us.data = envelope_final_temp;
    % image.data = filtered_img;
    reg_us.transmit_f_number = 0;
    reg_us.receive_f_number = rx_f_number;
    reg_us.transmit_apodization_window = 'none';
    reg_us.receive_apodization_window = 'Tukey 50%';
    reg_us.show(dynamic_range);
    saveas(gcf, windowed_path_jpg);
    cnr_value = tools.contrast_score(scan, pht, reg_us, length(pw_indices{f}));
    disp(mean(cnr_value));
    windowed_path = ['Results weighting strategy 1/simulation/MEAN-FILTER/Window-41 Uniform/TwoMinus/windowed_15', '.hdf5'];
    reg_us.write_file(windowed_path);
    %-- reshape
    envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  length(pw_indices)]));

    unwindowed_path_jpg = ['Results weighting strategy 1/simulation/MEAN-FILTER/Window-41 Uniform/TwoMinus/unwindowed_15', '.jpg'];
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
    % image.show(dynamic_range);
    % saveas(gcf, unwindowed_path_jpg);

    unwindowed_path = ['Results weighting strategy 1/simulation/MEAN-FILTER/Window-41 Uniform/TwoMinus/unwindowed_15', '.hdf5'];
    image.write_file(unwindowed_path);

    % contrast_unwindowed = tools.contrast_score(path_scan, path_pht, image, flag_simu, flag_display, 9);
    % 
    % contrast_windowed = tools.contrast_score(path_scan, path_pht, reg_us, flag_simu, flag_display, 9);
    % disp(contrast_windowed);
    % disp(contrast_unwindowed);
    % disp(contrast_windowed-contrast_unwindowed);

end

