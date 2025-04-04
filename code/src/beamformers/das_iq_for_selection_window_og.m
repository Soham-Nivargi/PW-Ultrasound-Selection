function [beamform_next, beamform_comp] = das_iq_for_selection_window_og(scan,pht,dataset,pw_indice)

    assert(~isempty(dataset.modulation_frequency)&&dataset.modulation_frequency~=0,'The supplied dataset is not IQ');

    rx_f_number = 1.75;
    rx_aperture = scan.z/rx_f_number;
    rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey50');

    %-- angular apodization -> no apodization
    angular_apodization = ones(scan.pixels,dataset.firings); 

    %-- beamforming loop
    % beamformed_data = zeros(scan.pixels,length(arr));
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;

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


        % figure();
        % imshow(mask{k});
        title(sprintf('Occlusion %2d', k));
    end

    for i=1:9
        win = tukeywin(sum(mask{i}(:,coordinates{i}(1))), 0.25) * tukeywin(sum(mask{i}(coordinates{i}(2),:)), 0.25)';
        w{i} = zeros(size(X));
        w{i}(mask{i}==1) = win;    
        % figure();
        % imshow(w{i});

        mask_all = mask_all | mask{i};
    end

    mask_all = ~mask_all;

    temp_us = us_image('DAS-IQ beamforming');
    temp_us.author = 'Soham Nivargi';
    temp_us.affiliation = 'IIT Bombay';
    temp_us.algorithm = 'Delay-and-Sum (IQ version)';
    temp_us.scan = scan;
    temp_us.transmit_f_number = 0;
    temp_us.receive_f_number = rx_f_number;
    temp_us.transmit_apodization_window = 'none';
    temp_us.receive_apodization_window = 'Tukey 25%';
    
    contrast_init = zeros(1,9);
    

    for f=1:length(pw_indice) 
        contrast_reg = contrast_init;
        reg_image = zeros(rows,cols);
        reg_image_comp = zeros(rows, cols);
        j = 0;

        for pw=pw_indice{f}
            temp_data = zeros(scan.pixels, 1);
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
                % beamformed_data(:,f) = beamformed_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                temp_data(:,1) = temp_data(:,1)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            end
            envelope_temp_data = reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);
    
            temp_us.number_plane_waves = j;
            temp_us.data = abs(envelope_temp_data + reg_image);
            contrast_curr = tools.contrast_score(scan, pht, temp_us, j);
             
            weights = contrast_curr - contrast_reg;
    
            weights = 1 ./ (1 + exp(-6*weights));
            
            window = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);
            for n=1:9
                window = window + weights(n)*w{n};
            end
    
            window = min(window, 1);
            window = (3+window)/4;

            sigma = 20;  % Increased blur strength
            filter_size = 41;  % Larger filter for smooth blending of 40-50 pixel regions
            h = fspecial('gaussian', filter_size, sigma);

            % Apply the filter
            img_blurred = imfilter(window, h, 'replicate');

            % Display original and blurred images
            % figure;
            % subplot(1,2,1), imshow(window), title('Original Image');
            % subplot(1,2,2), imshow(img_blurred), title('Heavily Blurred Image');

            window = img_blurred;

            if j==0
                reg_image = envelope_temp_data;
                reg_image_comp = abs(envelope_temp_data);
            else
                reg_image = reg_image.*(2-window) + envelope_temp_data.*(window);
                reg_image_comp = ((reg_image_comp*j).*(2-window) + abs(envelope_temp_data.*(window)))/(j+1);
            end

            temp_us.number_plane_waves = j+1;
            temp_us.data = abs(reg_image);
            
            contrast_reg =  tools.contrast_score(scan, pht, temp_us, j+1);
            j = j+1;
        end
        beamform_next{f} = reg_image;
        beamform_comp{f} = reg_image_comp;
    end
end
