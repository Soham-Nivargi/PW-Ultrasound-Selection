function image = das_iq_windowed(scan,dataset,pw_indices, path_scan, path_pht, flag_simu, flag_display)

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
    rows = numel(scan.z_axis);
    cols = numel(scan.x_axis);

    win = tukeywin(floor(rows/2), 0.25) * tukeywin(floor(cols/2)-1, 0.25)';
    padded_win = padarray(win, [rows-floor(rows/4) cols-floor(cols/4)], 0);
    coordinates = {[74,181], [74, 345], [74,507], [196,181], [196, 345], [196,507], [318,507], [318, 345], [318,181]};
    w = cell(1,9);
    for i=1:9
        w{i} = padded_win(rows-coordinates{i}(2):2*rows-coordinates{i}(2)-1, cols-coordinates{i}(1):2*cols-coordinates{i}(1)-1);
        w{i} = w{i} / max(w{i}(:));
    end
    
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    wb = waitbar(0,'DAS beamforming');

    for f=1:length(pw_indices)
        reg_image = zeros([numel(scan.z_axis) numel(scan.x_axis)  1]);
        j=1;
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

                %-- beamformed data
                beamformed_data(:,f) = beamformed_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                temp_data(:,f) = temp_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
                window_data(:,f) = window_data(:,f)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx);%.*interp1(time_vector,zeros(:,nrx,pw),delay,'spline',0);
            end
            envelope_temp_data = reshape(temp_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);
            figure();
            imshow(envelope_temp_data, []);
            title('OG image');

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

            contrast_temp = tools.contrast_score(path_scan, path_pht, temp_us, flag_simu, flag_display);
            
            [~, idx_u] = sort(contrast_temp, 'descend');
            disp(idx_u(1));
            window = w{idx_u(1)};
            
            window = min(window, 1);
            disp(max(window(:)));

            figure();
            imshow(window);
            title('Window');

            reg_image = reg_image.*((2-window)/3) + envelope_temp_data.*((1+window)/3); %envelope_temp_data .*(w/2)

            filename1 = ['six_regions/trying_us_imag_', num2str(i), '.jpg'];
            dynamic_range = 60;
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

            j = j+1;
            

            
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])   
        end
    end

    dynamic_range = 60;
    envelope_final_temp = abs(reg_image);
    windowed_path_jpg = ['six_regions/trying_window_wala_image', '.jpg'];
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
    
    windowed_path = ['six_regions/trying_window_wala_image', '.hdf5'];
    reg_us.write_file(windowed_path);
    %-- reshape
    envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  length(pw_indices)]));

    unwindowed_path_jpg = ['six_regions/non_window_wala_image', '.jpg'];
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

    unwindowed_path = ['six_regions/non_window_wala_image', '.hdf5'];
    image.write_file(unwindowed_path);

end

