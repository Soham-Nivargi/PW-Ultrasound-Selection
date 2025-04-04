function [beamform_max, ind_max] = das_iq_for_greedy_unwin(scan,pht,dataset,pw_indices,iter_beamform)

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
    % beamformed_data = zeros(scan.pixels,length(pw_indices));

    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    wb = waitbar(0,'DAS beamforming'); 
    beamform_max = iter_beamform;
    ind_max = 1;
    mean_max = 0;

    dynamic_range = 60;
    temp_us = us_image('DAS-IQ beamforming');
    temp_us.author = 'Soham Nivargi';
    temp_us.affiliation = 'IIT Bombay';
    temp_us.algorithm = 'Delay-and-Sum (IQ version)';
    temp_us.scan = scan;
    temp_us.transmit_f_number = 0;
    temp_us.receive_f_number = 1.75;
    temp_us.transmit_apodization_window = 'none';
    temp_us.receive_apodization_window = 'Tukey 50%';

    for k=1:length(pw_indices)
        waitbar(k/length(pw_indices),wb,sprintf('DAS-IQ beamforming %0.0f%%',k/length(pw_indices)*100));

        reg_image = iter_beamform;
        disp(['Candidate set: ',num2str(pw_indices{k})]);
        j=1;

        % loop that computes the beamformed data of each image
        for pw=pw_indices{k}
            temp_data = zeros(scan.pixels, length(pw_indices{k}));
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
            reg_image = reg_image + envelope_temp_data;
     
            temp_us.number_plane_waves = j+1;
            temp_us.data = abs(reg_image);

            contrast_reg =  tools.contrast_score(scan, pht, temp_us, j+1);
            mean_reg = mean(contrast_reg);
            j = j+1; 
            clc;
            disp([num2str(pw),' / ',num2str(length(pw_indices{k}))])   
        end
        
        if(mean_reg > mean_max)
            beamform_max = reg_image;
            ind_max = k;
            mean_max = mean_reg;
        end

    end
end

