function image = das_iq_original(scan,dataset,arr) %, path_reconstruted_img, path_image

    assert(~isempty(dataset.modulation_frequency)&&dataset.modulation_frequency~=0,'The supplied dataset is not IQ');

    %-- select the plane waves that will be used in each frame
    if nargin < 3
        arr{1} = 1:dataset.firings;
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
    beamformed_data = zeros(scan.pixels,length(arr));
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;
    % wb = waitbar(0,'DAS beamforming');

    for f=1:length(arr) 
        % waitbar(f/length(arr),wb,sprintf('DAS-IQ beamforming %0.0f%%',f/length(arr)*100));
        for pw=arr{f}
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
            end
            clc;
            disp([num2str(pw),' / ',num2str(length(arr{f}))])   
        end
    end
    % close(wb);

    %-- reshape
    envelope_beamformed_data = abs(reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  length(arr)]));

    %-- declare an us_image object to store the beamformed data
    image = us_image('DAS-IQ beamforming');
    image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
    image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
    image.algorithm = 'Delay-and-Sum (IQ version)';
    image.scan = scan;
    image.number_plane_waves = cellfun('length',arr);
    image.data = envelope_beamformed_data;
    image.transmit_f_number = 0;
    image.receive_f_number = rx_f_number;
    image.transmit_apodization_window = 'none';
    image.receive_apodization_window = 'Tukey 25%';
    path_image = ['Results/Paper_sampling_strategy/Reference/all.jpg'];
    dynamic_range = 60;
    image.show(dynamic_range);

    saveas(gcf, path_image);
    path_reconstructed_image = ['Results/Paper_sampling_strategy/Reference/all.hdf5'];
    % -- Save results
    image.write_file(path_reconstructed_image);
end
