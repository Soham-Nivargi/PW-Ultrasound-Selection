function envelope_beamformed_data = das_iq_one_image(scan,dataset,pw) %, path_reconstruted_img, path_image

    %-- Function which implements the conventional Delay And Sum (DAS) beamform technique with apodization in reception
    %-- The corresponding code is dedicated to the reconstrucion of dataset (rawdata) saved in IQ format

    %-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

    %-- $Date: 2016/03/01 $  

    assert(~isempty(dataset.modulation_frequency)&&dataset.modulation_frequency~=0,'The supplied dataset is not IQ');

    %-- receive apodization: 
    %-- dynamically expanding receive aperture with Tukey 25% apodization
    rx_f_number = 1.75;
    rx_aperture = scan.z/rx_f_number;
    rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
    receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'tukey50');

    %-- angular apodization -> no apodization
    angular_apodization = ones(scan.pixels,dataset.firings); 

    %-- beamforming loop
    beamformed_data = zeros(scan.pixels,1);
    time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
    w0 = 2*pi*dataset.modulation_frequency;

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
        beamformed_data(:) = beamformed_data(:)+phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
    end
    clc;  
    %-- reshape
    envelope_beamformed_data = reshape(beamformed_data,[numel(scan.z_axis) numel(scan.x_axis)  1]);
end
