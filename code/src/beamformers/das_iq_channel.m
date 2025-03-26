function das_iq_channel(scan,dataset,pw)
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
end