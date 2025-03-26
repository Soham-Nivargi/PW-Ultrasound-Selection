function [beamform_curr, pw_indices, queue] = selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh,queue,beamform_curr)
    
    left_mid = floor((queue{1}(1)+queue{1}(2))/2);
    right_mid = floor((queue{1}(2)+queue{1}(3))/2);
    disp(max(abs(beamform_curr(:))));
    
    pw_indice{3} = [left_mid];
    pw_indice{4} = [right_mid];
    pw_indice{2} = [right_mid,left_mid];
    pw_indice{1} = [left_mid,right_mid];
    %-- receive apodization: 
    %-- dynamically expanding receive aperture with Tukey 25% apodization

    beamform_res = das_iq_for_selection_win(scan,pht,dataset,pw_indice,beamform_curr,length(pw_indices));
    
    temp_us = us_image('DAS-IQ beamforming');
    temp_us.author = 'Soham Nivargi';
    temp_us.affiliation = 'IIT Bombay';
    temp_us.algorithm = 'Delay-and-Sum (IQ version)';
    temp_us.scan = scan;
    temp_us.transmit_f_number = 0;
    temp_us.receive_f_number = 1.75;
    temp_us.transmit_apodization_window = 'none';
    temp_us.receive_apodization_window = 'Tukey 25%';

    for i=1:length(beamform_res)
        temp_us.number_plane_waves = length(pw_indices)+length(pw_indice{i});
        temp_us.data = abs(beamform_res{i});

        contrast_reg{i} = tools.contrast_score(scan, pht, temp_us, temp_us.number_plane_waves);

        if max(abs(beamform_res{i}(:)))>=1
           beamform_res1{i} = beamform_res{i}/max(abs(beamform_res{i}(:)));
        else
           beamform_res1{i} = beamform_res{i};
        end
    end
    if(sum(contrast_reg{1}-contrast_reg{2})>0)
        k=1;
        beamform_curr = beamform_res{1};
    else
        k=2;
        beamform_curr = beamform_res{2};
    end
    left_sim = ssim(abs(beamform_res1{k}), abs(beamform_res1{3}));
    right_sim = ssim(abs(beamform_res1{k}), abs(beamform_res1{4}));
    pw_indices = [pw_indices, pw_indice{k}];
    if left_sim < thresh
        if(right_sim<thresh)
            if left_sim>right_sim
                queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
            else
                queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
            end
        else
            queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
        end
    else
        if(right_sim<thresh)
            queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
        end
    end
end