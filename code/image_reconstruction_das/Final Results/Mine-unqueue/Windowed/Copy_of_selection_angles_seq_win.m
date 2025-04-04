function [beamform_curr, pw_indices, thresh] = Copy_of_selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh,left,mid,right,beamform_curr)
    
    if thresh<=0
        return;
    end
    left_mid = floor((left+mid)/2);
    right_mid = floor((mid+right)/2);
    % disp(max(abs(beamform_curr(:))));

    if ismember(left_mid, pw_indices) || ismember(right_mid, pw_indices)
        return
    end
    
    pw_indice{3} = [mid,left_mid];
    pw_indice{4} = [mid,right_mid];
    pw_indice{2} = [mid,right_mid,left_mid];
    pw_indice{1} = [mid,left_mid,right_mid];
    %-- receive apodization: 
    %-- dynamically expanding receive aperture with Tukey 25% apodization

    [beamform_res, beamform_comp] = das_iq_for_selection_window_og(scan,pht,dataset,pw_indice);
    
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
    end
    if(sum(contrast_reg{1}-contrast_reg{2})>0)
        k=1;
        % beamform_curr = beamform_res{1};
    else
        k=2;
        % beamform_curr = beamform_res{2};
    end

    left_sim = ssim(beamform_comp{k}, beamform_comp{3});
    right_sim = ssim(beamform_comp{k}, beamform_comp{4});
    pw_indices = [pw_indices, pw_indice{k}(2:end)];

    
            if left_sim>right_sim
                % queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                % queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                [beamform_curr, pw_indices, thresh] = Copy_of_selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh-2,...
                    left, left_mid, mid, beamform_curr);
                [beamform_curr, pw_indices, thresh] = Copy_of_selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh,...
                    mid, right_mid, right, beamform_curr);
            else
                % queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
                % queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
                [beamform_curr, pw_indices, thresh] = Copy_of_selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh-2,...
                    mid, right_mid, right, beamform_curr);
                [beamform_curr, pw_indices, thresh] = Copy_of_selection_angles_seq_win(scan,dataset,pht,pw_indices,thresh,...
                    left, left_mid, mid, beamform_curr);
            end
end