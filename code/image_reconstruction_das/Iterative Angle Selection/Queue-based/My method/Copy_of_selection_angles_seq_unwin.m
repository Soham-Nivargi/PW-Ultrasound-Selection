function [beamform_curr, pw_indices, queue] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,queue,beamform_curr)
    % image = 0;
    left_mid = floor((queue{1}(1)+queue{1}(2))/2);
    right_mid = floor((queue{1}(2)+queue{1}(3))/2);
    
    if ismember(left_mid, pw_indices{1}) || ismember(right_mid, pw_indices{1})
        return
    end

    pw_indices{1} = [pw_indices{1},left_mid,right_mid];
    
    beamform_left = das_iq_one_image(scan,dataset,left_mid);
    beamform_right = das_iq_one_image(scan,dataset,right_mid);
    beamform_mid = das_iq_one_image(scan,dataset,queue{1}(2));

    beamform_next{1} = (abs(beamform_mid)+abs(beamform_left))/2;

    beamform_next{2} = (abs(beamform_mid)+abs(beamform_right))/2;
    
    beamform_comp = (abs(beamform_left)+abs(beamform_mid)+abs(beamform_right))/3;

    beamform_curr = beamform_curr + beamform_left + beamform_right;

    left_sim = ssim(beamform_comp, beamform_next{1});
    right_sim = ssim(beamform_comp, beamform_next{2});

    if left_sim>right_sim
        queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
        queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
    else
        queue{end+1} = [queue{1}(2),right_mid,queue{1}(3)];
        queue{end+1} = [queue{1}(1),left_mid,queue{1}(2)];
    end
end