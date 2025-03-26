function [beamform_curr, pw_indices, queue] = Copy_of_selection_angles_seq_unwin(scan,dataset,pw_indices,queue,beamform_curr)
    % image = 0;
    left_mid = floor((queue{1}(1)+queue{1}(2))/2);
    right_mid = floor((queue{1}(2)+queue{1}(3))/2);
    disp(max(abs(beamform_curr(:))));
    
    pw_indices{2} = [pw_indices{1},left_mid];
    pw_indices{3} = [pw_indices{1},right_mid];
    pw_indices{1} = [pw_indices{1},left_mid,right_mid];
    
    beamform_left = das_iq_one_image(scan,dataset,left_mid);
    beamform_right = das_iq_one_image(scan,dataset,right_mid);

    beamform_next{1} = (abs(beamform_curr)*(length(pw_indices{2})-1) + ...
        abs(beamform_left))/length(pw_indices{2});

    beamform_next{2} = (abs(beamform_curr)*(length(pw_indices{3})-1) + ...
        abs(beamform_right))/length(pw_indices{3});
    
    beamform_comp = (abs(beamform_curr)*(length(pw_indices{1})-2) + abs(beamform_left) + ...
        abs(beamform_right))/length(pw_indices{1});

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