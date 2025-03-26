function [beamform_curr, pw_indices] = selection_angles_seq(scan,dataset,pht,pw_indices,thresh,left,mid,right,beamform_curr)
    % image = 0;
    left_mid = floor((left+mid)/2);
    right_mid = floor((right+mid)/2);
    disp(max(abs(beamform_curr(:))));
    
    pw_indices{2} = [pw_indices{1},left_mid];
    pw_indices{3} = [pw_indices{1},right_mid];
    pw_indices{1} = [pw_indices{1},left_mid,right_mid];
    
    beamform_left = das_iq_one_image(scan,dataset,left_mid);
    beamform_right = das_iq_one_image(scan,dataset,right_mid);

    beamform_next{1} = (beamform_curr*(length(pw_indices{2})-1) + beamform_left)/length(pw_indices{2});
    beamform_next{2} = (beamform_curr*(length(pw_indices{3})-1) + beamform_right)/length(pw_indices{3});
    
    if max(abs(beamform_next{1}(:)))>=1
        beamform_next{1} = beamform_next{1}/max(abs(beamform_next{1}(:)));
    end
    
    if max(abs(beamform_next{2}(:)))>=1
        beamform_next{2} = beamform_next{2}/max(abs(beamform_next{2}(:)));
    end

    beamform_curr = (beamform_curr*(length(pw_indices{1})-2) + beamform_left + beamform_right)/length(pw_indices{1});
    if max(abs(beamform_curr(:)))>=1
        beamform_curr = beamform_curr/max(abs(beamform_curr(:)));
    end

    left_sim = ssim(abs(beamform_curr), abs(beamform_next{1}));
    right_sim = ssim(abs(beamform_curr), abs(beamform_next{2}));

    if left_sim < thresh
        [beamform_curr, pw_indices] = selection_angles_seq(scan,dataset,pht,pw_indices,thresh, ...
            left,left_mid,mid,beamform_curr);
    end

    if right_sim < thresh
        [beamform_curr, pw_indices] = selection_angles_seq(scan,dataset,pht,pw_indices,thresh, ...
            mid,right_mid,right,beamform_curr);
    end
end