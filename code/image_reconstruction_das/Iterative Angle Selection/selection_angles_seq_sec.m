function [beamform_curr, pw_indices, queue] = selection_angles_seq_sec(scan,dataset,pht,pw_indices,thresh,queue,beamform_curr)
    % image = 0;
    left_mid = floor((queue{1}(1)+queue{1}(2))/2);
    right_mid = floor((queue{1}(2)+queue{1}(3))/2);
    disp(max(abs(beamform_curr(:))));
    
    pw_indices{2} = [pw_indices{1},left_mid];
    pw_indices{3} = [pw_indices{1},right_mid];
    pw_indices{1} = [pw_indices{1},left_mid,right_mid];
    
    beamform_left = das_iq_one_image(scan,dataset,left_mid);
    beamform_right = das_iq_one_image(scan,dataset,right_mid);

    beamform_next{1} = beamform_curr + beamform_left;
    beamform_next{2} = beamform_curr + beamform_right;
    
    if max(abs(beamform_next{1}(:)))>=1
        beamform_l = beamform_next{1}/max(abs(beamform_next{1}(:)));
    else
        beamform_l = beamform_next{1};
    end
    
    if max(abs(beamform_next{2}(:)))>=1
        beamform_r = beamform_next{2}/max(abs(beamform_next{2}(:)));
    else
        beamform_r = beamform_next{2};
    end

    beamform_curr = beamform_curr + beamform_left + beamform_right;

    if max(abs(beamform_curr(:)))>=1
        beamform_curr1 = beamform_curr/max(abs(beamform_curr(:)));
    else
        beamform_curr1 = beamform_curr;
    end

    left_sim = ssim(abs(beamform_curr1), abs(beamform_l));
    right_sim = ssim(abs(beamform_curr1), abs(beamform_r));

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