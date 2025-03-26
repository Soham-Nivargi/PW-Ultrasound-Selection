function [beamform_curr, pw_indices, queue] = selection_angles_seq_unwin(scan,dataset,pw_indices,thresh,queue,beamform_curr)
    
    mid = floor(mean(queue{1}));
    
    pw_indices{2} = [queue{1}(1),mid];
    pw_indices{3} = [queue{1}(2),mid];
    pw_indices{4} = [queue{1}(1),queue{1}(2),mid];
    pw_indices{1} = [pw_indices{1}, mid];
    
    beamform_left = das_iq_one_image(scan,dataset,queue{1}(1));
    beamform_right = das_iq_one_image(scan,dataset,queue{1}(2));
    beamform_mid = das_iq_one_image(scan,dataset,mid);

    beamform_next{1} = (abs(beamform_mid)+abs(beamform_left))/2;

    beamform_next{2} = (abs(beamform_mid)+abs(beamform_right))/2;
    
    beamform_comp = (abs(beamform_left)+abs(beamform_mid)+abs(beamform_right))/3;

    left_sim = ssim(beamform_comp, beamform_next{1});
    right_sim = ssim(beamform_comp, beamform_next{2});
    
    beamform_curr = beamform_curr + beamform_mid;

    if left_sim < thresh
        if(right_sim<thresh)
            if left_sim>right_sim
                queue{end+1} = [queue{1}(1),mid];
                queue{end+1} = [mid,queue{1}(2)];
            else
                queue{end+1} = [mid,queue{1}(2)];
                queue{end+1} = [queue{1}(1),mid];
            end
        else
            queue{end+1} = [queue{1}(1),mid];
        end
    else
        if(right_sim<thresh)
            queue{end+1} = [mid,queue{1}(2)];
        end
    end
end