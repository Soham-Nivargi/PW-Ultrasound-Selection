function [beamform_curr, pw_indices, queue] = selection_exact_paper_angles(scan,dataset,pht,pw_indices,thresh,queue,beamform_curr)
    
    mid = floor(mean(queue{1}));
    % if()
    pw_indices{5} = setdiff(pw_indices{1}, queue{1});
    pw_indices{2} = [pw_indices{5},queue{1}(1),mid];
    pw_indices{3} = [pw_indices{5},queue{1}(2),mid];
    % pw_indices{4} = [queue{1}(1),queue{1}(2),mid];
    pw_indices{1} = [pw_indices{1}, mid];
    
    beamform_left = das_iq_paper_unwindowed(scan,dataset,pw_indices{2});
    beamform_right = das_iq_paper_unwindowed(scan,dataset,pw_indices{3});
    beamform_curr = das_iq_paper_unwindowed(scan,dataset,pw_indices{1});
    
    % beamform_curr = beamform_curr + das_iq_one_image(scan, dataset, mid);

    if max(abs(beamform_left(:)))>=1
        beamform_l = beamform_left/max(abs(beamform_left(:)));
    else
        beamform_l = beamform_left;
    end
    
    if max(abs(beamform_right(:)))>=1
        beamform_r = beamform_right/max(abs(beamform_right(:)));
    else
        beamform_r = beamform_right;
    end

    if max(abs(beamform_curr(:)))>=1
        beamform_now1 = beamform_curr/max(abs(beamform_curr(:)));
    else
        beamform_now1 = beamform_curr;
    end

    left_sim = ssim(abs(beamform_now1), abs(beamform_l));
    right_sim = ssim(abs(beamform_now1), abs(beamform_r));

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